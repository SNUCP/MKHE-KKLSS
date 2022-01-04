package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"

type Evaluator struct {
	params Parameters
	ksw    *KeySwitcher
	conv   *FastBasisExtender
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a small pool of polynomials
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters) *Evaluator {
	eval := new(Evaluator)
	eval.params = params
	eval.ksw = NewKeySwitcher(params)
	eval.conv = NewFastBasisExtender(params.RingP(), params.RingQ(), params.RingQMul(), params.RingR())

	return eval
}

func (eval *Evaluator) newCiphertextBinary(op0, op1 *Ciphertext) (ctOut *Ciphertext) {
	idset := op0.IDSet().Union(op1.IDSet())
	return NewCiphertext(eval.params, idset)
}

// evaluateInPlaceBinary applies the provided function in place on el0 and el1 and returns the result in elOut.
func (eval *Evaluator) evaluateInPlace(ct0, ct1, ctOut *Ciphertext, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {
	idset0 := ct0.IDSet()
	idset1 := ct1.IDSet()

	evaluate(ct0.Value["0"], ct1.Value["0"], ctOut.Value["0"])
	for id := range ctOut.IDSet().Value {
		if !idset0.Has(id) {
			ctOut.Value[id].Copy(ct1.Value[id])
		} else if !idset1.Has(id) {
			ctOut.Value[id].Copy(ct0.Value[id])
		} else {
			evaluate(ct0.Value[id], ct1.Value[id], ctOut.Value[id])
		}
	}
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *Evaluator) add(op0, op1 *Ciphertext, ctOut *Ciphertext) {
	eval.evaluateInPlace(op0, op1, ctOut, eval.params.RingQ().Add)
}

// AddNew adds op0 to op1 and returns the result in a newly created element.
func (eval *Evaluator) AddNew(op0, op1 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.add(op0, op1, ctOut)
	return
}

// Sub subtracts op1 from op0 and returns the result in ctOut.
func (eval *Evaluator) sub(op0, op1 *Ciphertext, ctOut *Ciphertext) {

	eval.evaluateInPlace(op0, op1, ctOut, eval.params.RingQ().Sub)

	//negate polys which is not contained in op0
	idset0 := op0.IDSet()
	for id := range ctOut.IDSet().Value {
		if !idset0.Has(id) {
			eval.params.RingQ().Neg(ctOut.Value[id], ctOut.Value[id])
		}
	}
}

// SubNew subtracts op1 from op0 and returns the result in a newly created element.
func (eval *Evaluator) SubNew(op0, op1 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.sub(op0, op1, ctOut)

	return
}

// MulRelinNew multiplies ct0 by ct1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelinNew(op0, op1 *Ciphertext, rlkSet *RelinearizationKeySet) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.mulRelin(op0, op1, rlkSet, ctOut)
	return
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) mulRelin(ct0, ct1 *Ciphertext, rlkSet *RelinearizationKeySet, ctOut *Ciphertext) {

	params := eval.params

	ct0NTT := mkrlwe.NewCiphertextNTT(eval.params.Parameters, ct0.IDSet(), eval.params.MaxLevel())
	ct1NTT := mkrlwe.NewCiphertextNTT(eval.params.Parameters, ct1.IDSet(), eval.params.MaxLevel())
	ctOutNTT := mkrlwe.NewCiphertextNTT(eval.params.Parameters, ctOut.IDSet(), eval.params.MaxLevel())

	for id := range ct0NTT.Value {
		params.RingQ().NTT(ct0.Value[id], ct0NTT.Value[id])
	}

	for id := range ct1NTT.Value {
		eval.conv.RescaleNTT(ct1.Value[id], ct1NTT.Value[id])
	}

	eval.ksw.MulAndRelinBFV(ct0NTT, ct1NTT, rlkSet, ctOutNTT)

	for id := range ctOutNTT.Value {
		params.RingQ().InvNTT(ctOutNTT.Value[id], ctOut.Value[id])
	}
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *Evaluator) RotateNew(ct0 *Ciphertext, rotidx int, rkSet *mkrlwe.RotationKeySet) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.IDSet())
	eval.rotate(ct0, rotidx, rkSet, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *Evaluator) rotate(ct0 *Ciphertext, rotidx int, rkSet *mkrlwe.RotationKeySet, ctOut *Ciphertext) {

	// normalize rotidx
	for rotidx >= eval.params.N()/2 {
		rotidx -= eval.params.N() / 2
	}

	for rotidx < 0 {
		rotidx += eval.params.N() / 2
	}

	if rotidx == 0 {
		ctOut.Ciphertext.Copy(ct0.Ciphertext)
		return
	}

	_, in := eval.params.CRS[rotidx]

	ringQ := eval.params.RingQ()
	ctTmp := ct0.CopyNew()

	for id := range ctOut.Value {
		ringQ.NTT(ctTmp.Value[id], ctTmp.Value[id])
		ctTmp.Value[id].IsNTT = true
		ctOut.Value[id].IsNTT = true
	}

	if in {
		eval.ksw.Rotate(ctTmp.Ciphertext, rotidx, rkSet, ctOut.Ciphertext)
	} else {
		for k := 1; rotidx > 0; k *= 2 {
			if rotidx%2 != 0 {
				eval.ksw.Rotate(ctTmp.Ciphertext, k, rkSet, ctOut.Ciphertext)
				ctTmp.Ciphertext.Copy(ctOut.Ciphertext)
			}
			rotidx /= 2
		}
	}

	for id := range ctOut.Value {
		ringQ.InvNTT(ctOut.Value[id], ctOut.Value[id])
		ctOut.Value[id].IsNTT = false
	}

}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *Evaluator) ConjugateNew(ct0 *Ciphertext, ckSet *mkrlwe.ConjugationKeySet) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.IDSet())
	eval.conjugate(ct0, ckSet, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *Evaluator) conjugate(ct0 *Ciphertext, ckSet *mkrlwe.ConjugationKeySet, ctOut *Ciphertext) {

	ringQ := eval.params.RingQ()
	ctTmp := ct0.CopyNew()

	for id := range ctOut.Value {
		ringQ.NTT(ctTmp.Value[id], ctTmp.Value[id])
		ctTmp.Value[id].IsNTT = true
		ctOut.Value[id].IsNTT = true
	}

	eval.ksw.Conjugate(ctTmp.Ciphertext, ckSet, ctOut.Ciphertext)

	for id := range ctOut.Value {
		ringQ.InvNTT(ctOut.Value[id], ctOut.Value[id])
		ctOut.Value[id].IsNTT = false
	}
}
