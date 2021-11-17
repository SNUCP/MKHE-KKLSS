package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"

type Evaluator struct {
	params Parameters
	kswRP  *mkrlwe.KeySwitcher
	kswQP  *mkrlwe.KeySwitcher
	conv   *FastBasisExtender
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a small pool of polynomials
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters) *Evaluator {
	eval := new(Evaluator)
	eval.params = params
	eval.kswRP = mkrlwe.NewKeySwitcher(params.paramsRP)
	eval.kswQP = mkrlwe.NewKeySwitcher(params.paramsQP)
	eval.conv = NewFastBasisExtender(params.RingP(), params.RingQ(), params.RingQ1(), params.RingR())

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

func (eval *Evaluator) modUpAndNTT(ctQ, ctR *mkrlwe.Ciphertext) {
	for id := range ctQ.Value {
		eval.conv.ModUpQtoRAndNTT(ctQ.Value[id], ctR.Value[id])
	}
}

func (eval *Evaluator) modUpAndRescaleNTT(ctQ, ctR *mkrlwe.Ciphertext) {
	for id := range ctQ.Value {
		eval.conv.RescaleNTT(ctQ.Value[id], ctR.Value[id])
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
func (eval *Evaluator) MulRelinNew(op0, op1 *Ciphertext, rlkSet *mkrlwe.RelinearizationKeySet) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.mulRelin(op0, op1, rlkSet, ctOut)
	return
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) mulRelin(ct0, ct1 *Ciphertext, rlkSet *mkrlwe.RelinearizationKeySet, ctOut *Ciphertext) {

	ct0R := mkrlwe.NewCiphertextNTT(eval.params.paramsRP, ct0.IDSet(), eval.params.paramsRP.MaxLevel())
	ct1R := mkrlwe.NewCiphertextNTT(eval.params.paramsRP, ct1.IDSet(), eval.params.paramsRP.MaxLevel())
	ctOutR := mkrlwe.NewCiphertextNTT(eval.params.paramsRP, ctOut.IDSet(), eval.params.paramsRP.MaxLevel())

	eval.modUpAndNTT(ct0.Ciphertext, ct0R)
	eval.modUpAndRescaleNTT(ct1.Ciphertext, ct1R)

	eval.kswRP.MulAndRelin(ct0R, ct1R, rlkSet, ctOutR)

	for id := range ctOutR.Value {
		eval.conv.Quantize(ctOutR.Value[id], ctOut.Value[id], eval.params.T())
	}
}
