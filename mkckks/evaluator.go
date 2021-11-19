package mkckks

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"

import "math"
import "unsafe"
import "errors"

type Evaluator struct {
	params    Parameters
	ksw       *mkrlwe.KeySwitcher
	ctxtPool  *mkrlwe.Ciphertext
	polyQPool *ring.Poly
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a small pool of polynomials
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params Parameters) *Evaluator {
	eval := new(Evaluator)
	eval.params = params

	if params.PCount() != 0 {
		eval.ksw = mkrlwe.NewKeySwitcher(params.Parameters)
	}

	eval.ctxtPool = mkrlwe.NewCiphertext(params.Parameters, mkrlwe.NewIDSet(), params.MaxLevel())

	ringQ := params.RingQ()
	eval.polyQPool = ringQ.NewPoly()
	eval.polyQPool.IsNTT = true

	return eval
}

func (eval *Evaluator) getConstAndScale(level int, constant interface{}) (cReal, cImag, scale float64) {

	// Converts to float64 and determines if a scaling is required (which is the case if either real or imag have a rational part)
	scale = 1
	switch constant := constant.(type) {
	case complex128:
		cReal = real(constant)
		cImag = imag(constant)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

	case float64:
		cReal = constant
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.params.RingQ().Modulus[level])
			}
		}

	case uint64:
		cReal = float64(constant)
		cImag = float64(0)

	case int64:
		cReal = float64(constant)
		cImag = float64(0)

	case int:
		cReal = float64(constant)
		cImag = float64(0)
	}

	return
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *Evaluator) DropLevelNew(ct0 *Ciphertext, levels int) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *Evaluator) DropLevel(ct0 *Ciphertext, levels int) {
	level := ct0.Level()

	for id := range ct0.Value {
		ct0.Value[id].Coeffs = ct0.Value[id].Coeffs[:level+1-levels]
	}
}

// MultByConst multiplies ct0 by the input constant and returns the result in ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *Evaluator) MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinInt(ct0.Level(), ctOut.Level())

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	// Component wise multiplication of the following vector with the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	ringQ := eval.params.RingQ()
	var scaledConst, scaledConstReal, scaledConstImag uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		bredParams := ringQ.BredParams[i]
		mredParams := ringQ.MredParams[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.MRed(scaledConstImag, ringQ.NttPsi[i][1], qi, mredParams)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConst = ring.MForm(scaledConst, qi, bredParams)

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]

			for j := 0; j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.MRed(x[0], scaledConst, qi, mredParams)
				z[1] = ring.MRed(x[1], scaledConst, qi, mredParams)
				z[2] = ring.MRed(x[2], scaledConst, qi, mredParams)
				z[3] = ring.MRed(x[3], scaledConst, qi, mredParams)
				z[4] = ring.MRed(x[4], scaledConst, qi, mredParams)
				z[5] = ring.MRed(x[5], scaledConst, qi, mredParams)
				z[6] = ring.MRed(x[6], scaledConst, qi, mredParams)
				z[7] = ring.MRed(x[7], scaledConst, qi, mredParams)
			}
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConst = ring.MForm(scaledConst, qi, bredParams)
		}

		for u := range ct0.Value {
			p0tmp := ct0.Value[u].Coeffs[i]
			p1tmp := ctOut.Value[u].Coeffs[i]
			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.MRed(x[0], scaledConst, qi, mredParams)
				z[1] = ring.MRed(x[1], scaledConst, qi, mredParams)
				z[2] = ring.MRed(x[2], scaledConst, qi, mredParams)
				z[3] = ring.MRed(x[3], scaledConst, qi, mredParams)
				z[4] = ring.MRed(x[4], scaledConst, qi, mredParams)
				z[5] = ring.MRed(x[5], scaledConst, qi, mredParams)
				z[6] = ring.MRed(x[6], scaledConst, qi, mredParams)
				z[7] = ring.MRed(x[7], scaledConst, qi, mredParams)
			}
		}
	}

	ctOut.Scale = ct0.Scale * scale
}

func (eval *Evaluator) evaluateInPlace(c0, c1, ctOut *Ciphertext, evaluate func(int, *ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *mkrlwe.Ciphertext

	idset0 := c0.IDSet()
	idset1 := c1.IDSet()
	idset := idset0.Union(idset1)

	eval.ctxtPool.PadCiphertext(idset)

	level := utils.MinInt(utils.MinInt(c0.Level(), c1.Level()), ctOut.Level())

	c0Scale := c0.ScalingFactor()
	c1Scale := c1.ScalingFactor()
	ctOutScale := ctOut.ScalingFactor()

	if ctOut.Level() > level {
		eval.DropLevel(&Ciphertext{ctOut.El(), ctOutScale}, ctOut.Level()-utils.MinInt(c0.Level(), c1.Level()))
	}

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.ctxtPool.El()

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{tmp1, ctOutScale})

		} else if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{c0.El(), c0Scale})

			ctOut.SetScalingFactor(c1Scale)

			tmp1 = c1.El()

		} else {

			tmp1 = c1.El()
		}

		tmp0 = c0.El()

	} else if ctOut == c1 {

		if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.ctxtPool.El()

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{tmp0, ctOutScale})

		} else if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{ctOut.El(), ctOutScale})

			ctOut.SetScalingFactor(c0Scale)

			tmp0 = c0.El()

		} else {

			tmp0 = c0.El()
		}

		tmp1 = c1.El()

	} else {

		if c1Scale > c0Scale && math.Floor(c1Scale/c0Scale) > 1 {

			tmp0 = eval.ctxtPool.El()

			eval.MultByConst(&Ciphertext{c0.El(), c0Scale}, math.Floor(c1Scale/c0Scale), &Ciphertext{tmp0, ctOutScale})

			tmp1 = c1.El()

		} else if c0Scale > c1Scale && math.Floor(c0Scale/c1Scale) > 1 {

			tmp1 = eval.ctxtPool.El()

			eval.MultByConst(&Ciphertext{c1.El(), c1Scale}, math.Floor(c0Scale/c1Scale), &Ciphertext{tmp1, ctOutScale})

			tmp0 = c0.El()

		} else {
			tmp0 = c0.El()
			tmp1 = c1.El()
		}
	}

	evaluate(level, tmp0.Value["0"], tmp1.Value["0"], ctOut.Value["0"])
	for id := range ctOut.IDSet().Value {
		if !idset0.Has(id) {
			ring.CopyValuesLvl(level, tmp1.Value[id], ctOut.Value[id])
		} else if !idset1.Has(id) {
			ring.CopyValuesLvl(level, tmp0.Value[id], ctOut.Value[id])
		} else {
			evaluate(level, tmp0.Value[id], tmp1.Value[id], ctOut.Value[id])
		}
	}

}

func (eval *Evaluator) newCiphertextBinary(op0, op1 *Ciphertext) (ctOut *Ciphertext) {

	maxScale := utils.MaxFloat64(op0.ScalingFactor(), op1.ScalingFactor())
	minLevel := utils.MinInt(op0.Level(), op1.Level())
	idset := op0.IDSet().Union(op1.IDSet())

	return NewCiphertext(eval.params, idset, minLevel, maxScale)
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *Evaluator) add(op0, op1 *Ciphertext, ctOut *Ciphertext) {
	eval.evaluateInPlace(op0, op1, ctOut, eval.params.RingQ().AddLvl)

}

// AddNew adds op0 to op1 and returns the result in a newly created element.
func (eval *Evaluator) AddNew(op0, op1 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.add(op0, op1, ctOut)
	return
}

// Sub subtracts op1 from op0 and returns the result in ctOut.
func (eval *Evaluator) sub(op0, op1 *Ciphertext, ctOut *Ciphertext) {

	eval.evaluateInPlace(op0, op1, ctOut, eval.params.RingQ().SubLvl)

	level := utils.MinInt(utils.MinInt(op0.Level(), op1.Level()), ctOut.Level())

	//negate polys which is not contained in op0
	idset0 := op0.IDSet()
	for id := range ctOut.IDSet().Value {
		if !idset0.Has(id) {
			eval.params.RingQ().NegLvl(level, ctOut.Value[id], ctOut.Value[id])
		}
	}

}

// SubNew subtracts op1 from op0 and returns the result in a newly created element.
func (eval *Evaluator) SubNew(op0, op1 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.sub(op0, op1, ctOut)

	return
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "minScale <= 0", ct.Scale = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != ctOut.Level()
func (eval *Evaluator) Rescale(ctIn *Ciphertext, minScale float64, ctOut *Ciphertext) (err error) {

	ringQ := eval.params.RingQ()

	if minScale <= 0 {
		return errors.New("cannot Rescale: minScale is 0")
	}

	if ctIn.Scale == 0 {
		return errors.New("cannot Rescale: ciphertext scale is 0")
	}

	if ctIn.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	ctOut.Scale = ctIn.Scale

	var nbRescales int
	// Divides the scale by each moduli of the modulus chain as long as the scale isn't smaller than minScale/2
	// or until the output Level() would be zero
	for ctOut.Scale/float64(ringQ.Modulus[ctIn.Level()-nbRescales]) >= minScale/2 && ctIn.Level()-nbRescales >= 0 {
		ctOut.Scale /= (float64(ringQ.Modulus[ctIn.Level()-nbRescales]))
		nbRescales++
	}

	if nbRescales > 0 {
		level := ctIn.Level()
		for i := range ctOut.Value {
			ringQ.DivRoundByLastModulusManyNTTLvl(level, nbRescales, ctIn.Value[i], eval.polyQPool, ctOut.Value[i])
			ctOut.Value[i].Coeffs = ctOut.Value[i].Coeffs[:level+1-nbRescales]
		}
	} else {
		if ctIn != ctOut {
			ctOut = ctIn.CopyNew()
		}
	}

	return nil
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.Scale = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval *Evaluator) RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error) {

	ctOut = NewCiphertext(eval.params, ct0.IDSet(), ct0.Level(), ct0.Scale)

	return ctOut, eval.Rescale(ct0, threshold, ctOut)
}

// MulRelinNew multiplies ct0 by ct1 with relinearization and returns the result in a newly created element.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) MulRelinNew(op0, op1 *Ciphertext, rlkSet *mkrlwe.RelinearizationKeySet) (ctOut *Ciphertext) {
	//ctOut = NewCiphertext(eval.params, op0.IDSet().Union(op1.IDSet()), utils.MinInt(op0.Level(), op1.Level()), 0)
	ctOut = eval.newCiphertextBinary(op0, op1)
	ctOut.Scale = 0
	eval.mulRelin(op0, op1, rlkSet, ctOut)
	return
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// The procedure will panic if either op0.Degree or op1.Degree > 1.
// The procedure will panic if ctOut.Degree != op0.Degree + op1.Degree.
// The procedure will panic if the evaluator was not created with an relinearization key.
func (eval *Evaluator) mulRelin(op0, op1 *Ciphertext, rlkSet *mkrlwe.RelinearizationKeySet, ctOut *Ciphertext) {

	level := utils.MinInt(utils.MinInt(op0.Level(), op1.Level()), ctOut.Level())

	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	ctOut.Scale = op0.ScalingFactor() * op1.ScalingFactor()
	eval.ksw.MulAndRelin(op0.Ciphertext, op1.Ciphertext, rlkSet, ctOut.Ciphertext)
	eval.Rescale(ctOut, eval.params.Scale(), ctOut)
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *Evaluator) RotateNew(ct0 *Ciphertext, rotidx int, rkSet *mkrlwe.RotationKeySet) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.IDSet(), ct0.Level(), ct0.Scale)
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

	if in {
		eval.ksw.Rotate(ct0.Ciphertext, rotidx, rkSet, ctOut.Ciphertext)
		return
	}

	ctTmp := ct0.CopyNew()
	for k := 1; rotidx > 0; k *= 2 {
		if rotidx%2 != 0 {
			eval.ksw.Rotate(ctTmp.Ciphertext, k, rkSet, ctOut.Ciphertext)
			ctTmp.Ciphertext.Copy(ctOut.Ciphertext)
		}
		rotidx /= 2
	}

}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *Evaluator) ConjugateNew(ct0 *Ciphertext, ckSet *mkrlwe.ConjugationKeySet) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.IDSet(), ct0.Level(), ct0.Scale)
	eval.conjugate(ct0, ckSet, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *Evaluator) conjugate(ct0 *Ciphertext, ckSet *mkrlwe.ConjugationKeySet, ctOut *Ciphertext) {
	eval.ksw.Conjugate(ct0.Ciphertext, ckSet, ctOut.Ciphertext)
}
