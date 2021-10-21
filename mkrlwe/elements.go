package mkrlwe

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/utils"

type Ciphertext struct {
	Value0 *ring.Poly
	Value  map[string]*ring.Poly
}

type PolyVector struct {
	Value []ring.Poly
}

// NewCiphertext returns a new Element with zero values
func NewCiphertext(params rlwe.Parameters, idset []string, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value0 = ring.NewPoly(params.N(), level+1)
	el.Value = make(map[string]*ring.Poly)

	for _, id := range idset {
		el.Value[id] = ring.NewPoly(params.N(), level+1)
	}

	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set
func NewCiphertextNTT(params rlwe.Parameters, idset []string, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value0 = ring.NewPoly(params.N(), level+1)
	el.Value0.IsNTT = true

	el.Value = make(map[string]*ring.Poly)

	for _, id := range idset {
		el.Value[id] = ring.NewPoly(params.N(), level+1)
		el.Value[id].IsNTT = true
	}

	return el
}

// SetValue sets the input slice of polynomials as the value of the target element
func (el *Ciphertext) SetValue(value0 *ring.Poly, value map[string]*ring.Poly) {
	el.Value0 = value0
	el.Value = value
}

// NumID returns the number of IDs engaged in the target elements
func (el *Ciphertext) NumID() int {
	return len(el.Value)
}

// Level returns the level of the target elements
func (el *Ciphertext) Level() int {
	return len(el.Value0.Coeffs) - 1
}

// SwitchCiphertextRingDegreeNTT changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (else a nil pointer).
// The ctIn must be in the NTT domain and ctOut will be in the NTT domain.
func SwitchCiphertextRingDegreeNTT(ctIn *Ciphertext, ringQLargeDim *ring.Ring, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value0.Coeffs[0]), len(ctOut.Value0.Coeffs[0])

	if NIn > NOut {
		r := ringQLargeDim
		gap := NIn / NOut
		pool := make([]uint64, NIn)
		for i := range ctOut.Value {
			for j := range ctOut.Value[i].Coeffs {
				tmp0, tmp1 := ctOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
				ring.InvNTT(tmp1, pool, NIn, r.NttPsiInv[j], r.NttNInv[j], r.Modulus[j], r.MredParams[j])
				for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+1, w1+gap {
					pool[w0] = pool[w1]
				}
				ring.NTT(pool, tmp0, NOut, r.NttPsi[j], r.Modulus[j], r.MredParams[j], r.BredParams[j])
			}
		}

		for j := range ctOut.Value0.Coeffs {
			tmp0, tmp1 := ctOut.Value0.Coeffs[j], ctIn.Value0.Coeffs[j]
			ring.InvNTT(tmp1, pool, NIn, r.NttPsiInv[j], r.NttNInv[j], r.Modulus[j], r.MredParams[j])
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+1, w1+gap {
				pool[w0] = pool[w1]
			}
			ring.NTT(pool, tmp0, NOut, r.NttPsi[j], r.Modulus[j], r.MredParams[j], r.BredParams[j])
		}

	} else {
		for i := range ctOut.Value {
			ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value[i], ctOut.Value[i])
		}
		ring.MapSmallDimensionToLargerDimensionNTT(ctIn.Value0, ctOut.Value0)
	}
}

// SwitchCiphertextRingDegree changes the ring degree of ctIn to the one of ctOut.
// Maps Y^{N/n} -> X^{N} or X^{N} -> Y^{N/n}.
// If the ring degree of ctOut is larger than the one of ctIn, then the ringQ of ctIn
// must be provided (else a nil pointer).
func SwitchCiphertextRingDegree(ctIn *Ciphertext, ctOut *Ciphertext) {

	NIn, NOut := len(ctIn.Value0.Coeffs[0]), len(ctOut.Value0.Coeffs[0])

	gapIn, gapOut := NOut/NIn, 1
	if NIn > NOut {
		gapIn, gapOut = 1, NIn/NOut
	}

	for i := range ctOut.Value {
		for j := range ctOut.Value[i].Coeffs {
			tmp0, tmp1 := ctOut.Value[i].Coeffs[j], ctIn.Value[i].Coeffs[j]
			for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
				tmp0[w0] = tmp1[w1]
			}
		}
	}

	for j := range ctOut.Value0.Coeffs {
		tmp0, tmp1 := ctOut.Value0.Coeffs[j], ctIn.Value0.Coeffs[j]
		for w0, w1 := 0, 0; w0 < NOut; w0, w1 = w0+gapIn, w1+gapOut {
			tmp0[w0] = tmp1[w1]
		}
	}
}

// CopyNew creates a new element as a copy of the target element.
func (el *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)

	ctxCopy.Value0 = el.Value0.CopyNew()
	for id := range el.Value {
		ctxCopy.Value[id] = el.Value[id]
	}

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Ciphertext) Copy(ctxCopy *Ciphertext) {

	if el != ctxCopy {
		for id := range ctxCopy.Value {
			el.Value[id].Copy(ctxCopy.Value[id])
		}
		el.Value0.Copy(ctxCopy.Value0)
	}

}

// El returns a pointer to this Element
func (el *Ciphertext) El() *Ciphertext {
	return el
}

// RLWEElement returns a pointer to this Element
func (el *Ciphertext) MKRLWEElement() *Ciphertext {
	return el
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients
func PopulateElementRandom(prng utils.PRNG, params rlwe.Parameters, el *Ciphertext) {

	ringQ, err := ring.NewRing(params.N(), params.Q())
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)

	for id := range el.Value {
		sampler.Read(el.Value[id])
	}
	sampler.Read(el.Value0)

}
