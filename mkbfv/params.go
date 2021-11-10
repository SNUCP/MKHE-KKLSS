package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/bfv"

// Parameters represents a parameter set for the BFV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
// R = Q*QMul
type Parameters struct {
	paramsQP mkrlwe.Parameters
	paramsRP mkrlwe.Parameters
	ringQMul *ring.Ring
	ringR    *ring.Ring
	ringT    *ring.Ring
	CRS      [2]*mkrlwe.SwitchingKey
}

// NewParameters instantiate a set of MKCKKS parameters from the generic CKKS parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(bfvParams bfv.Parameters) (params Parameters) {

	ringP := bfvParams.RingP()
	ringQ := bfvParams.RingQ()
	ringQMul := bfvParams.RingQMul()

	r := make([]uint64, 0)
	r = append(r, ringQ.Modulus...)
	r = append(r, ringQMul.Modulus...)

	ringR, err := ring.NewRing(bfvParams.N(), r)
	if err != nil {
		panic("cannot NewParameters: ring R cannot be generated")
	}

	paramsRPRLWE, err := rlwe.NewParameters(bfvParams.N(), ringR.Modulus, ringP.Modulus, bfvParams.Sigma())
	if err != nil {
		panic("cannot NewParameters: paramsRP cannot be generated")
	}

	params.paramsRP = mkrlwe.NewParameters(paramsRPRLWE, 3)
	params.paramsQP = mkrlwe.NewParameters(bfvParams.Parameters, 3)
	params.ringQMul = bfvParams.RingQMul()
	params.ringR = ringR
	params.ringT = bfvParams.RingT()
	params.CRS[0] = params.paramsRP.CRS[0]
	params.CRS[1] = params.paramsRP.CRS[1]

	return params
}

// RingQMul returns a pointer to the ring of the extended basis for multiplication
func (p Parameters) RingQMul() *ring.Ring {
	return p.ringQMul
}

func (p Parameters) RingQ() *ring.Ring {
	return p.paramsQP.RingQ()
}

func (p Parameters) RingP() *ring.Ring {
	return p.paramsQP.RingP()
}

func (p Parameters) RingQP() *rlwe.RingQP {
	return p.paramsQP.RingQP()
}

func (p Parameters) RingR() *ring.Ring {
	return p.paramsRP.RingQ()
}

func (p Parameters) RingRP() *rlwe.RingQP {
	return p.paramsRP.RingQP()
}

// T returns the plaintext coefficient modulus t
func (p Parameters) T() uint64 {
	return p.ringT.Modulus[0]
}

// RingT returns a pointer to the plaintext ring
func (p Parameters) RingT() *ring.Ring {
	return p.ringT
}

func (p Parameters) QCount() int {
	return p.paramsQP.QCount()
}

func (p Parameters) PCount() int {
	return p.paramsQP.PCount()
}

func (p Parameters) RCount() int {
	return p.paramsRP.QCount()
}

func (p Parameters) QMulCount() int {
	return len(p.RingQMul().Modulus)
}
