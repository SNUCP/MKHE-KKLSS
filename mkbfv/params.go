package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
type ParametersLiteral struct {
	LogN    int // Log Ring degree (power of 2)
	Q       []uint64
	QMul    []uint64
	P       []uint64
	LogQ    []int   `json:",omitempty"`
	LogQMul []int   `json:",omitempty"`
	LogP    []int   `json:",omitempty"`
	Sigma   float64 // Gaussian sampling standard deviation
	T       uint64  // Plaintext modulus
}

// Parameters represents a parameter set for the BFV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
// R = Q*QMul
type Parameters struct {
	paramsQP    mkrlwe.Parameters
	paramsQMulP mkrlwe.Parameters
	paramsRP    mkrlwe.Parameters
	ringT       *ring.Ring
}

// NewParameters instantiate a set of MKCKKS parameters from the generic CKKS parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(pl ParametersLiteral) (params Parameters) {

	if len(pl.Q) != len(pl.QMul) {
		panic("cannot NewParametersFromLiteral: length of Q & QMul is not equal")
	}

	N := (1 << pl.LogN)
	R := make([]uint64, 0)
	R = append(R, pl.Q...)
	R = append(R, pl.QMul...)

	ringT, err := ring.NewRing(N, []uint64{pl.T})
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring T cannot be generated")
	}

	params.ringT = ringT

	rlweParamsQP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring QP cannot be generated")
	}

	rlweParamsQMulP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.QMul, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring QMulP cannot be generated")
	}

	rlweParamsRP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: R, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic(err)
		panic("cannot NewParametersFromLiteral: ring RP cannot be generated")

	}

	params.paramsQP = mkrlwe.NewParameters(rlweParamsQP, 2)
	params.paramsQMulP = mkrlwe.NewParameters(rlweParamsQMulP, 2)
	params.paramsRP = mkrlwe.NewParameters(rlweParamsRP, 2)

	conv := NewFastBasisExtender(
		params.paramsQP.RingP(), params.paramsQP.RingQ(),
		params.paramsQMulP.RingQ(),
		params.paramsRP.RingQ(),
	)

	idxs := []int{
		0, -1, //CRS for relin key
		-2, //CRS for conj key
	}

	// CRS for rot keys
	for i := 0; i < params.paramsQP.LogN()-1; i++ {
		idxs = append(idxs, 1<<i)
	}

	// apply GadgetTransform
	for _, idx := range idxs {
		conv.GadgetTransform(params.paramsQP.CRS[idx], params.paramsQMulP.CRS[idx], params.paramsRP.CRS[idx])
	}
	return params
}

func (p Parameters) RingQMul() *ring.Ring {
	return p.paramsQMulP.RingQ()
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

func (p Parameters) RingQMulP() *rlwe.RingQP {
	return p.paramsQMulP.RingQP()
}

func (p Parameters) RingR() *ring.Ring {
	return p.paramsRP.RingQ()
}

func (p Parameters) RingRP() *rlwe.RingQP {
	return &rlwe.RingQP{p.RingR(), p.RingP()}
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
	return p.paramsQMulP.QCount()
}

func (p Parameters) N() int {
	return p.paramsQP.N()
}
