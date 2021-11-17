package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
type ParametersLiteral struct {
	LogN  int // Log Ring degree (power of 2)
	Q     []uint64
	Q1    []uint64
	P     []uint64
	LogQ  []int   `json:",omitempty"`
	LogQ1 []int   `json:",omitempty"`
	LogP  []int   `json:",omitempty"`
	Sigma float64 // Gaussian sampling standard deviation
	T     uint64  // Plaintext modulus
}

// Parameters represents a parameter set for the BFV cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
// R = Q*QMul
type Parameters struct {
	paramsQP  mkrlwe.Parameters
	paramsQ1P mkrlwe.Parameters
	paramsRP  mkrlwe.Parameters
	ringT     *ring.Ring
	CRS       [2]*mkrlwe.SwitchingKey
}

// NewParameters instantiate a set of MKCKKS parameters from the generic CKKS parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(pl ParametersLiteral) (params Parameters) {

	if len(pl.Q) != len(pl.Q1) {
		panic("cannot NewParametersFromLiteral: length of Q & Q1 is not equal")
	}

	N := (1 << pl.LogN)
	R := make([]uint64, 0)
	R = append(R, pl.Q...)
	R = append(R, pl.Q1...)

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

	rlweParamsQ1P, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q1, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring Q1P cannot be generated")
	}

	rlweParamsRP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: R, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic(err)
		panic("cannot NewParametersFromLiteral: ring RP cannot be generated")

	}

	params.paramsQP = mkrlwe.NewParameters(rlweParamsQP, 3)
	params.paramsQ1P = mkrlwe.NewParameters(rlweParamsQ1P, 3)
	params.paramsRP = mkrlwe.NewParameters(rlweParamsRP, 3)

	params.CRS[0] = params.paramsRP.CRS[0]
	params.CRS[1] = params.paramsRP.CRS[1]

	conv := NewFastBasisExtender(
		params.paramsQP.RingP(), params.paramsQP.RingQ(),
		params.paramsQ1P.RingQ(),
		params.paramsRP.RingQ(),
	)

	// apply GadgetTransform
	conv.GadgetTransform(params.paramsQP.CRS[0], params.paramsQ1P.CRS[0], params.CRS[0])
	conv.GadgetTransform(params.paramsQP.CRS[1], params.paramsQ1P.CRS[1], params.CRS[1])

	return params
}

func (p Parameters) RingQ1() *ring.Ring {
	return p.paramsQ1P.RingQ()
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

func (p Parameters) RingQ1P() *rlwe.RingQP {
	return p.paramsQ1P.RingQP()
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

func (p Parameters) Q1Count() int {
	return p.paramsQ1P.QCount()
}

func (p Parameters) N() int {
	return p.paramsQP.N()
}
