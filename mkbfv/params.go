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
	ringQMul    *ring.Ring
	ringR       *ring.Ring
	ringT       *ring.Ring
	CRS         [2]*mkrlwe.SwitchingKey
}

// NewParameters instantiate a set of MKCKKS parameters from the generic CKKS parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(pl ParametersLiteral) (params Parameters) {

	if len(pl.Q) != len(pl.QMul) {
		panic("cannot NewParametersFromLiteral: length of Q & QMul is not equal")
	}

	if len(pl.Q)%3 != 0 {
		panic("cannot NewParametersFromLiteral: length of Q is not a multiple of 3")
	}

	N := (1 << pl.LogN)
	R := make([]uint64, 0)
	R = append(R, pl.Q...)
	R = append(R, pl.QMul...)

	ringQMul, err := ring.NewRing(N, pl.QMul)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring QMul cannot be generated")
	}

	ringR, err := ring.NewRing(N, R)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring R cannot be generated")
	}

	ringT, err := ring.NewRing(N, []uint64{pl.T})
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring T cannot be generated")
	}

	params.ringQMul = ringQMul
	params.ringT = ringT
	params.ringR = ringR

	rlweParamsQP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q, P: pl.P, Sigma: pl.Sigma},
	)

	rlweParamsQMulP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.QMul, P: pl.P, Sigma: pl.Sigma},
	)

	rlweParamsRP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: R, P: pl.P, Sigma: pl.Sigma},
	)

	params.paramsQP = mkrlwe.NewParameters(rlweParamsQP, 3)
	params.paramsQMulP = mkrlwe.NewParameters(rlweParamsQMulP, 3)
	params.paramsRP = mkrlwe.NewParameters(rlweParamsRP, 3)

	params.CRS[0] = params.paramsRP.CRS[0]
	params.CRS[1] = params.paramsRP.CRS[1]

	conv := NewFastBasisExtender(params.paramsQP.RingP(), params.paramsQP.RingQ(),
		params.paramsQMulP.RingQ(), params.paramsRP.RingQ(),
	)

	// apply GadgetTransform
	conv.GadgetTransform(params.paramsQP.CRS[0], params.paramsQMulP.CRS[0], params.CRS[0])
	conv.GadgetTransform(params.paramsQP.CRS[1], params.paramsQMulP.CRS[1], params.CRS[1])

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

func (p Parameters) RingQMulP() *rlwe.RingQP {
	return p.paramsQMulP.RingQP()
}

func (p Parameters) RingR() *ring.Ring {
	return p.ringR
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
	return len(p.RingR().Modulus)
}

func (p Parameters) QMulCount() int {
	return len(p.RingQMul().Modulus)
}

func (p Parameters) N() int {
	return p.paramsQP.N()
}
