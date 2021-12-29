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
	mkrlwe.Parameters
	paramsRP mkrlwe.Parameters
	ringR    *ring.Ring
	ringQMul *ring.Ring
	ringT    *ring.Ring
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

	ringR, err := ring.NewRing(N, R)
	if err != nil {

		panic("cannot NewParametersFromLiteral: ring R cannot be generated")
	}
	params.ringR = ringR

	ringQMul, err := ring.NewRing(N, pl.QMul)
	if err != nil {

		panic("cannot NewParametersFromLiteral: ring QMul cannot be generated")
	}
	params.ringQMul = ringQMul

	rlweParamsQP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: pl.Q, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic("cannot NewParametersFromLiteral: ring QP cannot be generated")
	}
	params.Parameters = mkrlwe.NewParameters(rlweParamsQP, 2)

	rlweParamsRP, err := rlwe.NewParametersFromLiteral(
		rlwe.ParametersLiteral{LogN: pl.LogN, Q: R, P: pl.P, Sigma: pl.Sigma},
	)
	if err != nil {
		panic(err)
		panic("cannot NewParametersFromLiteral: ring RP cannot be generated")

	}
	params.paramsRP = mkrlwe.NewParameters(rlweParamsRP, 2)

	conv := NewFastBasisExtender(params.RingP(), params.RingQ(), ringQMul, ringR)
	conv.GadgetTransform(params.CRS[0], params.CRS[-3], params.paramsRP.CRS[0])
	conv.GadgetTransform(params.CRS[-1], params.CRS[-4], params.paramsRP.CRS[-1])

	return
}

func (p Parameters) RingQMul() *ring.Ring {
	return p.ringQMul
}

func (p Parameters) RingQMulP() *rlwe.RingQP {
	return &rlwe.RingQP{p.RingQMul(), p.RingP()}
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
