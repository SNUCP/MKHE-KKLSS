package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "math/big"
import "math"

type PolyQPVector struct {
	Value []rlwe.PolyQP
}

// Dim returns dimension of the target polynomial vector
func (polvec *PolyQPVector) Dim() int {
	return len(polvec.Value)
}

// SetEntry sets entry of given index with given poly
func (polvec *PolyQPVector) SetEntry(idx int, pol rlwe.PolyQP) {
	if idx >= polvec.Dim() {
		return
	}

	polvec.Value[idx] = pol
}

// Equals returns true if the receiver PolyQPVector is equal to the provided other PolyQPVector.
// This method checks for equality of its two sub-polynomials.
func (polvec *PolyQPVector) Equals(other *PolyQPVector) bool {
	if polvec.Dim() != other.Dim() {
		return false
	}

	if polvec == other {
		return true
	}

	for i := 0; i < polvec.Dim(); i++ {
		if !polvec.Value[i].Equals(other.Value[i]) {
			return false
		}
	}

	return true
}

// CopyNew creates an exact copy of the target polynomial.
func (polvec *PolyQPVector) CopyNew() *PolyQPVector {
	if polvec == nil {
		return nil
	}

	ret := new(PolyQPVector)
	dim := polvec.Dim()
	ret.Value = make([]rlwe.PolyQP, dim)

	for i := 0; i < dim; i++ {
		ret.Value[i].Q = polvec.Value[i].Q.CopyNew()
		ret.Value[i].P = polvec.Value[i].P.CopyNew()
	}

	return ret
}

type RingQP struct {
	rlwe.RingQP
}

// NewRingQP returns a new RingQP
func NewRingQP(r rlwe.RingQP) *RingQP {
	ret := new(RingQP)
	ret.RingP = r.RingP
	ret.RingQ = r.RingQ
	return ret
}

// NewPoly creates a new polynomial vector of given dimension with all coefficients set to 0.
func (r *RingQP) NewPolyVector(dim int) *PolyQPVector {
	polvec := new(PolyQPVector)
	polvec.Value = make([]rlwe.PolyQP, dim)

	for i := 0; i < dim; i++ {
		polvec.Value[i] = r.NewPoly()
	}

	return polvec
}

// NewPolyVecLvl creates a new polynomial vector of given dimension with all coefficients set to 0.
func (r *RingQP) NewPolyVecLvl(dim, levelQ, levelP int) *PolyQPVector {

	polvec := new(PolyQPVector)
	polvec.Value = make([]rlwe.PolyQP, dim)

	for i := 0; i < dim; i++ {
		polvec.Value[i] = r.NewPolyLvl(levelQ, levelP)
	}

	return polvec
}

// AddLvl adds p1 to p2 coefficient-wise and writes the result on p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) AddLvl(levelQ, levelP int, p1, p2, pOut *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot AddLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != pOut.Dim() {
		panic("cannot AddLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.AddLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.AddLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// AddNoModLvl adds p1 to p2 coefficient-wise and writes the result on p3 without modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) AddNoModLvl(levelQ, levelP int, p1, p2, pOut *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot AddNoModLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != pOut.Dim() {
		panic("cannot AddNoModLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.AddNoModLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.AddNoModLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// SubLvl subtracts p2 to p1 coefficient-wise and writes the result on p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) SubLvl(levelQ, levelP int, p1, p2, pOut *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot SubLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != pOut.Dim() {
		panic("cannot SubLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.SubLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.SubLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) NTTLvl(levelQ, levelP int, p, pOut *PolyQPVector) {

	if p.Dim() != pOut.Dim() {
		panic("cannot NTTLvl: input poly vectors have different dimensions")
	}

	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.NTTLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.NTTLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) InvNTTLvl(levelQ, levelP int, p, pOut *PolyQPVector) {

	if p.Dim() != pOut.Dim() {
		panic("cannot InvNTTLvl: input poly vectors have different dimensions")
	}

	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.InvNTTLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.InvNTTLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// NTTLazyLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
// Output values are in the range [0, 2q-1].
func (r *RingQP) NTTLazyLvl(levelQ, levelP int, p, pOut *PolyQPVector) {

	if p.Dim() != pOut.Dim() {
		panic("cannot NTTLazyLvl: input poly vectors have different dimensions")
	}

	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.NTTLazyLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.NTTLazyLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// MFormLvl switches p1 to the Montgomery domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MFormLvl(levelQ, levelP int, p, pOut *PolyQPVector) {

	if p.Dim() != pOut.Dim() {
		panic("cannot MFormLvl: input poly vectors have different dimensions")
	}

	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.MFormLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.MFormLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// InvMFormLvl switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) InvMFormLvl(levelQ, levelP int, p, pOut *PolyQPVector) {

	if p.Dim() != pOut.Dim() {
		panic("cannot InvMFormLvl: input poly vectors have different dimensions")
	}

	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.InvMFormLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.InvMFormLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MulCoeffsMontgomeryLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulCoeffsMontgomeryLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != p3.Dim() {
		panic("cannot MulCoeffsMontgomeryLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		r.RingQ.MulCoeffsMontgomeryLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, p3.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryLvl(levelP, p1.Value[i].P, p2.Value[i].P, p3.Value[i].P)
	}
}

// MulCoeffsMontgomeryConstantLvl multiplies p1 by p2 coefficient-wise with a constant-time Montgomery modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
// Result is within [0, 2q-1].
func (r *RingQP) MulCoeffsMontgomeryConstantLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulCoeffsMontgomeryConstantLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != p3.Dim() {
		panic("cannot MulCoeffsMontgomeryConstantLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		r.RingQ.MulCoeffsMontgomeryConstantLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, p3.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryConstantLvl(levelP, p1.Value[i].P, p2.Value[i].P, p3.Value[i].P)
	}
}

// MulCoeffsMontgomeryConstantAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and adds the result on p3.
// Result is within [0, 2q-1]
func (r *RingQP) MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulCoeffsMontgomeryConstantAndAddNoModLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != p3.Dim() {
		panic("cannot MulCoeffsMontgomeryConstantAndAddNoModLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		r.RingQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, p3.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelP, p1.Value[i].P, p2.Value[i].P, p3.Value[i].P)
	}
}

// MulCoeffsMontgomeryAndSubLvl multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MulCoeffsMontgomeryAndSubLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulCoeffsMontgomeryAndSubLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != p3.Dim() {
		panic("cannot MulCoeffsMontgomeryAndSubLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		r.RingQ.MulCoeffsMontgomeryAndSubLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, p3.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryAndSubLvl(levelP, p1.Value[i].P, p2.Value[i].P, p3.Value[i].P)
	}
}

// MulCoeffsMontgomeryAndAddLvl multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MulCoeffsMontgomeryAndAddLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulCoeffsMontgomeryAndAddLvl: input poly vectors have different dimensions")
	}

	if p2.Dim() != p3.Dim() {
		panic("cannot MulCoeffsMontgomeryAndAddLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		r.RingQ.MulCoeffsMontgomeryAndAddLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, p3.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryAndAddLvl(levelP, p1.Value[i].P, p2.Value[i].P, p3.Value[i].P)
	}
}

// PermuteNTTWithIndexLvl applies the automorphism X^{5^j} on p1 and writes the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) PermuteNTTWithIndexLvl(levelQ, levelP int, p1 *PolyQPVector, index []uint64, p2 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot PermuteNTTWithIndexLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		ring.PermuteNTTWithIndexLvl(levelQ, p1.Value[i].Q, index, p2.Value[i].Q)
		ring.PermuteNTTWithIndexLvl(levelP, p1.Value[i].P, index, p2.Value[i].P)
	}

}

// PermuteNTTWithIndexAndAddNoModLvl applies the automorphism X^{5^j} on p1 and adds the result on p2.
// Index of automorphism must be provided.
// Method is not in place.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) PermuteNTTWithIndexAndAddNoModLvl(levelQ, levelP int, p1 *PolyQPVector, index []uint64, p2 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot PermuteNTTWithIndexAndAddNoModLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, p1.Value[i].Q, index, p2.Value[i].Q)
		ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, p1.Value[i].P, index, p2.Value[i].P)
	}
}

// MulPolyMontgomeryLvl multiplies each entry of p1 by a polynomial b and writes the result on p2.
// Inputs should be in NTT form, and one of them should be in MForm
// reduction algorithm uses Montgomery reduction
func (r *RingQP) MulPolyMontgomeryLvl(levelQ, levelP int, p1 *PolyQPVector, b *rlwe.PolyQP, p2 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulPoly: input and output poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.MulCoeffsMontgomeryLvl(levelQ, b.Q, p2.Value[i].Q, p2.Value[i].Q)
		r.RingP.MulCoeffsMontgomeryLvl(levelP, b.P, p2.Value[i].P, p2.Value[i].P)
	}
}

// MulScalarBigintLvl multiplies each entry of p1 by a bigint b and writes the result on p2.
// reduction algorithm uses Montgomery reduction
func (r *RingQP) MulScalarBigintLvl(levelQ, levelP int, p1 *PolyQPVector, b *big.Int, p2 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot MulScalarBigintLvl: input and output poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.MulScalarBigintLvl(levelQ, p1.Value[i].Q, b, p2.Value[i].Q)
		r.RingP.MulScalarBigintLvl(levelP, p1.Value[i].P, b, p2.Value[i].P)
	}
}

// InternalProductLvl performs internal product to input ringQ polynomial & polyVector
func (r *RingQP) InternalProductLvl(ks *rlwe.KeySwitcher, levelQ int, cx *ring.Poly, polvec *PolyQPVector, c0Q, c0P, c1Q, c1P *ring.Poly) {

	ringQ := ks.RingQ()
	ringP := ks.RingP()
	ringQP := ks.RingQP()

	c2QP := ks.Pool[0]

	var cxNTT, cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxNTT = cx
		cxInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		cxNTT = ks.PoolInvNTT
		cxInvNTT = cx
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	alpha := len(polvec.Value[0].P.Coeffs)
	levelP := alpha - 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1
	reduce := 0
	c0QP := rlwe.PolyQP{c0Q, c0P}

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {
		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, polvec.Value[i], c2QP, c0QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, polvec.Value[i], c2QP, c0QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
	}
}

// CopyValuesLvl copies the values of p1 on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) CopyValuesLvl(levelQ, levelP int, p1, p2 *PolyQPVector) {

	if p1.Dim() != p2.Dim() {
		panic("cannot CopyValuesLvl: input poly vectors have different dimensions")
	}

	dim := p1.Dim()

	for i := 0; i < dim; i++ {

		ring.CopyValuesLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q)
		ring.CopyValuesLvl(levelP, p1.Value[i].P, p2.Value[i].P)
	}
}

// Copy copies the input polyQP on the target polyQP.
func (p *PolyQPVector) Copy(polFrom *PolyQPVector) {
	dim := p.Dim()
	for i := 0; i < dim; i++ {
		p.Value[i].Q.Copy(polFrom.Value[i].Q)
		p.Value[i].P.Copy(polFrom.Value[i].P)
	}
}
