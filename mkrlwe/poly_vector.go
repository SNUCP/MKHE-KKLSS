package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"

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
	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.AddLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.AddLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// AddNoModLvl adds p1 to p2 coefficient-wise and writes the result on p3 without modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) AddNoModLvl(levelQ, levelP int, p1, p2, pOut *PolyQPVector) {

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.AddNoModLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.AddNoModLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// SubLvl subtracts p2 to p1 coefficient-wise and writes the result on p3.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) SubLvl(levelQ, levelP int, p1, p2, pOut *PolyQPVector) {

	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.SubLvl(levelQ, p1.Value[i].Q, p2.Value[i].Q, pOut.Value[i].Q)
		r.RingP.SubLvl(levelP, p1.Value[i].P, p2.Value[i].P, pOut.Value[i].P)
	}
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) NTTLvl(levelQ, levelP int, p, pOut *PolyQPVector) {
	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.NTTLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.NTTLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) InvNTTLvl(levelQ, levelP int, p, pOut *PolyQPVector) {
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
	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.NTTLazyLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.NTTLazyLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// MFormLvl switches p1 to the Montgomery domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MFormLvl(levelQ, levelP int, p, pOut *PolyQPVector) {
	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.MFormLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.MFormLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// InvMFormLvl switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) InvMFormLvl(levelQ, levelP int, p, pOut *PolyQPVector) {
	dim := p.Dim()

	for i := 0; i < dim; i++ {
		r.RingQ.InvMFormLvl(levelQ, p.Value[i].Q, pOut.Value[i].Q)
		r.RingP.InvMFormLvl(levelP, p.Value[i].P, pOut.Value[i].P)
	}
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery modular reduction.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) MulCoeffsMontgomeryLvl(levelQ, levelP int, p1, p2, p3 *PolyQPVector) {
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
	dim := p1.Dim()

	for i := 0; i < dim; i++ {
		ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, p1.Value[i].Q, index, p2.Value[i].Q)
		ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, p1.Value[i].P, index, p2.Value[i].P)
	}

}

// CopyValuesLvl copies the values of p1 on p2.
// The operation is performed at levelQ for the ringQ and levelP for the ringP.
func (r *RingQP) CopyValuesLvl(levelQ, levelP int, p1, p2 *PolyQPVector) {
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
