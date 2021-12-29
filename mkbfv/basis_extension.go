package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "mk-lattigo/mkrlwe"

import "math/big"

type FastBasisExtender struct {
	ringP    *ring.Ring
	ringQ    *ring.Ring
	ringQMul *ring.Ring
	ringR    *ring.Ring

	convQQMul *mkrlwe.FastBasisExtender
	convPQMul *mkrlwe.FastBasisExtender

	polypoolQ    *ring.Poly
	polypoolP    *ring.Poly
	polypoolQMul *ring.Poly

	mFormQMul  *ring.Poly
	QMulBigInt *big.Int
}

func NewFastBasisExtender(ringP, ringQ, ringQMul, ringR *ring.Ring) (conv *FastBasisExtender) {

	if len(ringQ.Modulus) != len(ringQMul.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringQMul has different level")
	}

	if 2*len(ringQ.Modulus) != len(ringR.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringR has different level")
	}

	conv = new(FastBasisExtender)
	conv.ringP = ringP
	conv.ringQ = ringQ
	conv.ringQMul = ringQMul
	conv.ringR = ringR
	conv.convQQMul = mkrlwe.NewFastBasisExtender(ringQ, ringQMul)
	conv.convPQMul = mkrlwe.NewFastBasisExtender(ringP, ringQMul)

	conv.polypoolQ = ringQ.NewPoly()
	conv.polypoolP = ringP.NewPoly()
	conv.polypoolQMul = ringQMul.NewPoly()

	conv.mFormQMul = ringQ.NewPoly()
	ringQ.AddScalarBigint(conv.mFormQMul, ringQMul.ModulusBigint, conv.mFormQMul)
	ringQ.MForm(conv.mFormQMul, conv.mFormQMul)

	conv.QMulBigInt = ringQMul.ModulusBigint

	return conv
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQPtoRP(polyQP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringQ.InvNTT(polyQP.Q, conv.polypoolQ)

	conv.convQQMul.ModUpQtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul)

	conv.ringQMul.NTT(conv.polypoolQMul, conv.polypoolQMul)

	polyRP.P.Copy(polyQP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], polyQP.Q.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], conv.polypoolQMul.Coeffs[i])
	}
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQMulPtoRP(polyQMulP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringQMul.InvNTT(polyQMulP.Q, conv.polypoolQMul)

	conv.convQQMul.ModUpPtoQ(levelQMul, levelQ, conv.polypoolQMul, conv.polypoolQ)

	conv.ringQ.NTT(conv.polypoolQ, conv.polypoolQ)

	polyRP.P.Copy(polyQMulP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], polyQMulP.Q.Coeffs[i])
	}

}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQtoRAndNTT(polyQ, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.convQQMul.ModUpQtoP(levelQ, levelQMul, polyQ, conv.polypoolQMul)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], polyQ.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], conv.polypoolQMul.Coeffs[i])
	}

	conv.ringR.NTT(polyR, polyR)
}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQMultoRAndNTT(polyQMul, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.convQQMul.ModUpPtoQ(levelQMul, levelQ, polyQMul, conv.polypoolQ)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], polyQMul.Coeffs[i])
	}

	conv.ringR.NTT(polyR, polyR)
}

// Assume input polyR is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModDownRPtoQP(polyRP, polyQP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := len(conv.ringQMul.Modulus) - 1
	levelP := len(conv.ringP.Modulus) - 1

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyRP.Q.Coeffs[i])
		copy(conv.polypoolQMul.Coeffs[i], polyRP.Q.Coeffs[i+levelQ+1])
	}

	for i := 0; i < levelP+1; i++ {
		copy(conv.polypoolP.Coeffs[i], polyRP.P.Coeffs[i])
	}

	conv.convQQMul.ModDownQPtoQ(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, polyQP.Q)
	conv.convPQMul.ModDownQPtoQ(levelP, levelQMul, conv.polypoolP, conv.polypoolQMul, polyQP.P)

}

func (conv *FastBasisExtender) Quantize(polyR *ring.Poly, polyQ *ring.Poly, t uint64) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(conv.polypoolQMul.Coeffs[i], polyR.Coeffs[i+levelQ+1])
	}

	conv.ringQ.InvNTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQMul.InvNTT(conv.polypoolQMul, conv.polypoolQMul)

	conv.ringQ.MulScalar(conv.polypoolQ, t, conv.polypoolQ)
	conv.ringQMul.MulScalar(conv.polypoolQMul, t, conv.polypoolQMul)

	conv.convQQMul.ModDownQPtoQ(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, polyQ)
}

func (conv *FastBasisExtender) RescaleNTT(polyQ *ring.Poly, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringQ.MulCoeffsMontgomery(polyQ, conv.mFormQMul, conv.polypoolQ)
	conv.ringQMul.MulScalar(conv.polypoolQMul, 0, conv.polypoolQMul)
	conv.convQQMul.ModDownQPtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, conv.polypoolQMul)

	conv.ModUpQMultoRAndNTT(conv.polypoolQMul, polyR)
}

func (conv *FastBasisExtender) GadgetTransform(swk1, swk2, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swk1.Value)

	for i := 0; i < beta; i++ {
		conv.ringQ.InvMForm(swk1.Value[i].Q, swk1.Value[i].Q)
		conv.ringQ.InvMForm(swk2.Value[i].Q, swk2.Value[i].Q)

		conv.ringP.InvMForm(swk1.Value[i].P, swk1.Value[i].P)
		conv.ringP.InvMForm(swk2.Value[i].P, swk2.Value[i].P)

		conv.ringQ.InvNTT(swk1.Value[i].Q, swk1.Value[i].Q)
		conv.ringQ.InvNTT(swk2.Value[i].Q, swk2.Value[i].Q)
	}

	for i := 0; i < beta; i++ {
		conv.ModUpQtoRAndNTT(swk1.Value[i].Q, swkRP.Value[i].Q)
		swkRP.Value[i].P.Copy(swk1.Value[i].P)
		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, conv.QMulBigInt, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, conv.QMulBigInt, swkRP.Value[i].P)
		conv.ringR.MForm(swkRP.Value[i].Q, swkRP.Value[i].Q)
		conv.ringP.MForm(swkRP.Value[i].P, swkRP.Value[i].P)

		conv.ModUpQtoRAndNTT(swk2.Value[i].Q, swkRP.Value[i+beta].Q)
		swkRP.Value[i+beta].P.Copy(swk2.Value[i].P)
		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, conv.QMulBigInt, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, conv.QMulBigInt, swkRP.Value[i+beta].P)
		conv.ringR.MForm(swkRP.Value[i+beta].Q, swkRP.Value[i+beta].Q)
		conv.ringP.MForm(swkRP.Value[i+beta].P, swkRP.Value[i+beta].P)
	}

	for i := 0; i < beta; i++ {
		conv.ringQ.NTT(swk1.Value[i].Q, swk1.Value[i].Q)
		conv.ringQ.NTT(swk2.Value[i].Q, swk2.Value[i].Q)

		conv.ringP.MForm(swk1.Value[i].P, swk1.Value[i].P)
		conv.ringP.MForm(swk2.Value[i].P, swk2.Value[i].P)

		conv.ringQ.MForm(swk1.Value[i].Q, swk1.Value[i].Q)
		conv.ringQ.MForm(swk2.Value[i].Q, swk2.Value[i].Q)
	}

}
