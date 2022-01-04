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

	polypoolQ    *ring.Poly
	polypoolQMul *ring.Poly
	mFormQMul    *ring.Poly
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

	conv.polypoolQ = ringQ.NewPoly()
	conv.polypoolQMul = ringQMul.NewPoly()

	conv.mFormQMul = ringQ.NewPoly()
	ringQ.AddScalarBigint(conv.mFormQMul, ringQMul.ModulusBigint, conv.mFormQMul)
	ringQ.MForm(conv.mFormQMul, conv.mFormQMul)

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

func (conv *FastBasisExtender) Quantize(polyR *ring.Poly, polyQ *ring.Poly, t uint64) {

	levelQ := len(conv.ringQ.Modulus) - 1

	/*
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
	*/

	for i := 0; i < levelQ+1; i++ {
		copy(polyQ.Coeffs[i], polyR.Coeffs[i])
	}

	conv.ringQ.InvNTT(polyQ, polyQ)
	conv.ringQ.MulScalar(polyQ, t, polyQ)

}

func (conv *FastBasisExtender) RescaleNTT(polyQ *ring.Poly, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringQ.MulCoeffsMontgomery(polyQ, conv.mFormQMul, conv.polypoolQ)
	conv.ringQMul.MulScalar(conv.polypoolQMul, 0, conv.polypoolQMul)
	conv.convQQMul.ModDownQPtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, conv.polypoolQMul)

	conv.ModUpQMultoRAndNTT(conv.polypoolQMul, polyR)

}

func (conv *FastBasisExtender) GadgetTransform(swkQP1, swkQP2, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swkQP1.Value)
	alpha := len(conv.ringQ.Modulus) / beta

	if beta != len(swkQP2.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if 2*beta != len(swkRP.Value) {
		panic("cannot GadgetTransform: swkRP have wrong dimension")
	}

	// apply InvMForm
	for i := 0; i < beta; i++ {
		conv.ringQ.InvMForm(swkQP1.Value[i].Q, swkQP1.Value[i].Q)
		conv.ringP.InvMForm(swkQP1.Value[i].P, swkQP1.Value[i].P)

		conv.ringQ.InvMForm(swkQP2.Value[i].Q, swkQP2.Value[i].Q)
		conv.ringP.InvMForm(swkQP2.Value[i].P, swkQP2.Value[i].P)
	}

	Q := conv.ringQ.ModulusBigint
	QMul := conv.ringQMul.ModulusBigint

	InvQ := big.NewInt(1)
	InvQ.Mul(InvQ, Q)
	InvQ.ModInverse(InvQ, QMul)

	InvQMul := big.NewInt(1)
	InvQMul.Mul(InvQMul, QMul)
	InvQMul.ModInverse(InvQMul, Q)

	//transform Q part
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(conv.ringQ.Modulus[i*alpha+j])))
		}

		InvQMulModQi := big.NewInt(1)
		InvQMulModQi.Mod(InvQMul, Qi)

		conv.modUpQPtoRP(swkQP1.Value[i], swkRP.Value[i])

		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, QMul, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, QMul, swkRP.Value[i].P)

		conv.modUpQPtoRP(swkQP2.Value[i], swkRP.Value[i+beta])

		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, QMul, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, QMul, swkRP.Value[i+beta].P)

	}

	// apply MForm
	for i := 0; i < beta; i++ {
		conv.ringQ.MForm(swkQP1.Value[i].Q, swkQP1.Value[i].Q)
		conv.ringP.MForm(swkQP1.Value[i].P, swkQP1.Value[i].P)

		conv.ringQ.MForm(swkQP2.Value[i].Q, swkQP2.Value[i].Q)
		conv.ringP.MForm(swkQP2.Value[i].P, swkQP2.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i].Q, swkRP.Value[i].Q)
		conv.ringP.MForm(swkRP.Value[i].P, swkRP.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i+beta].Q, swkRP.Value[i+beta].Q)
		conv.ringP.MForm(swkRP.Value[i+beta].P, swkRP.Value[i+beta].P)

	}

}

func (conv *FastBasisExtender) GadgetTransform2(swkQP1, swkQP2, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swkQP1.Value)
	alpha := len(conv.ringQ.Modulus) / beta

	if beta != len(swkQP2.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if 2*beta != len(swkRP.Value) {
		panic("cannot GadgetTransform: swkRP have wrong dimension")
	}

	// apply InvMForm
	for i := 0; i < beta; i++ {
		conv.ringQ.InvMForm(swkQP1.Value[i].Q, swkQP1.Value[i].Q)
		conv.ringP.InvMForm(swkQP1.Value[i].P, swkQP1.Value[i].P)

		conv.ringQ.InvMForm(swkQP2.Value[i].Q, swkQP2.Value[i].Q)
		conv.ringP.InvMForm(swkQP2.Value[i].P, swkQP2.Value[i].P)
	}

	Q := conv.ringQ.ModulusBigint
	QMul := conv.ringQMul.ModulusBigint

	InvQ := big.NewInt(1)
	InvQ.Mul(InvQ, Q)
	InvQ.ModInverse(InvQ, QMul)

	InvQMul := big.NewInt(1)
	InvQMul.Mul(InvQMul, QMul)
	InvQMul.ModInverse(InvQMul, Q)

	//transform Q part
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(conv.ringQ.Modulus[i*alpha+j])))
		}

		InvQMulModQi := big.NewInt(1)
		InvQMulModQi.Mod(InvQMul, Qi)

		conv.modUpQPtoRP(swkQP1.Value[i], swkRP.Value[i])

		/*
			conv.ringR.MulScalarBigint(swkRP.Value[i].Q, QMul, swkRP.Value[i].Q)
			conv.ringP.MulScalarBigint(swkRP.Value[i].P, QMul, swkRP.Value[i].P)
		*/

		conv.modUpQPtoRP(swkQP2.Value[i], swkRP.Value[i+beta])

		/*
			conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, QMul, swkRP.Value[i+beta].Q)
			conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, QMul, swkRP.Value[i+beta].P)
		*/
	}

	// apply MForm
	for i := 0; i < beta; i++ {
		conv.ringQ.MForm(swkQP1.Value[i].Q, swkQP1.Value[i].Q)
		conv.ringP.MForm(swkQP1.Value[i].P, swkQP1.Value[i].P)

		conv.ringQ.MForm(swkQP2.Value[i].Q, swkQP2.Value[i].Q)
		conv.ringP.MForm(swkQP2.Value[i].P, swkQP2.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i].Q, swkRP.Value[i].Q)
		conv.ringP.MForm(swkRP.Value[i].P, swkRP.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i+beta].Q, swkRP.Value[i+beta].Q)
		conv.ringP.MForm(swkRP.Value[i+beta].P, swkRP.Value[i+beta].P)

	}

}
