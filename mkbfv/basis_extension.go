package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "mk-lattigo/mkrlwe"

import "math/big"

type FastBasisExtender struct {
	ringP  *ring.Ring
	ringQ  *ring.Ring
	ringQ1 *ring.Ring
	ringR  *ring.Ring

	convQQ1 *mkrlwe.FastBasisExtender

	polypoolQ  *ring.Poly
	polypoolQ1 *ring.Poly
}

func NewFastBasisExtender(ringP, ringQ, ringQ1, ringR *ring.Ring) (conv *FastBasisExtender) {

	if len(ringQ.Modulus) != len(ringQ1.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringQ1 has different level")
	}

	if 2*len(ringQ.Modulus) != len(ringR.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringR has different level")
	}

	conv = new(FastBasisExtender)
	conv.ringP = ringP
	conv.ringQ = ringQ
	conv.ringQ1 = ringQ1
	conv.ringR = ringR
	conv.convQQ1 = mkrlwe.NewFastBasisExtender(ringQ, ringQ1)

	conv.polypoolQ = ringQ.NewPoly()
	conv.polypoolQ1 = ringQ1.NewPoly()

	return conv
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQPtoRP(polyQP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	conv.ringQ.InvNTT(polyQP.Q, conv.polypoolQ)

	conv.convQQ1.ModUpQtoP(levelQ, levelQ1, conv.polypoolQ, conv.polypoolQ1)

	conv.ringQ1.NTT(conv.polypoolQ1, conv.polypoolQ1)

	polyRP.P.Copy(polyQP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], polyQP.Q.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], conv.polypoolQ1.Coeffs[i])
	}
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQ1PtoRP(polyQ1P, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	conv.ringQ1.InvNTT(polyQ1P.Q, conv.polypoolQ1)

	conv.convQQ1.ModUpPtoQ(levelQ1, levelQ, conv.polypoolQ1, conv.polypoolQ)

	conv.ringQ.NTT(conv.polypoolQ, conv.polypoolQ)

	polyRP.P.Copy(polyQ1P.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], polyQ1P.Q.Coeffs[i])
	}

}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQtoRAndNTT(polyQ, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	conv.convQQ1.ModUpQtoP(levelQ, levelQ1, polyQ, conv.polypoolQ1)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], polyQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], conv.polypoolQ1.Coeffs[i])
	}

	conv.ringR.NTT(polyR, polyR)
}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQ1toRAndNTT(polyQ1, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	conv.convQQ1.ModUpPtoQ(levelQ1, levelQ, polyQ1, conv.polypoolQ)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], polyQ1.Coeffs[i])
	}

	conv.ringR.NTT(polyR, polyR)
}

func (conv *FastBasisExtender) Quantize(polyR *ring.Poly, polyQ *ring.Poly, t uint64) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(conv.polypoolQ1.Coeffs[i], polyR.Coeffs[i+levelQ+1])
	}

	conv.ringQ.InvNTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQ1.InvNTT(conv.polypoolQ1, conv.polypoolQ1)

	conv.ringQ.MulScalar(conv.polypoolQ, t, conv.polypoolQ)
	conv.ringQ1.MulScalar(conv.polypoolQ1, t, conv.polypoolQ1)

	conv.convQQ1.ModDownQPtoQ(levelQ, levelQ1, conv.polypoolQ, conv.polypoolQ1, polyQ)
}

func (conv *FastBasisExtender) RescaleNTT(polyQ *ring.Poly, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ

	conv.ModUpQtoRAndNTT(polyQ, polyR)

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(conv.polypoolQ1.Coeffs[i], polyR.Coeffs[i+levelQ+1])
	}

	conv.ringQ.InvNTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQ1.InvNTT(conv.polypoolQ1, conv.polypoolQ1)

	conv.ringQ.MulScalarBigint(conv.polypoolQ, conv.ringQ1.ModulusBigint, conv.polypoolQ)
	conv.ringQ1.MulScalarBigint(conv.polypoolQ1, conv.ringQ1.ModulusBigint, conv.polypoolQ1)
	conv.convQQ1.ModDownQPtoP(levelQ, levelQ1, conv.polypoolQ, conv.polypoolQ1, conv.polypoolQ1)

	conv.ModUpQ1toRAndNTT(conv.polypoolQ1, polyR)

}

func (conv *FastBasisExtender) GadgetTransform(swkQP, swkQ1P, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swkQP.Value)
	alpha := len(conv.ringQ.Modulus) / beta

	if beta != len(swkQ1P.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if 2*beta != len(swkRP.Value) {
		panic("cannot GadgetTransform: swkRP have wrong dimension")
	}

	// apply InvMForm
	for i := 0; i < beta; i++ {
		conv.ringQ.InvMForm(swkQP.Value[i].Q, swkQP.Value[i].Q)
		conv.ringP.InvMForm(swkQP.Value[i].P, swkQP.Value[i].P)

		conv.ringQ1.InvMForm(swkQ1P.Value[i].Q, swkQ1P.Value[i].Q)
		conv.ringP.InvMForm(swkQ1P.Value[i].P, swkQ1P.Value[i].P)
	}

	Q := conv.ringQ.ModulusBigint
	Q1 := conv.ringQ1.ModulusBigint

	InvQ := big.NewInt(1)
	InvQ.Mul(InvQ, Q)
	InvQ.ModInverse(InvQ, Q1)

	InvQ1 := big.NewInt(1)
	InvQ1.Mul(InvQ1, Q1)
	InvQ1.ModInverse(InvQ1, Q)

	//transform Q part
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(conv.ringQ.Modulus[i*alpha+j])))
		}

		InvQ1ModQi := big.NewInt(1)
		InvQ1ModQi.Mod(InvQ1, Qi)

		conv.modUpQPtoRP(swkQP.Value[i], swkRP.Value[i])
		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, InvQ1ModQi, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, InvQ1ModQi, swkRP.Value[i].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, Q1, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, Q1, swkRP.Value[i].P)

	}

	//transform Q1 part
	for i := 0; i < beta; i++ {
		Q1i := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Q1i.Mul(Q1i, big.NewInt(int64(conv.ringQ1.Modulus[i*alpha+j])))
		}

		InvQModQ1i := big.NewInt(1)
		InvQModQ1i.Mod(InvQ, Q1i)

		conv.modUpQ1PtoRP(swkQ1P.Value[i], swkRP.Value[i+beta])
		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, InvQModQ1i, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, InvQModQ1i, swkRP.Value[i+beta].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, Q, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, Q, swkRP.Value[i+beta].P)
	}

	// apply MForm
	for i := 0; i < beta; i++ {
		conv.ringQ.MForm(swkQP.Value[i].Q, swkQP.Value[i].Q)
		conv.ringP.MForm(swkQP.Value[i].P, swkQP.Value[i].P)

		conv.ringQ1.MForm(swkQ1P.Value[i].Q, swkQ1P.Value[i].Q)
		conv.ringP.MForm(swkQ1P.Value[i].P, swkQ1P.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i].Q, swkRP.Value[i].Q)
		conv.ringP.MForm(swkRP.Value[i].P, swkRP.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i+beta].Q, swkRP.Value[i+beta].Q)
		conv.ringP.MForm(swkRP.Value[i+beta].P, swkRP.Value[i+beta].P)

	}

}
