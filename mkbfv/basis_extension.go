package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "mk-lattigo/mkrlwe"

import "math/big"

type FastBasisExtender struct {
	ringP  *ring.Ring
	ringQ  *ring.Ring
	ringQ1 *ring.Ring
	ringQ2 *ring.Ring
	ringR  *ring.Ring

	convQQ1  *mkrlwe.FastBasisExtender
	convQQ2  *mkrlwe.FastBasisExtender
	convQ1Q2 *mkrlwe.FastBasisExtender

	polypoolQ  *ring.Poly
	polypoolQ1 *ring.Poly
	polypoolQ2 *ring.Poly

	quant_factor *big.Int
}

func NewFastBasisExtender(ringP, ringQ, ringQ1, ringQ2, ringR *ring.Ring, t uint64) (conv *FastBasisExtender) {

	if len(ringQ.Modulus) != len(ringQ1.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringQ1 has different level")
	}

	if len(ringQ1.Modulus) != len(ringQ2.Modulus) {
		panic("cannot NewFastBasisExtender ringQ1 and ringQ2 has different level")
	}

	if 3*len(ringQ.Modulus) != len(ringR.Modulus) {
		panic("cannot NewFastBasisExtender ringQ and ringR has different level")
	}

	conv = new(FastBasisExtender)
	conv.ringP = ringP
	conv.ringQ = ringQ
	conv.ringQ1 = ringQ1
	conv.ringQ2 = ringQ2
	conv.ringR = ringR
	conv.convQQ1 = mkrlwe.NewFastBasisExtender(ringQ, ringQ1)
	conv.convQQ2 = mkrlwe.NewFastBasisExtender(ringQ, ringQ2)
	conv.convQ1Q2 = mkrlwe.NewFastBasisExtender(ringQ1, ringQ2)

	conv.polypoolQ = ringQ.NewPoly()
	conv.polypoolQ1 = ringQ1.NewPoly()
	conv.polypoolQ2 = ringQ2.NewPoly()

	// gen quantization factor
	conv.quant_factor = big.NewInt(int64(t))
	conv.quant_factor.Mul(conv.quant_factor, ringQ1.ModulusBigint)
	conv.quant_factor.Mul(conv.quant_factor, ringQ2.ModulusBigint)
	conv.quant_factor.Div(conv.quant_factor, ringQ.ModulusBigint)

	return conv
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQPtoRP(polyQP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ
	levelQ2 := levelQ

	conv.ringQ.InvNTT(polyQP.Q, conv.polypoolQ)

	conv.convQQ1.ModUpQtoP(levelQ, levelQ1, conv.polypoolQ, conv.polypoolQ1)
	conv.convQQ2.ModUpQtoP(levelQ, levelQ2, conv.polypoolQ, conv.polypoolQ2)

	conv.ringQ1.NTT(conv.polypoolQ1, conv.polypoolQ1)
	conv.ringQ2.NTT(conv.polypoolQ2, conv.polypoolQ2)

	polyRP.P.Copy(polyQP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], polyQP.Q.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], conv.polypoolQ1.Coeffs[i])
	}

	for i := 0; i < levelQ2+1; i++ {
		copy(polyRP.Q.Coeffs[i+2*(levelQ+1)], conv.polypoolQ2.Coeffs[i])
	}
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQ1PtoRP(polyQ1P, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ
	levelQ2 := levelQ

	conv.ringQ1.InvNTT(polyQ1P.Q, conv.polypoolQ1)

	conv.convQQ1.ModUpPtoQ(levelQ1, levelQ, conv.polypoolQ1, conv.polypoolQ)
	conv.convQ1Q2.ModUpQtoP(levelQ1, levelQ2, conv.polypoolQ1, conv.polypoolQ2)

	conv.ringQ.NTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQ2.NTT(conv.polypoolQ2, conv.polypoolQ2)

	polyRP.P.Copy(polyQ1P.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], polyQ1P.Q.Coeffs[i])
	}

	for i := 0; i < levelQ2+1; i++ {
		copy(polyRP.Q.Coeffs[i+2*(levelQ+1)], conv.polypoolQ2.Coeffs[i])
	}
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQ2PtoRP(polyQ2P, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ
	levelQ2 := levelQ

	conv.ringQ2.InvNTT(polyQ2P.Q, conv.polypoolQ2)

	conv.convQQ2.ModUpPtoQ(levelQ2, levelQ, conv.polypoolQ2, conv.polypoolQ)
	conv.convQ1Q2.ModUpPtoQ(levelQ2, levelQ1, conv.polypoolQ2, conv.polypoolQ1)

	conv.ringQ.NTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQ1.NTT(conv.polypoolQ1, conv.polypoolQ1)

	polyRP.P.Copy(polyQ2P.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyRP.Q.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyRP.Q.Coeffs[i+levelQ+1], conv.polypoolQ1.Coeffs[i])
	}

	for i := 0; i < levelQ2+1; i++ {
		copy(polyRP.Q.Coeffs[i+2*(levelQ+1)], polyQ2P.Q.Coeffs[i])
	}
}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQtoRAndNTT(polyQ, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ
	levelQ2 := levelQ

	conv.convQQ1.ModUpQtoP(levelQ, levelQ1, polyQ, conv.polypoolQ1)
	conv.convQQ2.ModUpQtoP(levelQ, levelQ2, polyQ, conv.polypoolQ2)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], polyQ.Coeffs[i])
	}

	for i := 0; i < levelQ1+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], conv.polypoolQ1.Coeffs[i])
	}

	for i := 0; i < levelQ2+1; i++ {
		copy(polyR.Coeffs[i+2*(levelQ+1)], conv.polypoolQ2.Coeffs[i])
	}

	conv.ringR.NTT(polyR, polyR)
}

func (conv *FastBasisExtender) Quantize(polyR *ring.Poly, polyQ *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQ1 := levelQ
	levelQ2 := levelQ

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ1.Coeffs[i], polyR.Coeffs[i+levelQ+1])
	}

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ2.Coeffs[i], polyR.Coeffs[i+2*(levelQ+1)])
	}

	// Aplies inverse NTT
	conv.ringQ.InvNTT(conv.polypoolQ, conv.polypoolQ)
	conv.ringQ1.InvNTT(conv.polypoolQ1, conv.polypoolQ1)
	conv.ringQ2.InvNTT(conv.polypoolQ2, conv.polypoolQ2)

	// Multiply quant_factor

	conv.ringQ.MulScalarBigint(conv.polypoolQ, conv.quant_factor, conv.polypoolQ)
	conv.ringQ1.MulScalarBigint(conv.polypoolQ1, conv.quant_factor, conv.polypoolQ1)
	conv.ringQ2.MulScalarBigint(conv.polypoolQ2, conv.quant_factor, conv.polypoolQ2)

	// ModDown by Q2
	conv.convQQ2.ModDownQPtoP(levelQ, levelQ2, conv.polypoolQ, conv.polypoolQ2, conv.polypoolQ)
	conv.convQ1Q2.ModDownQPtoP(levelQ1, levelQ2, conv.polypoolQ1, conv.polypoolQ2, conv.polypoolQ1)

	// ModDown by Q1
	conv.convQQ1.ModDownQPtoP(levelQ, levelQ1, conv.polypoolQ, conv.polypoolQ1, polyQ)

}

func (conv *FastBasisExtender) GadgetTransform(swkQP, swkQ1P, swkQ2P, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swkQP.Value)
	alpha := len(conv.ringQ.Modulus) / beta

	if beta != len(swkQ1P.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if beta != len(swkQ2P.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if 3*beta != len(swkRP.Value) {
		panic("cannot GadgetTransform: swkRP have wrong dimension")
	}

	// apply InvMForm
	for i := 0; i < beta; i++ {
		conv.ringQ.InvMForm(swkQP.Value[i].Q, swkQP.Value[i].Q)
		conv.ringP.InvMForm(swkQP.Value[i].P, swkQP.Value[i].P)

		conv.ringQ1.InvMForm(swkQ1P.Value[i].Q, swkQ1P.Value[i].Q)
		conv.ringP.InvMForm(swkQ1P.Value[i].P, swkQ1P.Value[i].P)

		conv.ringQ2.InvMForm(swkQ2P.Value[i].Q, swkQ2P.Value[i].Q)
		conv.ringP.InvMForm(swkQ2P.Value[i].P, swkQ2P.Value[i].P)
	}

	Q := conv.ringQ.ModulusBigint
	Q1 := conv.ringQ1.ModulusBigint
	Q2 := conv.ringQ2.ModulusBigint

	InvQQ1 := big.NewInt(1)
	InvQQ1.Mul(InvQQ1, Q)
	InvQQ1.Mul(InvQQ1, Q1)
	InvQQ1.ModInverse(InvQQ1, Q2)

	InvQQ2 := big.NewInt(1)
	InvQQ2.Mul(InvQQ2, Q)
	InvQQ2.Mul(InvQQ2, Q2)
	InvQQ2.ModInverse(InvQQ2, Q1)

	InvQ1Q2 := big.NewInt(1)
	InvQ1Q2.Mul(InvQ1Q2, Q1)
	InvQ1Q2.Mul(InvQ1Q2, Q2)
	InvQ1Q2.ModInverse(InvQ1Q2, Q)

	//transform Q part
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(conv.ringQ.Modulus[i*alpha+j])))
		}

		InvQ1Q2ModQi := big.NewInt(1)
		InvQ1Q2ModQi.Mod(InvQ1Q2, Qi)

		conv.modUpQPtoRP(swkQP.Value[i], swkRP.Value[i])
		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, InvQ1Q2ModQi, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, InvQ1Q2ModQi, swkRP.Value[i].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, Q1, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, Q1, swkRP.Value[i].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, Q2, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, Q2, swkRP.Value[i].P)
	}

	//transform Q1 part
	for i := 0; i < beta; i++ {
		Q1i := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Q1i.Mul(Q1i, big.NewInt(int64(conv.ringQ1.Modulus[i*alpha+j])))
		}

		InvQQ2ModQ1i := big.NewInt(1)
		InvQQ2ModQ1i.Mod(InvQQ2, Q1i)

		conv.modUpQ1PtoRP(swkQ1P.Value[i], swkRP.Value[i+beta])
		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, InvQQ2ModQ1i, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, InvQQ2ModQ1i, swkRP.Value[i+beta].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, Q, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, Q, swkRP.Value[i+beta].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, Q2, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, Q2, swkRP.Value[i+beta].P)
	}

	//transform Q2 part
	for i := 0; i < beta; i++ {
		Q2i := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Q2i.Mul(Q2i, big.NewInt(int64(conv.ringQ2.Modulus[i*alpha+j])))
		}

		InvQQ1ModQ2i := big.NewInt(1)
		InvQQ1ModQ2i.Mod(InvQQ1, Q2i)

		conv.modUpQ2PtoRP(swkQ2P.Value[i], swkRP.Value[i+2*beta])
		conv.ringR.MulScalarBigint(swkRP.Value[i+2*beta].Q, InvQQ1ModQ2i, swkRP.Value[i+2*beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+2*beta].P, InvQQ1ModQ2i, swkRP.Value[i+2*beta].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i+2*beta].Q, Q, swkRP.Value[i+2*beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+2*beta].P, Q, swkRP.Value[i+2*beta].P)

		conv.ringR.MulScalarBigint(swkRP.Value[i+2*beta].Q, Q1, swkRP.Value[i+2*beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+2*beta].P, Q1, swkRP.Value[i+2*beta].P)
	}

	// apply MForm
	for i := 0; i < beta; i++ {
		conv.ringQ.MForm(swkQP.Value[i].Q, swkQP.Value[i].Q)
		conv.ringP.MForm(swkQP.Value[i].P, swkQP.Value[i].P)

		conv.ringQ1.MForm(swkQ1P.Value[i].Q, swkQ1P.Value[i].Q)
		conv.ringP.MForm(swkQ1P.Value[i].P, swkQ1P.Value[i].P)

		conv.ringQ2.MForm(swkQ2P.Value[i].Q, swkQ2P.Value[i].Q)
		conv.ringP.MForm(swkQ2P.Value[i].P, swkQ2P.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i].Q, swkRP.Value[i].Q)
		conv.ringP.MForm(swkRP.Value[i].P, swkRP.Value[i].P)

		conv.ringR.MForm(swkRP.Value[i+beta].Q, swkRP.Value[i+beta].Q)
		conv.ringP.MForm(swkRP.Value[i+beta].P, swkRP.Value[i+beta].P)

		conv.ringR.MForm(swkRP.Value[i+2*beta].Q, swkRP.Value[i+2*beta].Q)
		conv.ringP.MForm(swkRP.Value[i+2*beta].P, swkRP.Value[i+2*beta].P)

	}

}
