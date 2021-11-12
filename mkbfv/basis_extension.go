package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "mk-lattigo/mkrlwe"

import "math/big"

type FastBasisExtender struct {
	baseconverter *mkrlwe.FastBasisExtender
	ringP         *ring.Ring
	ringQ         *ring.Ring
	ringQMul      *ring.Ring
	ringR         *ring.Ring

	QMulHalf *big.Int

	polypoolQ    *ring.Poly
	polypoolQMul *ring.Poly
}

func NewFastBasisExtender(ringP, ringQ, ringQMul, ringR *ring.Ring) (conv *FastBasisExtender) {
	conv = new(FastBasisExtender)
	conv.ringP = ringP
	conv.ringQ = ringQ
	conv.ringQMul = ringQMul
	conv.ringR = ringR
	conv.baseconverter = mkrlwe.NewFastBasisExtender(ringQ, ringQMul)

	conv.polypoolQ = ringQ.NewPoly()
	conv.polypoolQMul = ringQMul.NewPoly()

	conv.QMulHalf = big.NewInt(1)
	conv.QMulHalf.Rsh(conv.ringQMul.ModulusBigint, 1)

	return conv
}

// assume input polyQP is in NTTForm
func (conv *FastBasisExtender) modUpQPtoRP(polyQP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := len(conv.ringQMul.Modulus) - 1

	conv.ringQ.InvNTTLvl(levelQ, polyQP.Q, conv.polypoolQ)
	conv.baseconverter.ModUpQtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul)
	conv.ringQMul.NTTLvl(levelQMul, conv.polypoolQMul, conv.polypoolQMul)

	polyRP.P.Copy(polyQP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(polyQP.Q.Coeffs[i], polyRP.Q.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(conv.polypoolQMul.Coeffs[i], polyRP.Q.Coeffs[i+levelQ+1])
	}
}

// assume input polyQMulP is in NTTForm
func (conv *FastBasisExtender) modUpQMulPtoRP(polyQMulP, polyRP rlwe.PolyQP) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := len(conv.ringQMul.Modulus) - 1

	conv.ringQMul.InvNTTLvl(levelQMul, polyQMulP.Q, conv.polypoolQMul)
	conv.baseconverter.ModUpPtoQ(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul)
	conv.ringQ.NTTLvl(levelQ, conv.polypoolQ, conv.polypoolQ)

	polyRP.P.Copy(polyQMulP.P)

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], polyRP.Q.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyQMulP.Q.Coeffs[i], polyRP.Q.Coeffs[i+levelQ+1])
	}
}

// Assume input polyQ is in InvNTT Form and output is in NTT
func (conv *FastBasisExtender) ModUpQtoRAndNTT(polyQ, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := len(conv.ringQMul.Modulus) - 1
	levelR := len(conv.ringR.Modulus) - 1

	conv.baseconverter.ModUpQtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul)

	for i := 0; i < levelQ+1; i++ {
		copy(polyQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(conv.polypoolQMul.Coeffs[i], polyR.Coeffs[i+levelQ+1])
	}

	conv.ringR.NTTLvl(levelR, polyR, polyR)
}

func (conv *FastBasisExtender) Quantize(polyR *ring.Poly, t uint64, polyQ *ring.Poly) {
	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := len(conv.ringQMul.Modulus) - 1

	if levelQ != levelQMul {
		panic("cannot ModDownRtoQ: levelQ & levelQMul are different")
	}

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], conv.polypoolQ.Coeffs[i])
	}

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[levelQ+1+i], conv.polypoolQMul.Coeffs[i])
	}

	// Aplies inverse NTT
	conv.ringQ.InvNTTLazyLvl(levelQ, conv.polypoolQ, conv.polypoolQ)
	conv.ringQMul.InvNTTLazyLvl(levelQMul, conv.polypoolQMul, conv.polypoolQMul)

	// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
	conv.baseconverter.ModDownQPtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, conv.polypoolQMul)

	// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
	conv.ringQMul.AddScalarBigint(conv.polypoolQMul, conv.QMulHalf, conv.polypoolQMul)
	conv.baseconverter.ModUpPtoQ(levelQMul, levelQ, conv.polypoolQMul, conv.polypoolQ)
	conv.ringQ.SubScalarBigint(conv.polypoolQ, conv.QMulHalf, conv.polypoolQ)

	// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
	conv.ringQ.MulScalar(conv.polypoolQ, t, polyQ)

}

func (conv *FastBasisExtender) GadgetTransform(swkQP, swkQMulP, swkRP *mkrlwe.SwitchingKey) {
	beta := len(swkQP.Value)
	alpha := len(conv.ringQ.Modulus) / beta

	if beta != len(swkQMulP.Value) {
		panic("cannot GadgetTransform: swkQP & swkQMulP don't have the same dimension")
	}

	if 2*beta != len(swkRP.Value) {
		panic("cannot GadgetTransform: swkRP have wrong dimension")
	}

	Q := conv.ringQ.ModulusBigint
	QMul := conv.ringQMul.ModulusBigint

	InvQ := big.NewInt(0)
	InvQ.ModInverse(Q, QMul)

	InvQMul := big.NewInt(0)
	InvQMul.ModInverse(QMul, Q)

	//transform Q part
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		InvQMulModQi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(conv.ringQ.Modulus[i*alpha+j])))
		}
		InvQMulModQi.Mul(InvQMulModQi, InvQMul)
		InvQMulModQi.Mod(InvQMulModQi, Qi)

		conv.modUpQPtoRP(swkQP.Value[i], swkRP.Value[i])
		conv.ringR.MulScalarBigint(swkRP.Value[i].Q, InvQMulModQi, swkRP.Value[i].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i].P, InvQMulModQi, swkRP.Value[i].P)
	}

	//transform QMul part
	for i := 0; i < beta; i++ {
		QMuli := big.NewInt(1)
		InvQModQMuli := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			QMuli.Mul(QMuli, big.NewInt(int64(conv.ringQMul.Modulus[i*alpha+j])))
		}
		InvQModQMuli.Mul(InvQModQMuli, InvQ)
		InvQModQMuli.Mod(InvQModQMuli, QMuli)

		conv.modUpQMulPtoRP(swkQMulP.Value[i], swkRP.Value[i+beta])
		conv.ringR.MulScalarBigint(swkRP.Value[i+beta].Q, InvQModQMuli, swkRP.Value[i+beta].Q)
		conv.ringP.MulScalarBigint(swkRP.Value[i+beta].P, InvQModQMuli, swkRP.Value[i+beta].P)
	}

}
