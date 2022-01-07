package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "mk-lattigo/mkrlwe"

type FastBasisExtender struct {
	ringP    *ring.Ring
	ringQ    *ring.Ring
	ringQMul *ring.Ring
	ringR    *ring.Ring

	convQQMul *mkrlwe.FastBasisExtender

	polypoolQ    *ring.Poly
	polypoolQMul *ring.Poly
	polypoolR    *ring.Poly
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
	conv.polypoolR = ringR.NewPoly()

	conv.mFormQMul = ringQ.NewPoly()
	ringQ.AddScalarBigint(conv.mFormQMul, ringQMul.ModulusBigint, conv.mFormQMul)
	ringQ.MForm(conv.mFormQMul, conv.mFormQMul)

	return conv
}

// assume input and output are in InvNTTForm
func (conv *FastBasisExtender) ModUpQtoR(polyQ, polyR *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.convQQMul.ModUpQtoP(levelQ, levelQMul, polyQ, conv.polypoolQMul)

	for i := 0; i < levelQ+1; i++ {
		copy(polyR.Coeffs[i], polyQ.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(polyR.Coeffs[i+levelQ+1], conv.polypoolQMul.Coeffs[i])
	}
}

// assume input and output are in InvNTTForm
func (conv *FastBasisExtender) Quantize(polyR, polyQ *ring.Poly, t uint64) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringR.MulScalar(polyR, t, conv.polypoolR)
	conv.ringR.InvNTT(conv.polypoolR, conv.polypoolR)

	for i := 0; i < levelQ+1; i++ {
		copy(conv.polypoolQ.Coeffs[i], conv.polypoolR.Coeffs[i])
		copy(conv.polypoolQMul.Coeffs[i], conv.polypoolR.Coeffs[i+levelQ+1])
	}

	conv.convQQMul.ModDownQPtoQ(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, polyQ)
}

// assume input polyQ is in InvNTTForm
func (conv *FastBasisExtender) Rescale(polyQ *ring.Poly, polyQOut *ring.Poly) {

	levelQ := len(conv.ringQ.Modulus) - 1
	levelQMul := levelQ

	conv.ringQ.MulCoeffsMontgomery(polyQ, conv.mFormQMul, conv.polypoolQ)
	conv.ringQMul.MulScalar(conv.polypoolQMul, 0, conv.polypoolQMul)
	conv.convQQMul.ModDownQPtoP(levelQ, levelQMul, conv.polypoolQ, conv.polypoolQMul, conv.polypoolQMul)
	conv.convQQMul.ModUpPtoQ(levelQMul, levelQ, conv.polypoolQMul, polyQOut)
}
