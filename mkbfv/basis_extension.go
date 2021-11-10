package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "mk-lattigo/mkrlwe"

type FastBasisExtender struct {
	params             Parameters
	baseconverterQQMul *mkrlwe.FastBasisExtender
	polypoolQMul       *ring.Poly
}

func NewFastBasisExtender(params Parameters) (baseconverter *FastBasisExtender) {
	baseconverter = new(FastBasisExtender)

	baseconverter.params = params
	baseconverter.baseconverterQQMul = mkrlwe.NewFastBasisExtender(params.RingQ(), params.RingQMul())
	baseconverter.polypoolQMul = params.RingQMul().NewPoly()

	return baseconverter
}

func (conv *FastBasisExtender) ModUpQtoR(polyQ, polyR *ring.Poly) {

	levelQ := conv.params.QCount() - 1
	levelQMul := conv.params.QMulCount() - 1

	conv.baseconverterQQMul.ModUpQtoP(levelQ, levelQMul, polyQ, conv.polypoolQMul)

	for i := 0; i < levelQ+1; i++ {
		copy(polyQ.Coeffs[i], polyR.Coeffs[i])
	}

	for i := 0; i < levelQMul+1; i++ {
		copy(conv.polypoolQMul.Coeffs[i+levelQ+1], polyR.Coeffs[i+levelQ+1])
	}
}
