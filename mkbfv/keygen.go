package mkbfv

import "mk-lattigo/mkrlwe"
import "math/big"

type KeyGenerator struct {
	*mkrlwe.KeyGenerator
	params   Parameters
	baseconv *FastBasisExtender
}

func NewKeyGenerator(params Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.params = params
	keygen.KeyGenerator = mkrlwe.NewKeyGenerator(params.Parameters)
	keygen.baseconv = NewFastBasisExtender(
		params.RingP(), params.RingQ(),
		params.RingQMul(), params.RingR(),
	)

	return keygen
}

func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *mkrlwe.SecretKey) (rlk *RelinearizationKey) {

	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQP := params.RingQP()
	ringQ := params.RingQ()
	ringP := params.RingP()

	id := sk.ID

	//rlk = (b, d, v)
	rlk = NewRelinearizationKey(params, id)
	beta := params.Beta(levelQ)

	//set CRS
	a1 := params.CRS[0]
	a2 := params.CRS[-3]

	u := params.CRS[-1]

	//generate vector b = -sa + e in MForm
	b1 := rlk.Value[0].Value[0]
	b2 := rlk.Value[1].Value[0]
	e1 := params.RingQP().NewPoly()
	e2 := params.RingQP().NewPoly()

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a1.Value[i], sk.Value, b1.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])
		keygen.GenGaussianError(e1)
		ringQP.SubLvl(levelQ, levelP, e1, b1.Value[i], b1.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])

		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a2.Value[i], sk.Value, b2.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
		keygen.GenGaussianError(e2)
		ringQP.SubLvl(levelQ, levelP, e2, b2.Value[i], b2.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
	}

	//generate vector d = -ra + sg + e in MForm
	d1 := rlk.Value[0].Value[1]
	d2 := rlk.Value[1].Value[1]

	keygen.GenBFVSwitchingKey(sk, d1, d2)
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a1.Value[i], r.Value, d1.Value[i])
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a2.Value[i], r.Value, d2.Value[i])
	}

	//generate vector v = -su - rg + e in MForm
	v := rlk.Value[0].Value[2]

	keygen.GenSwitchingKey(r, v)

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u.Value[i], sk.Value, v.Value[i])
		ringQ.NegLvl(levelQ, v.Value[i].Q, v.Value[i].Q)
		ringP.NegLvl(levelP, v.Value[i].P, v.Value[i].P)

	}

	return
}

//For an input secretkey s, gen gs + e in MForm
func (keygen *KeyGenerator) GenBFVSwitchingKey(sk *mkrlwe.SecretKey, swk1, swk2 *mkrlwe.SwitchingKey) {
	params := keygen.params
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	alpha := params.Alpha()
	beta := params.Beta(levelQ)

	Q := params.RingQ().ModulusBigint
	QMul := params.RingQMul().ModulusBigint
	QQMul := big.NewInt(1)
	QQMul.Mul(Q, QMul)

	P := params.RingP().ModulusBigint

	//gen swk1
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(params.RingQ().Modulus[i*alpha+j])))
		}
		Gi := big.NewInt(1)
		Gi.Div(QQMul, Qi)

		Ti := big.NewInt(1)
		Ti.Div(QQMul, Qi)
		Ti.ModInverse(Ti, Qi)

		Gi.Mul(Gi, params.RingT().ModulusBigint)
		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, P)
		Gi.Div(Gi, QMul)

		params.RingQP().InvMFormLvl(levelQ, levelP, sk.Value, swk1.Value[i])
		params.RingQ().MulScalarBigint(swk1.Value[i].Q, Gi, swk1.Value[i].Q)
		params.RingP().MulScalarBigint(swk1.Value[i].P, Gi, swk1.Value[i].P)

		e := params.RingQP().NewPoly()
		keygen.GenGaussianError(e)

		params.RingQP().AddLvl(levelQ, levelP, swk1.Value[i], e, swk1.Value[i])
		params.RingQP().MFormLvl(levelQ, levelP, swk1.Value[i], swk1.Value[i])
	}

	//gen swk2
	for i := 0; i < beta; i++ {
		QMuli := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			QMuli.Mul(QMuli, big.NewInt(int64(params.RingQMul().Modulus[i*alpha+j])))
		}
		Gi := big.NewInt(1)
		Gi.Div(QQMul, QMuli)

		Ti := big.NewInt(1)
		Ti.Div(QQMul, QMuli)
		Ti.ModInverse(Ti, QMuli)

		Gi.Mul(Gi, params.RingT().ModulusBigint)
		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, P)
		Gi.Div(Gi, QMul)

		params.RingQP().InvMFormLvl(levelQ, levelP, sk.Value, swk2.Value[i])
		params.RingQ().MulScalarBigint(swk2.Value[i].Q, Gi, swk2.Value[i].Q)
		params.RingP().MulScalarBigint(swk2.Value[i].P, Gi, swk2.Value[i].P)

		e := params.RingQP().NewPoly()
		keygen.GenGaussianError(e)

		params.RingQP().AddLvl(levelQ, levelP, swk2.Value[i], e, swk2.Value[i])
		params.RingQP().MFormLvl(levelQ, levelP, swk2.Value[i], swk2.Value[i])
	}

}
