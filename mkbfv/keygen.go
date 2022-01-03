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

// GenRelinKey generates a new EvaluationKey that will be used to relinearize constant term of Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in  MontgomeryForm
func (keygen *KeyGenerator) GenConstRelinearizationKey() (rlk *mkrlwe.RelinearizationKey) {
	params := keygen.params
	ringQ := params.RingQ()
	ringP := params.RingP()

	sk := mkrlwe.NewSecretKey(params.Parameters, "0")
	r := mkrlwe.NewSecretKey(params.Parameters, "0")

	ringQ.AddScalar(sk.Value.Q, 1, sk.Value.Q)
	ringP.AddScalar(sk.Value.P, 1, sk.Value.P)

	ringQ.MForm(sk.Value.Q, sk.Value.Q)
	ringP.MForm(sk.Value.P, sk.Value.P)

	ringQ.MForm(r.Value.Q, r.Value.Q)
	ringP.MForm(r.Value.P, r.Value.P)

	rlk = keygen.GenRelinearizationKey(sk, r)
	return
}

func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *mkrlwe.SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQP := params.RingQP()
	ringQ := params.RingQ()
	ringP := params.RingP()

	id := sk.ID

	//rlk = (b, d, v)
	rlk = mkrlwe.NewRelinearizationKey(params.paramsRP, id)
	beta := params.Beta(levelQ)

	//set CRS
	a1 := params.CRS[0]
	a2 := params.CRS[-3]

	u1 := params.CRS[-1]
	u2 := params.CRS[-4]

	//generate vector b = -sa + e in MForm
	b := rlk.Value[0]
	b1 := mkrlwe.NewSwitchingKey(params.Parameters)
	b2 := mkrlwe.NewSwitchingKey(params.Parameters)
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
	keygen.baseconv.GadgetTransform(b1, b2, b)

	//generate vector d = -ra + sg + e in MForm
	d := rlk.Value[1]
	d1 := mkrlwe.NewSwitchingKey(params.Parameters)
	d2 := mkrlwe.NewSwitchingKey(params.Parameters)

	keygen.GenBFVSwitchingKey(sk, d1, d2)
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a1.Value[i], r.Value, d1.Value[i])
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a2.Value[i], r.Value, d2.Value[i])
	}
	keygen.baseconv.GadgetTransform(d1, d2, d)

	//generate vector v = -su - rg + e in MForm
	v := rlk.Value[2]
	v1 := mkrlwe.NewSwitchingKey(params.Parameters)
	v2 := mkrlwe.NewSwitchingKey(params.Parameters)

	keygen.GenBFVSwitchingKey(r, v1, v2)

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u1.Value[i], sk.Value, v1.Value[i])
		ringQ.NegLvl(levelQ, v1.Value[i].Q, v1.Value[i].Q)
		ringP.NegLvl(levelP, v1.Value[i].P, v1.Value[i].P)

		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u2.Value[i], sk.Value, v2.Value[i])
		ringQ.NegLvl(levelQ, v2.Value[i].Q, v2.Value[i].Q)
		ringP.NegLvl(levelP, v2.Value[i].P, v2.Value[i].P)
	}
	keygen.baseconv.GadgetTransform(v1, v2, v)

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
		Gi.Div(Q, Qi)

		Ti := big.NewInt(1)
		Ti.Div(QQMul, Qi)
		Ti.ModInverse(Ti, Qi)

		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, P)

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
