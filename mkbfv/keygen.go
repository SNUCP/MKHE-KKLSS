package mkbfv

import "mk-lattigo/mkrlwe"

import "math/big"

type KeyGenerator struct {
	*mkrlwe.KeyGenerator
	keygenRP *mkrlwe.KeyGenerator
	params   Parameters
	conv     *FastBasisExtender
}

func NewKeyGenerator(params Parameters) *KeyGenerator {
	keygen := new(KeyGenerator)
	keygen.KeyGenerator = mkrlwe.NewKeyGenerator(params.Parameters)
	keygen.keygenRP = mkrlwe.NewKeyGenerator(params.paramsRP)
	keygen.params = params
	keygen.conv = NewFastBasisExtender(params.RingP(), params.RingQ(), params.RingQMul(), params.RingR())

	return keygen
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in  MontgomeryForm
func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *mkrlwe.SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params
	paramsRP := params.paramsRP
	conv := keygen.conv
	id := sk.ID

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)

	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()

	//rlk = (b, d, v)
	rlk = mkrlwe.NewRelinearizationKey(paramsRP, id)

	//set CRS
	a1 := keygen.params.CRS[0]
	a2 := keygen.params.CRS[-1]
	u1 := keygen.params.CRS[-1]
	u2 := keygen.params.CRS[-4]

	//generate vector b = QMul * (-sa + e) in MForm
	b1 := mkrlwe.NewSwitchingKey(params.Parameters)
	b2 := mkrlwe.NewSwitchingKey(params.Parameters)
	b := rlk.Value[0]
	e := ringQP.NewPoly()

	// compute -sa + e
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a1.Value[i], sk.Value, b1.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])
		keygen.GenGaussianError(e)
		ringQP.SubLvl(levelQ, levelP, e, b1.Value[i], b1.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])

		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a2.Value[i], sk.Value, b2.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
		keygen.GenGaussianError(e)
		ringQP.SubLvl(levelQ, levelP, e, b2.Value[i], b2.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
	}
	conv.GadgetTransform(b1, b2, b)

	//generate vector d = QMul * (-ra + sg + e) in MForm
	d := rlk.Value[1]
	d1 := mkrlwe.NewSwitchingKey(params.Parameters)
	d2 := mkrlwe.NewSwitchingKey(params.Parameters)
	keygen.genBFVSwitchingKey(sk, d1, d2)

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a1.Value[i], r.Value, d1.Value[i])
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, a2.Value[i], r.Value, d2.Value[i])
	}
	conv.GadgetTransform(d1, d2, d)

	//generate vector v = -su - rg + e in MForm
	v := rlk.Value[2]
	v1 := mkrlwe.NewSwitchingKey(params.Parameters)
	v2 := mkrlwe.NewSwitchingKey(params.Parameters)
	keygen.genBFVSwitchingKey(r, v1, v2)

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u1.Value[i], sk.Value, v1.Value[i])
		ringQ.NegLvl(levelQ, v1.Value[i].Q, v1.Value[i].Q)
		ringP.NegLvl(levelP, v1.Value[i].P, v1.Value[i].P)

		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u2.Value[i], sk.Value, v2.Value[i])
		ringQ.NegLvl(levelQ, v2.Value[i].Q, v2.Value[i].Q)
		ringP.NegLvl(levelP, v2.Value[i].P, v2.Value[i].P)
	}
	conv.GadgetTransform(v1, v2, v)

	return
}

func (keygen *KeyGenerator) genBFVSwitchingKey(sk *mkrlwe.SecretKey, swk1, swk2 *mkrlwe.SwitchingKey) {
	params := keygen.params

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)
	alpha := params.Alpha()

	ringQ := params.RingQ()
	ringQMul := params.RingQMul()
	ringP := params.RingP()
	ringQP := params.RingQP()

	//g_i = QQMul/Q_i * [(QQMul/Q_i)^-1]Q_i
	QQMul := big.NewInt(1)
	QQMul.Mul(ringQ.ModulusBigint, ringQMul.ModulusBigint)

	e := ringQP.NewPoly()

	//gen swk1
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(ringQ.Modulus[i*alpha+j])))
		}

		Gi := big.NewInt(1)
		Gi.Div(QQMul, Qi)
		Ti := big.NewInt(1)
		Ti.ModInverse(Gi, Qi)
		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, ringP.ModulusBigint)
		Gi.Div(Gi, ringQMul.ModulusBigint)

		ringQP.InvMFormLvl(levelQ, levelP, sk.Value, swk1.Value[i])
		ringQ.MulScalarBigint(swk1.Value[i].Q, Gi, swk1.Value[i].Q)
		ringP.MulScalarBigint(swk1.Value[i].P, Gi, swk1.Value[i].P)

		keygen.GenGaussianError(e)
		ringQP.AddLvl(levelQ, levelP, swk1.Value[i], e, swk1.Value[i])

		ringQP.MFormLvl(levelQ, levelP, swk1.Value[i], swk1.Value[i])
	}

	//gen swk2
	for i := 0; i < beta; i++ {
		Qi := big.NewInt(1)
		for j := 0; j < alpha; j++ {
			Qi.Mul(Qi, big.NewInt(int64(ringQMul.Modulus[i*alpha+j])))
		}

		Gi := big.NewInt(1)
		Gi.Div(QQMul, Qi)
		Ti := big.NewInt(1)
		Ti.ModInverse(Gi, Qi)
		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, ringP.ModulusBigint)
		Gi.Div(Gi, ringQMul.ModulusBigint)

		ringQP.InvMFormLvl(levelQ, levelP, sk.Value, swk2.Value[i])
		ringQ.MulScalarBigint(swk2.Value[i].Q, Gi, swk2.Value[i].Q)
		ringP.MulScalarBigint(swk2.Value[i].P, Gi, swk2.Value[i].P)

		keygen.GenGaussianError(e)
		ringQP.AddLvl(levelQ, levelP, swk2.Value[i], e, swk2.Value[i])

		ringQP.MFormLvl(levelQ, levelP, swk2.Value[i], swk2.Value[i])
	}

}
