package mkbfv

import "github.com/ldsec/lattigo/v2/ring"
import "mk-lattigo/mkrlwe"
import "math/big"

type KeyGenerator struct {
	params      Parameters
	keygenQP    *mkrlwe.KeyGenerator
	keygenQMulP *mkrlwe.KeyGenerator
	keygenRP    *mkrlwe.KeyGenerator
	baseconv    *FastBasisExtender

	polypoolQ    *ring.Poly
	polypoolQMul *ring.Poly
}

func NewKeyGenerator(params Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.params = params
	keygen.keygenQP = mkrlwe.NewKeyGenerator(params.paramsQP)
	keygen.keygenQMulP = mkrlwe.NewKeyGenerator(params.paramsQMulP)
	keygen.keygenRP = mkrlwe.NewKeyGenerator(params.paramsRP)
	keygen.baseconv = NewFastBasisExtender(
		params.RingP(), params.RingQ(),
		params.RingQMul(), params.RingR(),
	)

	keygen.polypoolQ = params.RingQ().NewPoly()
	keygen.polypoolQMul = params.RingQMul().NewPoly()

	return keygen
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {

	params := keygen.params

	sk = new(SecretKey)
	sk.SecretKey = mkrlwe.NewSecretKey(params.paramsQP, id)
	sk.ValueQP = sk.SecretKey.Value
	sk.ValueQMulP = params.RingQMulP().NewPoly()
	sk.ValueRP = keygen.keygenRP.GenSecretKey(id).Value
	sk.ID = id

	levelQ := keygen.params.QCount() - 1
	levelQMul := levelQ

	sk.ValueQP.P.Copy(sk.ValueRP.P)
	for i := 0; i < levelQ+1; i++ {
		copy(sk.ValueQP.Q.Coeffs[i], sk.ValueRP.Q.Coeffs[i])
	}

	sk.ValueQMulP.P.Copy(sk.ValueRP.P)
	for i := 0; i < levelQMul+1; i++ {
		copy(sk.ValueQMulP.Q.Coeffs[i], sk.ValueRP.Q.Coeffs[i+levelQ+1])
	}

	return sk
}

func (keygen *KeyGenerator) GenPublicKey(sk *SecretKey) *mkrlwe.PublicKey {
	return keygen.keygenQP.GenPublicKey(sk.SecretKey)
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) GenKeyPair(id string) (sk *SecretKey, pk *mkrlwe.PublicKey) {
	sk = keygen.GenSecretKey(id)
	return sk, keygen.GenPublicKey(sk)
}

func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQP := params.RingQP()
	ringQMulP := params.RingQMulP()
	ringQ := params.RingQ()
	ringQMul := params.RingQMul()
	ringP := params.RingP()

	id := sk.ID

	//rlk = (b, d, v)
	rlk = mkrlwe.NewRelinearizationKey(params.paramsRP, id)
	beta := params.paramsQP.Beta(levelQ)

	//set CRS
	a1 := params.paramsQP.CRS[0]
	a2 := params.paramsQMulP.CRS[0]

	u1 := params.paramsQP.CRS[-1]
	u2 := params.paramsQMulP.CRS[-1]

	//generate vector b = -sa + e in MForm
	b := rlk.Value[0]
	b1 := mkrlwe.NewSwitchingKey(params.paramsQP)
	b2 := mkrlwe.NewSwitchingKey(params.paramsQMulP)
	e1 := params.RingQP().NewPoly()
	e2 := params.RingQMulP().NewPoly()

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a1.Value[i], sk.ValueQP, b1.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])
		keygen.keygenQP.GenGaussianError(e1)
		ringQP.SubLvl(levelQ, levelP, e1, b1.Value[i], b1.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])

		ringQMulP.MulCoeffsMontgomeryLvl(levelQ, levelP, a2.Value[i], sk.ValueQMulP, b2.Value[i])
		ringQMulP.InvMFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
		keygen.keygenQMulP.GenGaussianError(e2)
		ringQMulP.SubLvl(levelQ, levelP, e2, b2.Value[i], b2.Value[i])
		ringQMulP.MFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
	}
	keygen.baseconv.GadgetTransform(b1, b2, b)

	//generate vector d = -ra + sg + e in MForm
	d := rlk.Value[1]
	d1 := mkrlwe.NewSwitchingKey(params.paramsQP)
	d2 := mkrlwe.NewSwitchingKey(params.paramsQMulP)

	keygen.GenBFVSwitchingKey(sk, d1, d2)
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a1.Value[i], r.ValueQP, d1.Value[i])
		ringQMulP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a2.Value[i], r.ValueQMulP, d2.Value[i])
	}
	keygen.baseconv.GadgetTransform(d1, d2, d)

	//generate vector v = -su - rg + e in MForm
	v := rlk.Value[2]
	v1 := mkrlwe.NewSwitchingKey(params.paramsQP)
	v2 := mkrlwe.NewSwitchingKey(params.paramsQMulP)

	keygen.GenBFVSwitchingKey(r, v1, v2)

	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u1.Value[i], sk.ValueQP, v1.Value[i])
		ringQ.NegLvl(levelQ, v1.Value[i].Q, v1.Value[i].Q)
		ringP.NegLvl(levelP, v1.Value[i].P, v1.Value[i].P)

		ringQMulP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u2.Value[i], sk.ValueQMulP, v2.Value[i])
		ringQMul.NegLvl(levelQ, v2.Value[i].Q, v2.Value[i].Q)
		ringP.NegLvl(levelP, v2.Value[i].P, v2.Value[i].P)
	}
	keygen.baseconv.GadgetTransform(v1, v2, v)

	return
}

//For an input secretkey s, gen gs + e in MForm
func (keygen *KeyGenerator) GenBFVSwitchingKey(sk *SecretKey, swk1, swk2 *mkrlwe.SwitchingKey) {
	params := keygen.params
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	alpha := params.paramsQP.Alpha()
	beta := params.paramsQP.Beta(levelQ)

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

		params.RingQP().InvMFormLvl(levelQ, levelP, sk.ValueQP, swk1.Value[i])
		params.RingQ().MulScalarBigint(swk1.Value[i].Q, Gi, swk1.Value[i].Q)
		params.RingP().MulScalarBigint(swk1.Value[i].P, Gi, swk1.Value[i].P)

		e := params.RingQP().NewPoly()
		keygen.keygenQP.GenGaussianError(e)

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
		Gi.Div(QMul, QMuli)

		Ti := big.NewInt(1)
		Ti.Div(QQMul, QMuli)
		Ti.ModInverse(Ti, QMuli)

		Gi.Mul(Gi, Ti)
		Gi.Mul(Gi, P)

		params.RingQMulP().InvMFormLvl(levelQ, levelP, sk.ValueQMulP, swk2.Value[i])
		params.RingQMul().MulScalarBigint(swk2.Value[i].Q, Gi, swk2.Value[i].Q)
		params.RingP().MulScalarBigint(swk2.Value[i].P, Gi, swk2.Value[i].P)

		e := params.RingQMulP().NewPoly()
		keygen.keygenQMulP.GenGaussianError(e)

		params.RingQMulP().AddLvl(levelQ, levelP, swk2.Value[i], e, swk2.Value[i])
		params.RingQMulP().MFormLvl(levelQ, levelP, swk2.Value[i], swk2.Value[i])
	}

}

func (keygen *KeyGenerator) GenRotationKey(rotidx int, sk *SecretKey) (rk *mkrlwe.RotationKey) {
	params := keygen.params

	id := sk.ID

	skQP := mkrlwe.NewSecretKey(params.paramsQP, id)
	skQP.Value.Copy(sk.ValueQP)

	rk = keygen.keygenQP.GenRotationKey(rotidx, skQP)
	return rk
}

func (keygen *KeyGenerator) GenDefaultRotationKeys(sk *SecretKey, rtkSet *mkrlwe.RotationKeySet) {

	params := keygen.params

	id := sk.ID

	skQP := mkrlwe.NewSecretKey(params.paramsQP, id)
	skQP.Value.Copy(sk.ValueQP)

	keygen.keygenQP.GenDefaultRotationKeys(skQP, rtkSet)
}

func (keygen *KeyGenerator) GenConjugationKey(sk *SecretKey) (cjk *mkrlwe.ConjugationKey) {
	params := keygen.params

	id := sk.ID

	skQP := mkrlwe.NewSecretKey(params.paramsQP, id)
	skQP.Value.Copy(sk.ValueQP)

	cjk = keygen.keygenQP.GenConjugationKey(skQP)
	return cjk
}
