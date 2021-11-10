package mkbfv

//import "github.com/ldsec/lattigo/v2/ring"

import "mk-lattigo/mkrlwe"

import "math/big"

type KeyGenerator struct {
	params   Parameters
	keygenQP *mkrlwe.KeyGenerator
	keygenRP *mkrlwe.KeyGenerator
	baseconv *FastBasisExtender

	QInv    *mkrlwe.SwitchingKey
	QMulInv *mkrlwe.SwitchingKey
}

func NewKeyGenerator(params Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.params = params
	keygen.keygenQP = mkrlwe.NewKeyGenerator(params.paramsQP)
	keygen.keygenRP = mkrlwe.NewKeyGenerator(params.paramsRP)
	keygen.baseconv = NewFastBasisExtender(params)

	//gen QInv and QMulInv
	ringR := params.RingR()
	ringP := params.RingP()
	keygen.QInv = mkrlwe.NewSwitchingKey(params.paramsRP)
	keygen.QMulInv = mkrlwe.NewSwitchingKey(params.paramsRP)
	ks := mkrlwe.NewKeySwitcher(params.paramsRP)

	beta := params.paramsRP.Beta(params.RCount() - 1)
	Q := params.RingQ().ModulusBigint
	QMul := params.RingQMul().ModulusBigint
	tmpBigInt := big.NewInt(0)
	tmpPolyR := ringR.NewPoly()
	tmpPolyR.IsNTT = true

	tmpBigInt.ModInverse(Q, QMul)
	ringR.AddScalarBigint(tmpPolyR, tmpBigInt, tmpPolyR)
	ks.Decompose(params.RCount()-1, tmpPolyR, keygen.QInv)
	ringR.SubScalarBigint(tmpPolyR, tmpBigInt, tmpPolyR)

	for i := 0; i < beta; i++ {
		ringR.MulScalarBigint(keygen.QInv.Value[i].Q, QMul, keygen.QInv.Value[i].Q)
		ringP.MulScalarBigint(keygen.QInv.Value[i].P, QMul, keygen.QInv.Value[i].P)
	}

	tmpBigInt.ModInverse(QMul, Q)
	ringR.AddScalarBigint(tmpPolyR, tmpBigInt, tmpPolyR)
	ks.Decompose(params.RCount()-1, tmpPolyR, keygen.QMulInv)
	ringR.SubScalarBigint(tmpPolyR, tmpBigInt, tmpPolyR)

	for i := 0; i < beta; i++ {
		ringR.MulScalarBigint(keygen.QMulInv.Value[i].Q, Q, keygen.QMulInv.Value[i].Q)
		ringP.MulScalarBigint(keygen.QMulInv.Value[i].P, Q, keygen.QMulInv.Value[i].P)
	}

	return keygen
}

func (keygen *KeyGenerator) GenSecretKey(id string) *mkrlwe.SecretKey {
	return keygen.keygenQP.GenSecretKey(id)
}

func (keygen *KeyGenerator) GenPublicKey(sk *mkrlwe.SecretKey) *mkrlwe.PublicKey {
	return keygen.keygenQP.GenPublicKey(sk)
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) GenKeyPair(id string) (sk *mkrlwe.SecretKey, pk *mkrlwe.PublicKey) {
	sk = keygen.GenSecretKey(id)
	return sk, keygen.GenPublicKey(sk)
}

func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *mkrlwe.SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params
	paramsRP := params.paramsRP

	skR := mkrlwe.NewSecretKey(paramsRP, sk.ID)
	rR := mkrlwe.NewSecretKey(paramsRP, sk.ID)

	skR.Value.P.Copy(sk.Value.P)
	rR.Value.P.Copy(rR.Value.P)
	keygen.baseconv.ModUpQtoR(sk.Value.Q, skR.Value.Q)
	keygen.baseconv.ModUpQtoR(r.Value.Q, rR.Value.Q)

	rlk = keygen.keygenRP.GenRelinearizationKey(skR, rR)
	levelR := params.RCount() - 1
	levelP := params.PCount() - 1
	beta := paramsRP.Beta(params.RCount() - 1)
	ringRP := params.RingRP()

	e := ringRP.NewPoly()

	for i := 0; i < beta; i++ {
		keygen.keygenRP.GenGaussianError(e)
		ringRP.MFormLvl(levelR, levelP, e, e)
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QInv.Value[i], rlk.Value[0].Value[i])
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QMulInv.Value[i], rlk.Value[0].Value[i])

		keygen.keygenRP.GenGaussianError(e)
		ringRP.MFormLvl(levelR, levelP, e, e)
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QInv.Value[i], rlk.Value[1].Value[i])
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QMulInv.Value[i], rlk.Value[1].Value[i])

		keygen.keygenRP.GenGaussianError(e)
		ringRP.MFormLvl(levelR, levelP, e, e)
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QInv.Value[i], rlk.Value[2].Value[i])
		ringRP.MulCoeffsMontgomeryAndAddLvl(levelR, levelP, e, keygen.QMulInv.Value[i], rlk.Value[2].Value[i])

	}

	return rlk
}
