package mkbfv

import "mk-lattigo/mkrlwe"

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

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {

	params := keygen.params

	sk = new(SecretKey)
	sk.SecretKeyRP = keygen.keygenRP.GenSecretKey(id)
	sk.SecretKey = mkrlwe.NewSecretKey(params.Parameters, id)
	sk.ID = id

	levelQ := keygen.params.QCount() - 1

	sk.SecretKey.Value.P.Copy(sk.SecretKeyRP.Value.P)
	for i := 0; i < levelQ+1; i++ {
		copy(sk.SecretKey.Value.Q.Coeffs[i], sk.SecretKeyRP.Value.Q.Coeffs[i])
	}

	return sk
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) GenKeyPair(id string) (sk *SecretKey, pk *mkrlwe.PublicKey) {
	sk = keygen.GenSecretKey(id)
	return sk, keygen.GenPublicKey(sk.SecretKey)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in  MontgomeryForm
func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params
	paramsRP := params.paramsRP
	conv := keygen.conv
	swkRP := mkrlwe.NewSwitchingKey(paramsRP)
	id := sk.ID

	keygen.keygenRP.GenSwitchingKey(sk.SecretKeyRP, swkRP)

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
	a2 := keygen.params.CRS[-3]
	u1 := keygen.params.CRS[-1]
	u2 := keygen.params.CRS[-4]

	//generate vector b = QMul * (-sa + e) in MForm
	b1 := mkrlwe.NewSwitchingKey(params.Parameters)
	b2 := mkrlwe.NewSwitchingKey(params.Parameters)
	b := rlk.Value[0]
	tmp := ringQP.NewPoly()

	// compute -sa + e
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a1.Value[i], sk.Value, b1.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])
		keygen.GenGaussianError(tmp)
		ringQP.SubLvl(levelQ, levelP, tmp, b1.Value[i], b1.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b1.Value[i], b1.Value[i])

		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a2.Value[i], sk.Value, b2.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b2.Value[i], b2.Value[i])
		keygen.GenGaussianError(tmp)
		ringQP.SubLvl(levelQ, levelP, tmp, b2.Value[i], b2.Value[i])
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

func (keygen *KeyGenerator) genBFVSwitchingKey(sk *SecretKey, swk1, swk2 *mkrlwe.SwitchingKey) {
	params := keygen.params
	paramsRP := params.paramsRP
	conv := keygen.conv
	swkRP := mkrlwe.NewSwitchingKey(paramsRP)

	keygen.keygenRP.GenSwitchingKey(sk.SecretKeyRP, swkRP)

	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	levelR := paramsRP.QCount() - 1
	beta := params.Beta(levelQ)

	ringQP := params.RingQP()
	ringRP := params.RingRP()

	for i := 0; i < 2*beta; i++ {
		ringRP.InvMFormLvl(levelR, levelP, swkRP.Value[i], swkRP.Value[i])
		ringRP.InvNTTLvl(levelR, levelP, swkRP.Value[i], swkRP.Value[i])
	}

	// ModDown by QMul
	for i := 0; i < beta; i++ {
		conv.ModDownRPtoQP(swkRP.Value[i], swk1.Value[i])
		conv.ModDownRPtoQP(swkRP.Value[i+beta], swk2.Value[i])
	}

	// add error , apply NTT and MForm
	e := ringQP.NewPoly()
	for i := 0; i < beta; i++ {
		keygen.GenGaussianError(e)
		ringQP.InvNTTLvl(levelQ, levelP, e, e)
		ringQP.AddLvl(levelQ, levelP, swk1.Value[i], e, swk1.Value[i])
		ringQP.NTTLvl(levelQ, levelP, swk1.Value[i], swk1.Value[i])
		ringQP.MFormLvl(levelQ, levelP, swk1.Value[i], swk1.Value[i])

		keygen.GenGaussianError(e)
		ringQP.InvNTTLvl(levelQ, levelP, e, e)
		ringQP.AddLvl(levelQ, levelP, swk2.Value[i], e, swk2.Value[i])
		ringQP.NTTLvl(levelQ, levelP, swk2.Value[i], swk2.Value[i])
		ringQP.MFormLvl(levelQ, levelP, swk2.Value[i], swk2.Value[i])
	}

}
