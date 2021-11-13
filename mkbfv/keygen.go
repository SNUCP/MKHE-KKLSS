package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

import "mk-lattigo/mkrlwe"

type KeyGenerator struct {
	params      Parameters
	keygenQP    *mkrlwe.KeyGenerator
	keygenQMulP *mkrlwe.KeyGenerator
	baseconv    *FastBasisExtender

	polypoolQ    *ring.Poly
	polypoolQMul *ring.Poly
}

func NewKeyGenerator(params Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.params = params
	keygen.keygenQP = mkrlwe.NewKeyGenerator(params.paramsQP)
	keygen.keygenQMulP = mkrlwe.NewKeyGenerator(params.paramsQMulP)
	keygen.baseconv = NewFastBasisExtender(params.RingP(), params.RingQ(), params.RingQMul(), params.RingR())

	keygen.polypoolQ = params.RingQ().NewPoly()
	keygen.polypoolQMul = params.RingQMul().NewPoly()

	return keygen
}

// genSecretKeyFromSampler generates a new SecretKey sampled from the provided Sampler.
// output SecretKey is in MForm
func (keygen *KeyGenerator) genSecretKeyFromSampler(sampler ring.Sampler, id string) *SecretKey {
	ringQP := keygen.params.RingQP()
	ringQMulP := keygen.params.RingQMulP()
	ringRP := keygen.params.RingRP()
	params := keygen.params

	sk := new(SecretKey)
	sk.SecretKey = mkrlwe.NewSecretKey(keygen.params.paramsQP, id)
	sk.ValueQP = sk.SecretKey.Value
	sk.ValueQMulP = ringQMulP.NewPoly()
	sk.ValueRP = ringRP.NewPoly()
	sk.ID = id

	//gen sk.ValueQP
	levelQ, levelQMul, levelP := keygen.params.QCount()-1, keygen.params.QMulCount()-1, keygen.params.PCount()-1

	sampler.Read(sk.ValueQP.Q)
	ringQP.ExtendBasisSmallNormAndCenter(sk.ValueQP.Q, levelP, nil, sk.ValueQP.P)

	ringQP.NTTLvl(levelQ, levelP, sk.ValueQP, sk.ValueQP)
	ringQP.MFormLvl(levelQ, levelP, sk.ValueQP, sk.ValueQP)

	//gen sk.ValueQMulP
	sk.ValueQMulP.P.Copy(sk.ValueQP.P)
	params.RingQ().InvMForm(sk.ValueQP.Q, keygen.polypoolQ)
	params.RingQ().InvNTT(keygen.polypoolQ, keygen.polypoolQ)
	keygen.baseconv.baseconverter.ModUpQtoP(levelQ, levelQMul, keygen.polypoolQ, keygen.polypoolQMul)
	params.RingQMul().NTT(keygen.polypoolQMul, sk.ValueQMulP.Q)
	params.RingQMul().MForm(sk.ValueQMulP.Q, sk.ValueQMulP.Q)

	//gen sk.ValueRP
	sk.ValueRP.P.Copy(sk.ValueQP.P)
	for i := 0; i < levelQ; i++ {
		copy(sk.ValueRP.Q.Coeffs[i], sk.ValueQP.Q.Coeffs[i])
	}
	for i := 0; i < levelQMul; i++ {
		copy(sk.ValueRP.Q.Coeffs[i+levelQ+1], sk.ValueQP.Q.Coeffs[i])
	}

	return sk
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {

	ringQMulP := keygen.params.RingQMulP()
	ringRP := keygen.params.RingRP()
	params := keygen.params

	sk = new(SecretKey)
	sk.SecretKey = keygen.keygenQP.GenSecretKey(id)
	sk.ValueQP = sk.SecretKey.Value
	sk.ValueQMulP = ringQMulP.NewPoly()
	sk.ValueRP = ringRP.NewPoly()
	sk.ID = id

	levelQ, levelQMul := keygen.params.QCount()-1, keygen.params.QMulCount()-1

	//gen sk.ValueQMulP
	sk.ValueQMulP.P.Copy(sk.ValueQP.P)
	params.RingQ().InvMForm(sk.ValueQP.Q, keygen.polypoolQ)
	params.RingQ().InvNTT(keygen.polypoolQ, keygen.polypoolQ)
	keygen.baseconv.baseconverter.ModUpQtoP(levelQ, levelQMul, keygen.polypoolQ, keygen.polypoolQMul)
	params.RingQMul().NTT(keygen.polypoolQMul, sk.ValueQMulP.Q)
	params.RingQMul().MForm(sk.ValueQMulP.Q, sk.ValueQMulP.Q)

	//gen sk.ValueRP
	sk.ValueRP.P.Copy(sk.ValueQP.P)
	for i := 0; i < levelQ; i++ {
		copy(sk.ValueRP.Q.Coeffs[i], sk.ValueQP.Q.Coeffs[i])
	}
	for i := 0; i < levelQMul; i++ {
		copy(sk.ValueRP.Q.Coeffs[i+levelQ+1], sk.ValueQMulP.Q.Coeffs[i])
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

	id := sk.ID

	skQP := mkrlwe.NewSecretKey(params.paramsQP, id)
	rQP := mkrlwe.NewSecretKey(params.paramsQP, id)

	skQMulP := mkrlwe.NewSecretKey(params.paramsQMulP, id)
	rQMulP := mkrlwe.NewSecretKey(params.paramsQMulP, id)

	skQP.Value = sk.ValueQP
	rQP.Value = r.ValueQP

	skQMulP.Value = sk.ValueQMulP
	rQMulP.Value = r.ValueQMulP

	// gen rlk in mod QP and QMulP
	rlkQP := keygen.keygenQP.GenRelinearizationKey(skQP, rQP)
	rlkQMulP := keygen.keygenQMulP.GenRelinearizationKey(skQMulP, rQMulP)

	// apply GadgetTransform
	rlk = mkrlwe.NewRelinearizationKey(params.paramsRP, id)
	keygen.baseconv.GadgetTransform(rlkQP.Value[0], rlkQMulP.Value[0], rlk.Value[0])
	keygen.baseconv.GadgetTransform(rlkQP.Value[1], rlkQMulP.Value[1], rlk.Value[1])
	keygen.baseconv.GadgetTransform(rlkQP.Value[2], rlkQMulP.Value[2], rlk.Value[2])

	return rlk
}
