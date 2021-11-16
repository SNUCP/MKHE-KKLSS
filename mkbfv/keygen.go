package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

import "mk-lattigo/mkrlwe"

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
	keygen.baseconv = NewFastBasisExtender(params.RingP(), params.RingQ(), params.RingQMul(), params.RingR())

	keygen.polypoolQ = params.RingQ().NewPoly()
	keygen.polypoolQMul = params.RingQMul().NewPoly()

	return keygen
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {

	ringQMulP := keygen.params.RingQMulP()
	params := keygen.params

	sk = new(SecretKey)
	sk.SecretKey = mkrlwe.NewSecretKey(params.paramsQP, id)
	sk.ValueQP = sk.SecretKey.Value
	sk.ValueQMulP = ringQMulP.NewPoly()
	//sk.ValueRP = keygen.params.RingRP().NewPoly()
	sk.ValueRP = keygen.keygenRP.GenSecretKey(id).Value
	sk.ID = id

	levelQ := keygen.params.QCount() - 1
	levelQMul := keygen.params.QMulCount() - 1

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

	id := sk.ID

	/*
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
	*/
	skRP := mkrlwe.NewSecretKey(params.paramsRP, id)
	rRP := mkrlwe.NewSecretKey(params.paramsRP, id)

	skRP.Value.Copy(sk.ValueRP)
	rRP.Value.Copy(r.ValueRP)

	rlk = keygen.keygenRP.GenRelinearizationKey(skRP, rRP)
	return rlk
}
