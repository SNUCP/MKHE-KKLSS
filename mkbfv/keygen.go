package mkbfv

import "github.com/ldsec/lattigo/v2/ring"

import "mk-lattigo/mkrlwe"

type KeyGenerator struct {
	params    Parameters
	keygenQP  *mkrlwe.KeyGenerator
	keygenQ1P *mkrlwe.KeyGenerator
	keygenRP  *mkrlwe.KeyGenerator
	baseconv  *FastBasisExtender

	polypoolQ  *ring.Poly
	polypoolQ1 *ring.Poly
}

func NewKeyGenerator(params Parameters) (keygen *KeyGenerator) {
	keygen = new(KeyGenerator)
	keygen.params = params
	keygen.keygenQP = mkrlwe.NewKeyGenerator(params.paramsQP)
	keygen.keygenQ1P = mkrlwe.NewKeyGenerator(params.paramsQ1P)
	keygen.keygenRP = mkrlwe.NewKeyGenerator(params.paramsRP)
	keygen.baseconv = NewFastBasisExtender(
		params.RingP(), params.RingQ(),
		params.RingQ1(), params.RingR(),
	)

	keygen.polypoolQ = params.RingQ().NewPoly()
	keygen.polypoolQ1 = params.RingQ1().NewPoly()

	return keygen
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {

	params := keygen.params

	sk = new(SecretKey)
	sk.SecretKey = mkrlwe.NewSecretKey(params.paramsQP, id)
	sk.ValueQP = sk.SecretKey.Value
	sk.ValueQ1P = params.RingQ1P().NewPoly()
	sk.ValueRP = keygen.keygenRP.GenSecretKey(id).Value
	sk.ID = id

	levelQ := keygen.params.QCount() - 1
	levelQ1 := levelQ

	sk.ValueQP.P.Copy(sk.ValueRP.P)
	for i := 0; i < levelQ+1; i++ {
		copy(sk.ValueQP.Q.Coeffs[i], sk.ValueRP.Q.Coeffs[i])
	}

	sk.ValueQ1P.P.Copy(sk.ValueRP.P)
	for i := 0; i < levelQ1+1; i++ {
		copy(sk.ValueQ1P.Q.Coeffs[i], sk.ValueRP.Q.Coeffs[i+levelQ+1])
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

	skQ1P := mkrlwe.NewSecretKey(params.paramsQ1P, id)
	rQ1P := mkrlwe.NewSecretKey(params.paramsQ1P, id)

	skQP.Value.Copy(sk.ValueQP)
	rQP.Value.Copy(r.ValueQP)

	skQ1P.Value.Copy(sk.ValueQ1P)
	rQ1P.Value.Copy(r.ValueQ1P)

	// gen rlk in mod QP and QMulP
	rlkQP := keygen.keygenQP.GenRelinearizationKey(skQP, rQP)
	rlkQ1P := keygen.keygenQ1P.GenRelinearizationKey(skQ1P, rQ1P)

	// apply GadgetTransform
	rlk = mkrlwe.NewRelinearizationKey(params.paramsRP, id)
	keygen.baseconv.GadgetTransform(rlkQP.Value[0], rlkQ1P.Value[0], rlk.Value[0])
	keygen.baseconv.GadgetTransform(rlkQP.Value[1], rlkQ1P.Value[1], rlk.Value[1])
	keygen.baseconv.GadgetTransform(rlkQP.Value[2], rlkQ1P.Value[2], rlk.Value[2])

	return rlk
}
