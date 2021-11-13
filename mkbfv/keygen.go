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

func (keygen *KeyGenerator) GenRelinearizationKey(skQP, rQP *mkrlwe.SecretKey) (rlk *mkrlwe.RelinearizationKey) {

	params := keygen.params

	levelQ := params.QCount() - 1
	levelQMul := params.QMulCount() - 1

	id := skQP.ID

	// generate sk, r in mod QMulP
	skQMulP := mkrlwe.NewSecretKey(params.paramsQMulP, id)
	rQMulP := mkrlwe.NewSecretKey(params.paramsQMulP, id)

	skQMulP.Value.P.Copy(skQP.Value.P)
	rQMulP.Value.P.Copy(rQP.Value.P)

	params.RingQ().InvMForm(skQP.Value.Q, keygen.polypoolQ)
	params.RingQ().InvNTT(keygen.polypoolQ, keygen.polypoolQ)
	keygen.baseconv.baseconverter.ModUpQtoP(levelQ, levelQMul, keygen.polypoolQ, keygen.polypoolQMul)
	params.RingQMul().NTT(keygen.polypoolQMul, skQMulP.Value.Q)
	params.RingQMul().MForm(skQMulP.Value.Q, skQMulP.Value.Q)

	params.RingQ().InvMForm(rQP.Value.Q, keygen.polypoolQ)
	params.RingQ().InvNTT(keygen.polypoolQ, keygen.polypoolQ)
	keygen.baseconv.baseconverter.ModUpQtoP(levelQ, levelQMul, keygen.polypoolQ, keygen.polypoolQMul)
	params.RingQMul().NTT(keygen.polypoolQMul, rQMulP.Value.Q)
	params.RingQMul().MForm(rQMulP.Value.Q, rQMulP.Value.Q)

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
