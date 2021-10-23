package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"
import "math"
import "math/big"

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type KeyGenerator struct {
	params           Parameters
	poolQ            *ring.Poly
	poolQP           rlwe.PolyQP
	gaussianSamplerQ *ring.GaussianSampler
	uniformSamplerQ  *ring.UniformSampler
	uniformSamplerP  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) *KeyGenerator {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	keygen := new(KeyGenerator)
	keygen.params = params
	keygen.poolQ = params.RingQ().NewPoly()
	keygen.poolQP = params.RingQP().NewPoly()
	keygen.gaussianSamplerQ = ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma()))
	keygen.uniformSamplerQ = ring.NewUniformSampler(prng, params.RingQ())
	keygen.uniformSamplerP = ring.NewUniformSampler(prng, params.RingP())

	return keygen
}

// genSecretKeyFromSampler generates a new SecretKey sampled from the provided Sampler.
// output SecretKey is in MForm
func (keygen *KeyGenerator) genSecretKeyFromSampler(sampler ring.Sampler, id string) *SecretKey {
	ringQP := keygen.params.RingQP()
	sk := new(SecretKey)
	sk.Value = ringQP.NewPoly()
	sk.ID = id
	levelQ, levelP := keygen.params.QCount()-1, keygen.params.PCount()-1
	sampler.Read(sk.Value.Q)
	ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, nil, sk.Value.P)
	ringQP.NTTLvl(levelQ, levelP, sk.Value, sk.Value)
	ringQP.MFormLvl(levelQ, levelP, sk.Value, sk.Value)
	return sk
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *KeyGenerator) GenSecretKey(id string) (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0/3, id)
}

// GenSecretKey generates a new SecretKey with the error distribution.
func (keygen *KeyGenerator) GenSecretKeyGaussian(id string) (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.gaussianSamplerQ, id)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *KeyGenerator) GenSecretKeyWithDistrib(p float64, id string) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.params.RingQ(), p, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery, id)
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *KeyGenerator) GenSecretKeySparse(hw int, id string) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.params.RingQ(), hw, false)
	return keygen.genSecretKeyFromSampler(ternarySamplerMontgomery, id)
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *KeyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)
	ringQP := keygen.params.RingQP()
	levelQ, levelP := keygen.params.QCount()-1, keygen.params.PCount()-1

	id := sk.ID

	//pk[0] = [-as + e]
	//pk[1] = [a]
	pk = NewPublicKey(keygen.params, id)
	keygen.gaussianSamplerQ.Read(pk.Value[0].Q)
	ringQP.ExtendBasisSmallNormAndCenter(pk.Value[0].Q, levelP, nil, pk.Value[0].P)
	ringQP.NTTLvl(levelQ, levelP, pk.Value[0], pk.Value[0])

	keygen.uniformSamplerQ.Read(pk.Value[1].Q)
	keygen.uniformSamplerP.Read(pk.Value[1].P)

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, sk.Value, pk.Value[1], pk.Value[0])
	return pk
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey(sk.ID)
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *KeyGenerator) GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw, sk.ID)
	return sk, keygen.GenPublicKey(sk)
}

func (keygen *KeyGenerator) genGaussianSampleVector(e *PolyQPVector) {

	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1
	ringQP := keygen.params.RingQP()

	for i := 0; i < e.Dim(); i++ {
		keygen.gaussianSamplerQ.ReadLvl(levelQ, e.Value[i].Q)
		ringQP.ExtendBasisSmallNormAndCenter(e.Value[i].Q, levelP, nil, e.Value[i].P)
		ringQP.NTTLvl(levelQ, levelP, e.Value[i], e.Value[i])
	}

}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in MontgomeryForm
func (keygen *KeyGenerator) GenRelinearizationKey(sk *SecretKey) (rlk *RelinearizationKey) {

	if keygen.params.PCount() == 0 {
		panic("modulus P is empty")
	}
	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1
	ringP := keygen.params.RingP()
	ringQ := keygen.params.RingQ()
	ringQP := RingQP{*keygen.params.RingQP()}

	id := sk.ID

	//rlk = (b, d, v)
	rlk = NewRelinearizationKey(keygen.params, levelQ, levelP, id)
	alpha := levelP + 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	//generate temporary random vector r
	r := keygen.GenSecretKey(sk.ID)

	//temporary vector
	tmp := ringQP.NewPolyVector(beta)

	//generate gadget vector g
	g := ringQP.NewPolyVector(beta)
	for i := 0; i < beta; i++ {
		//compute pq_hat_i
		pq_hat := new(big.Int).Set(ringQ.ModulusBigint)
		pq_hat.Mul(pq_hat, ringP.ModulusBigint)
		for j := 0; j < alpha; j++ {
			index := i*alpha + j

			if index < levelQ+1 {
				pq_hat.Div(pq_hat, ring.NewUint(ringQ.Modulus[index]))
			} else {
				pq_hat.Div(pq_hat, ring.NewUint(ringP.Modulus[index-levelQ-1]))
			}
		}

		ringQ.AddScalarBigint(g.Value[i].Q, pq_hat, g.Value[i].Q)
		ringP.AddScalarBigint(g.Value[i].P, pq_hat, g.Value[i].P)
	}

	//generate vector b
	e := ringQP.NewPolyVector(beta)
	keygen.genGaussianSampleVector(e)
	ringQP.MulPolyMontgomeryLvl(levelQ, levelP, keygen.params.CRS[0], &sk.Value, rlk.Value[0])
	ringQP.AddLvl(levelQ, levelP, e, rlk.Value[0], rlk.Value[0])
	ringQP.MFormLvl(levelQ, levelP, rlk.Value[0], rlk.Value[0])

	//generate vector d
	ringQP.MulPolyMontgomeryLvl(levelQ, levelP, keygen.params.CRS[0], &r.Value, rlk.Value[1])
	ringQP.MulPolyMontgomeryLvl(levelQ, levelP, g, &sk.Value, tmp)
	ringQP.AddLvl(levelQ, levelP, tmp, rlk.Value[1], rlk.Value[1])
	keygen.genGaussianSampleVector(e)
	ringQP.AddLvl(levelQ, levelP, e, rlk.Value[1], rlk.Value[1])
	ringQP.MFormLvl(levelQ, levelP, rlk.Value[1], rlk.Value[1])

	//generate vector v
	ringQP.MulPolyMontgomeryLvl(levelQ, levelP, keygen.params.CRS[1], &sk.Value, rlk.Value[2])
	ringQP.MulPolyMontgomeryLvl(levelQ, levelP, g, &r.Value, tmp)
	ringQP.AddLvl(levelQ, levelP, tmp, rlk.Value[2], rlk.Value[2])
	keygen.genGaussianSampleVector(e)
	ringQP.AddLvl(levelQ, levelP, e, rlk.Value[2], rlk.Value[2])
	ringQP.MFormLvl(levelQ, levelP, rlk.Value[2], rlk.Value[2])

	return
}
