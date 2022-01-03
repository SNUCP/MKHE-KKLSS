package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"
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

	//set a to CRS[0][0]
	pk.Value[1].Q.Copy(keygen.params.CRS[0].Value[0].Q)
	pk.Value[1].P.Copy(keygen.params.CRS[0].Value[0].P)

	ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, sk.Value, pk.Value[1], pk.Value[0])
	return pk
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *KeyGenerator) GenKeyPair(id string) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey(id)
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *KeyGenerator) GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw, sk.ID)
	return sk, keygen.GenPublicKey(sk)
}

func (keygen *KeyGenerator) GenGaussianError(e rlwe.PolyQP) {

	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1
	ringQP := keygen.params.RingQP()

	keygen.gaussianSamplerQ.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.NTTLvl(levelQ, levelP, e, e)

}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize constant term of Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in  MontgomeryForm
func (keygen *KeyGenerator) GenConstRelinearizationKey() (rlk *RelinearizationKey) {
	params := keygen.params
	ringQ := params.RingQ()
	ringP := params.RingP()

	sk := NewSecretKey(params, "0")
	r := NewSecretKey(params, "0")

	ringQ.AddScalar(sk.Value.Q, 1, sk.Value.Q)
	ringP.AddScalar(sk.Value.P, 1, sk.Value.P)

	ringQ.MForm(sk.Value.Q, sk.Value.Q)
	ringP.MForm(sk.Value.P, sk.Value.P)

	ringQ.MForm(r.Value.Q, r.Value.Q)
	ringP.MForm(r.Value.P, r.Value.P)

	rlk = keygen.GenRelinearizationKey(sk, r)
	return
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
// RelinearizationKeys are triplet of polyvector in  MontgomeryForm
func (keygen *KeyGenerator) GenRelinearizationKey(sk, r *SecretKey) (rlk *RelinearizationKey) {

	if keygen.params.PCount() == 0 {
		panic("modulus P is empty")
	}
	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQP := params.RingQP()
	ringQ := params.RingQ()
	ringP := params.RingP()

	id := sk.ID

	//rlk = (b, d, v)
	rlk = NewRelinearizationKey(keygen.params, id)
	beta := params.Beta(levelQ)

	//set CRS
	a := keygen.params.CRS[0]
	u := keygen.params.CRS[-1]

	tmp := keygen.poolQP

	//generate vector b = -sa + e in MForm
	b := rlk.Value[0]
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, a.Value[i], sk.Value, b.Value[i])
		ringQP.InvMFormLvl(levelQ, levelP, b.Value[i], b.Value[i])
		keygen.GenGaussianError(tmp)
		ringQP.SubLvl(levelQ, levelP, tmp, b.Value[i], b.Value[i])
		ringQP.MFormLvl(levelQ, levelP, b.Value[i], b.Value[i])
	}

	//generate vector d = -ra + sg + e in MForm
	d := rlk.Value[1]
	keygen.GenSwitchingKey(sk, d)
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a.Value[i], r.Value, d.Value[i])
	}

	//generate vector v = -su - rg + e in MForm
	v := rlk.Value[2]
	keygen.GenSwitchingKey(r, v)
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u.Value[i], sk.Value, v.Value[i])
		ringQ.NegLvl(levelQ, v.Value[i].Q, v.Value[i].Q)
		ringP.NegLvl(levelP, v.Value[i].P, v.Value[i].P)
	}
	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
func (keygen *KeyGenerator) GenRotationKey(rotidx int, sk *SecretKey) (rk *RotationKey) {
	skIn := sk
	id := sk.ID
	skOut := NewSecretKey(keygen.params, id)

	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)
	ringQP := params.RingQP()

	// check CRS for given rot idx exists
	_, in := params.CRS[rotidx]
	if !in {
		panic("Cannot GenRotationKey: CRS for given rot idx is not generated")
	}

	// adjust rotidx
	for rotidx < 0 {
		rotidx += (params.N() / 2)
	}

	galEl := keygen.params.GaloisElementForColumnRotationBy(rotidx)
	index := ring.PermuteNTTIndex(galEl, uint64(params.N()))
	ring.PermuteNTTWithIndexLvl(params.QCount()-1, skIn.Value.Q, index, skOut.Value.Q)
	ring.PermuteNTTWithIndexLvl(params.PCount()-1, skIn.Value.P, index, skOut.Value.P)

	// rk  = Ps' + e
	rk = NewRotationKey(params, uint(rotidx), id)
	keygen.GenSwitchingKey(skOut, rk.Value)
	a := params.CRS[rotidx]

	// rk = -sa + Ps' + e
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a.Value[i], sk.Value, rk.Value.Value[i])
	}

	return rk
}

// GenRotationKeys generates a RotationKeys of rotidx power of 2 and add it to rtkSet
func (keygen *KeyGenerator) GenDefaultRotationKeys(sk *SecretKey, rtkSet *RotationKeySet) {
	for rotidx := 1; rotidx < keygen.params.N()/2; rotidx *= 2 {
		rtk := keygen.GenRotationKey(rotidx, sk)
		rtkSet.AddRotationKey(rtk)
	}
}

// GenConjugationKeys generates a ConjugationKeySet from a list of galois element corresponding to the desired conjugation
func (keygen *KeyGenerator) GenConjugationKey(sk *SecretKey) (cjk *ConjugationKey) {
	skIn := sk
	id := sk.ID
	skOut := NewSecretKey(keygen.params, id)

	params := keygen.params
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)
	ringQP := params.RingQP()

	galEl := keygen.params.GaloisElementForRowRotation()
	index := ring.PermuteNTTIndex(galEl, uint64(params.N()))
	ring.PermuteNTTWithIndexLvl(params.QCount()-1, skIn.Value.Q, index, skOut.Value.Q)
	ring.PermuteNTTWithIndexLvl(params.PCount()-1, skIn.Value.P, index, skOut.Value.P)

	// rk  = Ps' + e
	cjk = NewConjugationKey(params, id)
	keygen.GenSwitchingKey(skOut, cjk.Value)
	a := params.CRS[-2]

	// rk = -sa + Ps' + e
	for i := 0; i < beta; i++ {
		ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a.Value[i], sk.Value, cjk.Value.Value[i])
	}

	return cjk
}

//For an input secretkey s, gen gs + e in MForm
func (keygen *KeyGenerator) GenSwitchingKey(skIn *SecretKey, swk *SwitchingKey) {
	params := keygen.params
	ringQ := params.RingQ()
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	alpha := params.Alpha()
	beta := params.Beta(levelQ)

	var pBigInt *big.Int
	if levelP == keygen.params.PCount()-1 {
		pBigInt = keygen.params.RingP().ModulusBigint
	} else {
		P := keygen.params.RingP().Modulus
		pBigInt = new(big.Int).SetUint64(P[0])
		for i := 1; i < levelP+1; i++ {
			pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
		}
	}

	// Computes P * skIn
	ringQ.MulScalarBigintLvl(levelQ, skIn.Value.Q, pBigInt, keygen.poolQ)

	var index int
	for i := 0; i < beta; i++ {

		// e

		ringQP := params.RingQP()
		keygen.gaussianSamplerQ.ReadLvl(levelQ, swk.Value[i].Q)
		ringQP.ExtendBasisSmallNormAndCenter(swk.Value[i].Q, levelP, nil, swk.Value[i].P)
		ringQP.NTTLvl(levelQ, levelP, swk.Value[i], swk.Value[i])
		ringQP.MFormLvl(levelQ, levelP, swk.Value[i], swk.Value[i])
		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			// It handles the case where nb pj does not divide nb qi
			if index >= levelQ+1 {
				break
			}

			qi := ringQ.Modulus[index]
			p0tmp := keygen.poolQ.Coeffs[index]
			p1tmp := swk.Value[i].Q.Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}
		}
	}

}
