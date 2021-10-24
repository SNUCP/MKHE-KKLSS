package mkrlwe

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/utils"

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [1]*ring.Poly
	poolP [3]*ring.Poly

	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  *ring.UniformSampler
}

// Encryptor is a struct used to encrypt plaintext with public key
type Encryptor struct {
	encryptorBase
	pk            *PublicKey
	baseconverter *ring.FastBasisExtender
}

func newEncryptorBase(params Parameters) encryptorBase {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var poolP [3]*ring.Poly
	if params.PCount() != 0 {
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return encryptorBase{
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		poolQ:           [1]*ring.Poly{ringQ.NewPoly()},
		poolP:           poolP,
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSampler:  ring.NewUniformSampler(prng, ringQ),
	}
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {
	id := encryptor.pk.ID

	ringQ := encryptor.ringQ
	ringQP := encryptor.params.RingQP()

	levelQ := utils.MinInt(plaintext.Level(), ctOut.Level())
	levelP := 0

	poolQ0 := encryptor.poolQ[0]
	poolP0 := encryptor.poolP[0]
	poolP1 := encryptor.poolP[1]
	poolP2 := encryptor.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ciphertextNTT := ctOut.Value0.IsNTT

	u := rlwe.PolyQP{Q: poolQ0, P: poolP2}

	encryptor.ternarySampler.ReadLvl(levelQ, u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, nil, u.P)

	// (#Q + #P) NTT
	ringQP.NTTLvl(levelQ, levelP, u, u)
	ringQP.MFormLvl(levelQ, levelP, u, u)

	ct0QP := rlwe.PolyQP{Q: ctOut.Value0, P: poolP0}
	ct1QP := rlwe.PolyQP{Q: ctOut.Value[id], P: poolP1}

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, encryptor.pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, u, encryptor.pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.InvNTTLvl(levelQ, levelP, ct0QP, ct0QP)
	ringQP.InvNTTLvl(levelQ, levelP, ct1QP, ct1QP)

	e := rlwe.PolyQP{Q: poolQ0, P: poolP2}

	encryptor.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct0QP, e, ct0QP)

	encryptor.gaussianSampler.ReadLvl(levelQ, e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.AddLvl(levelQ, levelP, ct1QP, e, ct1QP)

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct0QP.Q)

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct1QP.Q)

	if ciphertextNTT {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ctOut.Value0, plaintext.Value, ctOut.Value0)
		}

		// 2*#Q NTT
		ringQ.NTTLvl(levelQ, ctOut.Value0, ctOut.Value0)
		ringQ.NTTLvl(levelQ, ctOut.Value[id], ctOut.Value[id])

		if plaintext.Value.IsNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(levelQ, ctOut.Value0, plaintext.Value, ctOut.Value0)
		}

	} else {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(levelQ, ctOut.Value0, plaintext.Value, ctOut.Value0)
		} else {
			ringQ.InvNTTLvl(levelQ, plaintext.Value, poolQ0)
			ringQ.AddLvl(levelQ, ctOut.Value0, poolQ0, ctOut.Value0)
		}
	}

	ctOut.Value[id].IsNTT = ctOut.Value0.IsNTT
	ctOut.Value0.Coeffs = ctOut.Value0.Coeffs[:levelQ+1]
	ctOut.Value[id].Coeffs = ctOut.Value[id].Coeffs[:levelQ+1]
}

// NewEncryptor instatiates a new generic RLWE Encryptor. The key argument can
// be either a *rlwe.PublicKey or a *rlwe.SecretKey.
func NewEncryptor(params Parameters, pk *PublicKey) *Encryptor {
	if pk.Value[0].Q.Degree() != params.N() || pk.Value[1].Q.Degree() != params.N() {
		panic("cannot newEncryptor: pk ring degree does not match params ring degree")
	}
	encryptorBase := newEncryptorBase(params)
	baseconverter := ring.NewFastBasisExtender(params.RingQ(), params.RingP())

	return &Encryptor{encryptorBase, pk, baseconverter}
}
