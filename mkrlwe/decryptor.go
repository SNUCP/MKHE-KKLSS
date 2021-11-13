package mkrlwe

import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/utils"

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type Decryptor struct {
	params Parameters
	ringQ  *ring.Ring
	pool   *ring.Poly
	sk     *SecretKey
}

// NewDecryptor instantiates a new generic RLWE Decryptor.
func NewDecryptor(params Parameters) *Decryptor {

	return &Decryptor{
		params: params,
		ringQ:  params.RingQ(),
		pool:   params.RingQ().NewPoly(),
	}
}

// PartialDecrypt partially decrypts the ct with single secretkey sk and update result inplace
func (decryptor *Decryptor) PartialDecrypt(ct *Ciphertext, sk *SecretKey) {
	ringQ := decryptor.ringQ
	id := sk.ID
	level := ct.Level()

	if !ct.Value[id].IsNTT {
		ringQ.NTTLvl(level, ct.Value[id], ct.Value[id])
	}

	ringQ.MulCoeffsMontgomeryLvl(level, ct.Value[id], sk.Value.Q, ct.Value[id])

	if !ct.Value[id].IsNTT {
		ringQ.InvNTTLvl(level, ct.Value[id], ct.Value[id])
	}

	ringQ.AddLvl(level, ct.Value["0"], ct.Value[id], ct.Value["0"])
	delete(ct.Value, id)
}

// Decrypt decrypts the ciphertext with given secretkey set and write the result in ptOut.
// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
// Output domain will match plaintext.Value.IsNTT value.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, skSet *SecretKeySet, plaintext *rlwe.Plaintext) {
	ringQ := decryptor.ringQ
	level := utils.MinInt(ciphertext.Level(), plaintext.Level())
	plaintext.Value.Coeffs = plaintext.Value.Coeffs[:level+1]

	ctTmp := ciphertext.CopyNew()
	idset := ctTmp.IDSet()
	for _, sk := range skSet.Value {
		if idset.Has(sk.ID) {
			decryptor.PartialDecrypt(ctTmp, sk)
		}
	}

	if len(ctTmp.Value) > 1 {
		panic("Cannot Decrypt: there is a missing secretkey")
	}

	ringQ.ReduceLvl(level, ctTmp.Value["0"], plaintext.Value)
}
