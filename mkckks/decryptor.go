package mkckks

import "github.com/ldsec/lattigo/v2/ckks"
import "mk-lattigo/mkrlwe"

type Decryptor struct {
	*mkrlwe.Decryptor
	encoder  ckks.Encoder
	params   Parameters
	ptxtPool *ckks.Plaintext
}

// NewDecryptor instantiates a Decryptor for the CKKS scheme.
func NewDecryptor(params Parameters) *Decryptor {
	ckksParams, _ := ckks.NewParameters(params.Parameters.Parameters, params.LogSlots(), params.Scale())

	ret := new(Decryptor)
	ret.Decryptor = mkrlwe.NewDecryptor(params.Parameters)
	ret.encoder = ckks.NewEncoder(ckksParams)
	ret.params = params
	ret.ptxtPool = ckks.NewPlaintext(ckksParams, params.MaxLevel(), params.Scale())
	ret.ptxtPool.Value.IsNTT = false
	return ret
}

// PartialDecrypt partially decrypts the ct with single secretkey sk and update result inplace
func (dec *Decryptor) PartialDecrypt(ct *Ciphertext, sk *mkrlwe.SecretKey) {
	dec.Decryptor.PartialDecrypt(ct.Ciphertext, sk)
}

// Decrypt decrypts the ciphertext with given secretkey set and write the result in ptOut.
// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
// Output domain will match plaintext.Value.IsNTT value.
func (dec *Decryptor) Decrypt(ciphertext *Ciphertext, skSet *mkrlwe.SecretKeySet) (msg *Message) {
	ctTmp := ciphertext.CopyNew()

	dec.Decryptor.Decrypt(ctTmp.Ciphertext, skSet, dec.ptxtPool.Plaintext)
	dec.ptxtPool.Scale = ctTmp.Scale
	msg = new(Message)
	msg.Value = dec.encoder.Decode(dec.ptxtPool, dec.params.logSlots)

	return
}
