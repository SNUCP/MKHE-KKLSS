package mkbfv

import "mk-lattigo/mkrlwe"
import "github.com/ldsec/lattigo/v2/bfv"

type Decryptor struct {
	*mkrlwe.Decryptor
	params   Parameters
	encoder  bfv.Encoder
	ptxtPool *bfv.Plaintext
}

// NewDecryptor instantiates a Decryptor for the CKKS scheme.
func NewDecryptor(params Parameters) *Decryptor {
	bfvParams, _ := bfv.NewParameters(params.paramsQP.Parameters, params.T())

	ret := new(Decryptor)
	ret.Decryptor = mkrlwe.NewDecryptor(params.paramsQP)
	ret.encoder = bfv.NewEncoder(bfvParams)
	ret.params = params
	ret.ptxtPool = bfv.NewPlaintext(bfvParams)

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
	dec.Decryptor.Decrypt(ciphertext.Ciphertext, skSet, dec.ptxtPool.Plaintext)
	msg = new(Message)
	dec.encoder.DecodeInt(dec.ptxtPool, msg.Value)

	return
}
