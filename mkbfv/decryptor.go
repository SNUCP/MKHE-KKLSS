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
func (dec *Decryptor) PartialDecrypt(ct *Ciphertext, sk *SecretKey) {
	dec.Decryptor.PartialDecrypt(ct.Ciphertext, sk.SecretKey)
}

// Decrypt decrypts the ciphertext with given secretkey set and write the result in ptOut.
// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
// Output domain will match plaintext.Value.IsNTT value.
func (dec *Decryptor) Decrypt(ciphertext *Ciphertext, skSet *SecretKeySet) (msg *Message) {

	ctTmp := ciphertext.CopyNew()
	idset := ctTmp.IDSet()
	for _, sk := range skSet.Value {
		if idset.Has(sk.ID) {
			dec.PartialDecrypt(ctTmp, sk)
		}
	}

	msg = NewMessage(dec.params)
	dec.ptxtPool.Value = ctTmp.Value["0"]
	dec.encoder.DecodeInt(dec.ptxtPool, msg.Value)

	return
}
