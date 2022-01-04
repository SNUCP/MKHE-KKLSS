package mkbfv

import "mk-lattigo/mkrlwe"

type Ciphertext struct {
	*mkrlwe.Ciphertext
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
// ciphertext is not in NTT Form
func NewCiphertext(params Parameters, idset *mkrlwe.IDSet) (ciphertext *Ciphertext) {
	return &Ciphertext{mkrlwe.NewCiphertext(params.Parameters, idset, params.MaxLevel())}
}

// NewCiphertext creates a new ciphertext parameterized by degree, level and scale.
// ciphertext is in NTT Form
func NewCiphertextNTT(params Parameters, idset *mkrlwe.IDSet) (ciphertext *Ciphertext) {
	return &Ciphertext{mkrlwe.NewCiphertextNTT(params.Parameters, idset, params.MaxLevel())}
}

// CopyNew creates a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() *Ciphertext {
	return &Ciphertext{ct.Ciphertext.CopyNew()}
}

type Message struct {
	Value []int64
}

func NewMessage(params Parameters) *Message {

	msg := new(Message)
	msg.Value = make([]int64, params.N())

	return msg
}

func (msg *Message) Slots() int {
	return len(msg.Value)
}
