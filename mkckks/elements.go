package mkckks

import "mk-lattigo/mkrlwe"

type Ciphertext struct {
	*mkrlwe.Ciphertext
	Scale float64
}

// NewCiphertext returns a new Element with zero values
func NewCiphertext(params Parameters, idset []string, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Ciphertext = mkrlwe.NewCiphertextNTT(params.Parameters, idset, level)
	el.Scale = params.scale

	return el
}

// ScalingFactor returns the scaling factor of the ciphertext
func (ct *Ciphertext) ScalingFactor() float64 {
	return ct.Scale
}

// SetScalingFactor sets the scaling factor of the ciphertext
func (ct *Ciphertext) SetScalingFactor(scale float64) {
	ct.Scale = scale
}

// Copy copies the given ciphertext ctp into the receiver ciphertext.
func (ct *Ciphertext) Copy(ctp *Ciphertext) {
	ct.Ciphertext.Copy(ctp.Ciphertext)
	ct.Scale = ctp.Scale
}

// CopyNew makes a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() (ctc *Ciphertext) {
	ctc = &Ciphertext{Ciphertext: ct.Ciphertext.CopyNew(), Scale: ct.Scale}
	return
}

type Message struct {
	Value []complex128
}

func NewMessage(params Parameters) *Message {

	msg := new(Message)
	msg.Value = make([]complex128, (1 << params.logSlots))

	return msg
}

func (msg *Message) Slots() int {
	return len(msg.Value)
}
