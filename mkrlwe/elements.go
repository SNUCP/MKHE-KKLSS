package mkrlwe

import "github.com/ldsec/lattigo/v2/ring"

type Ciphertext struct {
	Value0 *ring.Poly
	Value  map[string]*ring.Poly
}

// NewCiphertext returns a new Element with zero values
func NewCiphertext(params Parameters, idset []string, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value0 = ring.NewPoly(params.N(), level+1)
	el.Value = make(map[string]*ring.Poly)

	for _, id := range idset {
		el.Value[id] = ring.NewPoly(params.N(), level+1)
	}

	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set
func NewCiphertextNTT(params Parameters, idset []string, level int) *Ciphertext {
	el := NewCiphertext(params, idset, level)
	el.Value0 = ring.NewPoly(params.N(), level+1)
	el.Value0.IsNTT = true

	el.Value = make(map[string]*ring.Poly)

	for _, id := range idset {
		el.Value[id] = ring.NewPoly(params.N(), level+1)
		el.Value[id].IsNTT = true
	}

	return el
}

// NumID returns the number of IDs engaged in the target elements
func (el *Ciphertext) NumIDs() int {
	return len(el.Value)
}

// Level returns the level of the target elements
func (el *Ciphertext) Level() int {
	return len(el.Value0.Coeffs) - 1
}

// CopyNew creates a new element as a copy of the target element.
func (el *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)
	ctxCopy.Value = make(map[string]*ring.Poly)

	ctxCopy.Value0 = el.Value0.CopyNew()
	for id := range el.Value {
		ctxCopy.Value[id] = el.Value[id].CopyNew()
	}

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Ciphertext) Copy(ctxCopy *Ciphertext) {

	// clear target element
	for id := range ctxCopy.Value {
		delete(el.Value, id)
	}

	if el != ctxCopy {
		for id := range ctxCopy.Value {
			el.Value[id].Copy(ctxCopy.Value[id])
		}
		el.Value0.Copy(ctxCopy.Value0)
	}

}
