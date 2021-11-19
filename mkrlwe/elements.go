package mkrlwe

import "github.com/ldsec/lattigo/v2/ring"

type Ciphertext struct {
	Value map[string]*ring.Poly
}

// NewCiphertext returns a new Element with zero values
func NewCiphertext(params Parameters, idset *IDSet, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make(map[string]*ring.Poly)

	el.Value["0"] = ring.NewPoly(params.N(), level+1)

	for id := range idset.Value {
		el.Value[id] = ring.NewPoly(params.N(), level+1)
	}

	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set
func NewCiphertextNTT(params Parameters, idset *IDSet, level int) *Ciphertext {
	el := NewCiphertext(params, idset, level)

	for id := range el.Value {
		el.Value[id].IsNTT = true
	}

	return el
}

// NumID returns the number of IDs engaged in the target elements
func (el *Ciphertext) IDSet() *IDSet {
	idset := NewIDSet()

	for id := range el.Value {
		if id != "0" {
			idset.Add(id)
		}
	}

	return idset
}

// Level returns the level of the target elements
func (el *Ciphertext) Level() int {
	return len(el.Value["0"].Coeffs) - 1
}

// El returns a pointer to this Element
func (el *Ciphertext) El() *Ciphertext {
	return el
}

// CopyNew creates a new element as a copy of the target element.
func (el *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)
	ctxCopy.Value = make(map[string]*ring.Poly)

	ctxCopy.Value["0"] = el.Value["0"].CopyNew()
	for id := range el.Value {
		ctxCopy.Value[id] = el.Value[id].CopyNew()
	}

	return ctxCopy
}

// Copy copy value from the input element.
func (el *Ciphertext) Copy(ct *Ciphertext) {
	for id := range ct.Value {
		el.Value[id].Copy(ct.Value[id])
	}
}

// PadCiphertext pads a ciphertext with an input idset
func (el *Ciphertext) PadCiphertext(idset *IDSet) {

	isNTT := el.Value["0"].IsNTT
	oldIDs := el.IDSet()

	degree := el.Value["0"].Degree()
	level := el.Level()

	for id := range idset.Value {
		if !oldIDs.Has(id) {
			el.Value[id] = ring.NewPoly(degree, level+1)
			el.Value[id].IsNTT = isNTT
		}
	}
}
