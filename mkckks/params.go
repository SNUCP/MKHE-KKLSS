package mkckks

import "github.com/ldsec/lattigo/v2/ckks"
import "mk-lattigo/mkrlwe"

// Parameters represents a parameter set for the CKKS cryptosystem. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	mkrlwe.Parameters

	logSlots int
	scale    float64
}

// NewParameters instantiate a set of MKCKKS parameters from the generic CKKS parameters and the CKKS-specific ones.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParameters(ckksParams ckks.Parameters) Parameters {

	ret := new(Parameters)
	ret.Parameters = mkrlwe.NewParameters(ckksParams.Parameters)
	ret.logSlots = ckksParams.LogSlots()
	ret.scale = ckksParams.Scale()

	return *ret
}

// Scale returns the default plaintext/ciphertext scale
func (p Parameters) Scale() float64 {
	return p.scale
}

// Slots returns number of available plaintext slots
func (p Parameters) Slots() int {
	return 1 << p.logSlots
}

// LogSlots returns the log of the number of slots
func (p Parameters) LogSlots() int {
	return p.logSlots
}
