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
func NewParameters(ckksParams ckks.Parameters) *Parameters {

	ret := new(Parameters)
	ret.Parameters = *mkrlwe.NewParameters(ckksParams.Parameters)
	ret.logSlots = ckksParams.LogSlots()
	ret.scale = ckksParams.Scale()

	return ret
}
