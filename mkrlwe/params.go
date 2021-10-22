package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"

type Parameters struct {
	rlwe.Parameters
	CRS [2]*PolyQPVector
}

func NewParameters(rlweParams rlwe.Parameters, crs0, crs1 *PolyQPVector) *Parameters {
	ret := new(Parameters)
	ret.Parameters = rlweParams
	ret.CRS[0] = crs0
	ret.CRS[1] = crs1

	return ret
}
