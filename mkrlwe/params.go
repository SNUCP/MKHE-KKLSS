package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"
import "math"

type Parameters struct {
	rlwe.Parameters
	CRS [2]*PolyQPVector
}

// NewParameters takes rlwe Parameter as input, generate two CRSs
// and then return mkrlwe parameter
func NewParameters(params rlwe.Parameters) *Parameters {
	ret := new(Parameters)
	ret.Parameters = params
	ringQP := RingQP{*params.RingQP()}
	levelQ, levelP := params.QCount()-1, params.PCount()-1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	uniformSamplerQ := ring.NewUniformSampler(prng, params.RingQ())
	uniformSamplerP := ring.NewUniformSampler(prng, params.RingP())

	ret.CRS[0] = ringQP.NewPolyVector(beta)
	ret.CRS[1] = ringQP.NewPolyVector(beta)

	for i := 0; i < beta; i++ {
		uniformSamplerQ.Read(ret.CRS[0].Value[i].Q)
		uniformSamplerP.Read(ret.CRS[0].Value[i].P)

		uniformSamplerQ.Read(ret.CRS[1].Value[i].Q)
		uniformSamplerP.Read(ret.CRS[1].Value[i].P)
	}
	return ret
}
