package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"
import "math"

type Parameters struct {
	rlwe.Parameters
	CRS [2]*SwitchingKey
}

// NewParameters takes rlwe Parameter as input, generate two CRSs
// and then return mkrlwe parameter
func NewParameters(params rlwe.Parameters) Parameters {
	ret := new(Parameters)
	ret.Parameters = params
	ringQP := params.RingQP()
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	alpha := params.PCount() / 2
	beta := int(math.Ceil(float64(levelQ+1) / float64(alpha)))

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	uniformSamplerQ := ring.NewUniformSampler(prng, params.RingQ())
	uniformSamplerP := ring.NewUniformSampler(prng, params.RingP())

	ret.CRS[0] = new(SwitchingKey)
	ret.CRS[1] = new(SwitchingKey)

	ret.CRS[0].Value = make([]rlwe.PolyQP, beta)
	ret.CRS[1].Value = make([]rlwe.PolyQP, beta)

	for i := 0; i < beta; i++ {
		ret.CRS[0].Value[i] = ringQP.NewPoly()
		uniformSamplerQ.Read(ret.CRS[0].Value[i].Q)
		uniformSamplerP.Read(ret.CRS[0].Value[i].P)
		ringQP.MFormLvl(levelQ, levelP, ret.CRS[0].Value[i], ret.CRS[0].Value[i])

		ret.CRS[1].Value[i] = ringQP.NewPoly()
		uniformSamplerQ.Read(ret.CRS[1].Value[i].Q)
		uniformSamplerP.Read(ret.CRS[1].Value[i].P)
		ringQP.MFormLvl(levelQ, levelP, ret.CRS[1].Value[i], ret.CRS[1].Value[i])
	}
	return *ret
}

func (params Parameters) Alpha() int {
	return params.PCount() / 2
}

func (params Parameters) Beta(levelQ int) int {
	alpha := params.Alpha()
	beta := int(math.Ceil(float64(levelQ+1) / float64(alpha)))
	return beta
}
