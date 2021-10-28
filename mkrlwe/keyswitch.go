package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "math"

// KeySwitcher is a struct for RLWE key-switching.
type MKKeySwitcher struct {
	rlwe.KeySwitcher
	Parameters
}

func NewKeySwitcher(params Parameters) *MKKeySwitcher {
	ks := new(MKKeySwitcher)
	ks.KeySwitcher = *rlwe.NewKeySwitcher(params.Parameters)
	ks.Parameters = params
	return ks
}

// InternalProduct applies internal product of input poly a & bg
// the expected result is ab
// we assume one of a and bg is in MForm
func (ks *MKKeySwitcher) InternalProduct(levelQ int, a *ring.Poly, bg *SwitchingKey, c *ring.Poly) {
	params := ks.Parameters
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()

	var aNTT, aInvNTT *ring.Poly
	if a.IsNTT {
		aNTT = a
		aInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, aNTT, aInvNTT)
	} else {
		panic("Cannot InternalProduct: a should be in NTT")
	}

	alpha := len(bg.Value[0].P.Coeffs)
	levelP := alpha - 1
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))
	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(levelP) >> 1

	c0QP := ks.Pool[0]
	c1QP := ks.Pool[1]

	reduce := 0

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {

		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, aNTT, aInvNTT, c0QP.Q, c0QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, bg.Value[i], c0QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, bg.Value[i], c0QP, c1QP)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}

	// rescale by P
	if a.IsNTT {
		ks.Baseconverter.ModDownQPtoQNTT(levelQ, levelP, c1QP.Q, c1QP.P, c)
	} else {
		ringQ.InvNTTLazyLvl(levelQ, c1QP.Q, c1QP.Q)
		ringP.InvNTTLazyLvl(levelP, c1QP.P, c1QP.P)

		ks.Baseconverter.ModDownQPtoQ(levelQ, levelP, c1QP.Q, c1QP.P, c)
	}
}
