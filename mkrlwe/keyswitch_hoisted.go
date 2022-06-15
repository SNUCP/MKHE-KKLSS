package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"

// ExternalProduct applies internal product of input poly a & bg
// the expected result is ab
// we assume one of a and bg is in MForm
// we assume output is in InvNTT form
func (ks *KeySwitcher) ExternalProductHoisted(levelQ int, aHoisted []rlwe.PolyQP, bg *SwitchingKey, c *ring.Poly) {
	params := ks.Parameters
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()

	//alpha := params.Alpha()
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)

	//c0QP := ks.Pool[0]
	c1QP := ks.Pool[1]

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {
		//ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, params.Gamma(), aInvNTT, c0QP.Q, c0QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, bg.Value[i], aHoisted[i], c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, bg.Value[i], aHoisted[i], c1QP)
		}
	}

	// rescale by P

	ringQ.InvNTTLazyLvl(levelQ, c1QP.Q, c1QP.Q)
	ringP.InvNTTLazyLvl(levelP, c1QP.P, c1QP.P)

	ks.Baseconverter.ModDownQPtoQ(levelQ, levelP, c1QP.Q, c1QP.P, c)
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// Input ciphertext should be in NTT form
func (ks *KeySwitcher) MulAndRelinHoisted(op0, op1 *Ciphertext, op0Hoisted, op1Hoisted map[string]*SwitchingKey, rlkSet *RelinearizationKeySet, ctOut *Ciphertext) {

	level := ctOut.Level()

	if op0.Level() < level {
		panic("Cannot MulAndRelin: op0 and op1 have different levels")
	}

	if ctOut.Level() < level {
		panic("Cannot MulAndRelin: op0 and ctOut have different levels")
	}

	idset0 := op0.IDSet()
	idset1 := op1.IDSet()

	params := ks.Parameters
	ringQP := params.RingQP()
	ringQ := params.RingQ()

	levelP := params.PCount() - 1
	beta := params.Beta(level)

	x := ks.swkPool1
	y := ks.swkPool2

	//initialize x1, x2, y1, y2
	for i := 0; i < beta; i++ {
		x.Value[i].Q.Zero()
		x.Value[i].P.Zero()

		y.Value[i].Q.Zero()
		y.Value[i].P.Zero()
	}

	//gen x vector
	for id := range idset0.Value {
		if op0Hoisted == nil {
			ks.Decompose(level, op0.Value[id], ks.swkPool3)
			d := rlkSet.Value[id].Value[1]
			for i := 0; i < beta; i++ {
				ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, d.Value[i], ks.swkPool3.Value[i], x.Value[i])
			}
		} else {
			d := rlkSet.Value[id].Value[1]
			for i := 0; i < beta; i++ {
				ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, d.Value[i], op0Hoisted[id].Value[i], x.Value[i])
			}
		}
	}

	for i := 0; i < beta; i++ {
		ringQP.MFormLvl(level, levelP, x.Value[i], x.Value[i])
	}

	//gen y vector
	for id := range idset1.Value {
		if op1Hoisted == nil {
			ks.Decompose(level, op1.Value[id], ks.swkPool3)
			b := rlkSet.Value[id].Value[0]
			for i := 0; i < beta; i++ {
				ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, b.Value[i], ks.swkPool3.Value[i], y.Value[i])
			}
		} else {
			b := rlkSet.Value[id].Value[0]
			for i := 0; i < beta; i++ {
				ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, b.Value[i], op1Hoisted[id].Value[i], y.Value[i])
			}
		}

	}

	for i := 0; i < beta; i++ {
		ringQP.MFormLvl(level, levelP, y.Value[i], y.Value[i])
	}

	//ctOut_0 <- op0_0 * op1_0
	ringQ.NTTLvl(level, op0.Value["0"], ks.polyQPool[0])
	ringQ.NTTLvl(level, op1.Value["0"], ks.polyQPool[1])

	ringQ.MFormLvl(level, ks.polyQPool[0], ks.polyQPool[0])
	ringQ.MulCoeffsMontgomeryLvl(level, ks.polyQPool[0], ks.polyQPool[1], ctOut.Value["0"])

	//ctOut_j <- op0_0 * op1_j + op0_j * op1_0
	ringQ.MFormLvl(level, ks.polyQPool[1], ks.polyQPool[1])
	for id := range idset0.Value {
		ringQ.NTTLvl(level, op0.Value[id], ks.polyQPool[2])
		ringQ.MulCoeffsMontgomeryLvl(level, ks.polyQPool[1], ks.polyQPool[2], ctOut.Value[id])
	}

	for id := range idset1.Value {
		ringQ.NTTLvl(level, op1.Value[id], ks.polyQPool[2])
		if idset0.Has(id) {
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, ks.polyQPool[0], ks.polyQPool[2], ctOut.Value[id])
		} else {
			ringQ.MulCoeffsMontgomeryLvl(level, ks.polyQPool[0], ks.polyQPool[2], ctOut.Value[id])
		}
	}

	for id := range ctOut.Value {
		ringQ.InvNTTLvl(level, ctOut.Value[id], ctOut.Value[id])
	}

	//ctOut_j <- ctOut_j +  Inter(op1_j, x)
	for id := range idset1.Value {
		if op1Hoisted == nil {
			ks.ExternalProduct(level, op1.Value[id], x, ks.polyQPool[0])
		} else {
			ks.ExternalProductHoisted(level, op1Hoisted[id].Value, x, ks.polyQPool[0])
		}
		ringQ.AddLvl(level, ctOut.Value[id], ks.polyQPool[0], ctOut.Value[id])
	}

	//ctOut_0 <- ctOut_0 + Inter(Inter(op0_i, y), v_i)
	//ctOut_i <- ctOut_i + Inter(Inter(op0_i, y), u)

	u := params.CRS[-1]

	for id := range idset0.Value {

		v := rlkSet.Value[id].Value[2]

		if op0Hoisted == nil {
			ks.ExternalProduct(level, op0.Value[id], y, ks.polyQPool[0])
		} else {
			ks.ExternalProductHoisted(level, op0Hoisted[id].Value, y, ks.polyQPool[0])
		}

		ks.Decompose(level, ks.polyQPool[0], ks.swkPool3)

		ks.ExternalProductHoisted(level, ks.swkPool3.Value, v, ks.polyQPool[1])
		ringQ.AddLvl(level, ctOut.Value["0"], ks.polyQPool[1], ctOut.Value["0"])

		ks.ExternalProductHoisted(level, ks.swkPool3.Value, u, ks.polyQPool[2])
		ringQ.AddLvl(level, ctOut.Value[id], ks.polyQPool[2], ctOut.Value[id])
	}
}
