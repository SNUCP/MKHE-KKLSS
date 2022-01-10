package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"

// KeySwitcher is a struct for RLWE key-switching.
type KeySwitcher struct {
	rlwe.KeySwitcher
	Parameters
	Decomposer *Decomposer
	polyQPool  [3]*ring.Poly
	swkPool1   *SwitchingKey
	swkPool2   *SwitchingKey
	swkPool3   *SwitchingKey
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (ks *KeySwitcher) DecomposeSingleNTT(levelQ, levelP, alpha, beta, gamma int, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := ks.Parameters.RingQ()
	ringP := ks.Parameters.RingP()

	ks.Decomposer.DecomposeAndSplit(levelQ, levelP, alpha, beta, gamma, c2InvNTT, c2QiQ, c2QiP)

	// c2QiP = c2 mod qi mod pj
	ringQ.NTTLvl(levelQ, c2QiQ, c2QiQ)
	ringP.NTTLvl(levelP, c2QiP, c2QiP)
}

func NewKeySwitcher(params Parameters) *KeySwitcher {
	ks := new(KeySwitcher)
	ks.KeySwitcher = *rlwe.NewKeySwitcher(params.Parameters)
	ks.Parameters = params
	ks.Decomposer = NewDecomposer(params.RingQ(), params.RingP(), params.Gamma())

	ringQ := params.RingQ()
	ks.polyQPool = [3]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}

	ks.swkPool1 = NewSwitchingKey(params)
	ks.swkPool2 = NewSwitchingKey(params)
	ks.swkPool3 = NewSwitchingKey(params)

	return ks
}

func (ks *KeySwitcher) Decompose(levelQ int, a *ring.Poly, ad *SwitchingKey) {

	params := ks.Parameters
	ringQ := params.RingQ()

	var aNTT, aInvNTT *ring.Poly
	if a.IsNTT {
		aNTT = a
		aInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, aNTT, aInvNTT)
	} else {
		aInvNTT = a
	}

	alpha := params.Alpha()
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {
		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, params.Gamma(), aInvNTT, ad.Value[i].Q, ad.Value[i].P)
	}

	return
}

// InternalProduct applies internal product of input poly a & bg
// the expected result is ab
// we assume one of a and bg is in MForm
// we assume output is in InvNTT form
func (ks *KeySwitcher) InternalProduct(levelQ int, a *ring.Poly, bg *SwitchingKey, c *ring.Poly) {
	params := ks.Parameters
	ringQ := params.RingQ()
	ringP := params.RingP()
	ringQP := params.RingQP()

	var aInvNTT *ring.Poly
	if a.IsNTT {
		aInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, a, aInvNTT)
	} else {
		aInvNTT = a
		//panic("Cannot InternalProduct: a should be in NTT")
	}

	alpha := params.Alpha()
	levelP := params.PCount() - 1
	beta := params.Beta(levelQ)

	c0QP := ks.Pool[0]
	c1QP := ks.Pool[1]

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {
		ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, params.Gamma(), aInvNTT, c0QP.Q, c0QP.P)

		if i == 0 {
			ringQP.MulCoeffsMontgomeryLvl(levelQ, levelP, bg.Value[i], c0QP, c1QP)
		} else {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, bg.Value[i], c0QP, c1QP)
		}
	}

	// rescale by P

	ringQ.InvNTTLazyLvl(levelQ, c1QP.Q, c1QP.Q)
	ringP.InvNTTLazyLvl(levelP, c1QP.P, c1QP.P)

	ks.Baseconverter.ModDownQPtoQ(levelQ, levelP, c1QP.Q, c1QP.P, c)
}

// MulRelin multiplies op0 with op1 with relinearization and returns the result in ctOut.
// Input ciphertext should be in NTT form
func (ks *KeySwitcher) MulAndRelin(op0, op1 *Ciphertext, rlkSet *RelinearizationKeySet, ctOut *Ciphertext) {

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
		ks.Decompose(level, op0.Value[id], ks.swkPool3)
		d := rlkSet.Value[id].Value[1]
		for i := 0; i < beta; i++ {
			ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, d.Value[i], ks.swkPool3.Value[i], x.Value[i])
		}
	}

	for i := 0; i < beta; i++ {
		ringQP.MFormLvl(level, levelP, x.Value[i], x.Value[i])
	}

	//gen y vector
	for id := range idset1.Value {
		ks.Decompose(level, op1.Value[id], ks.swkPool3)
		b := rlkSet.Value[id].Value[0]
		for i := 0; i < beta; i++ {
			ringQP.MulCoeffsMontgomeryAndAddLvl(level, levelP, b.Value[i], ks.swkPool3.Value[i], y.Value[i])
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
		ks.InternalProduct(level, op1.Value[id], x, ks.polyQPool[0])
		ringQ.AddLvl(level, ctOut.Value[id], ks.polyQPool[0], ctOut.Value[id])
	}

	//ctOut_0 <- ctOut_0 + Inter(Inter(op0_i, y), v_i)
	//ctOut_i <- ctOut_i + Inter(Inter(op0_i, y), u)

	u := params.CRS[-1]

	for id := range idset0.Value {
		v := rlkSet.Value[id].Value[2]
		ks.InternalProduct(level, op0.Value[id], y, ks.polyQPool[0])

		ks.InternalProduct(level, ks.polyQPool[0], v, ks.polyQPool[1])
		ringQ.AddLvl(level, ctOut.Value["0"], ks.polyQPool[1], ctOut.Value["0"])

		ks.InternalProduct(level, ks.polyQPool[0], u, ks.polyQPool[2])
		ringQ.AddLvl(level, ctOut.Value[id], ks.polyQPool[2], ctOut.Value[id])
	}

}

// Rotate rotates ctIn with ctOut with RotationKeySet and returns the result in ctOut.
// Input ciphertext should be in InvNTT form
func (ks *KeySwitcher) Rotate(ctIn *Ciphertext, rotidx int, rkSet *RotationKeySet, ctOut *Ciphertext) {
	level := ctOut.Level()
	idset := ctIn.IDSet()
	params := ks.Parameters
	ringQ := params.RingQ()

	// check ctIn level
	if ctIn.Level() < level {
		panic("Cannot Rotate: ctIn and ctOut have different levels")
	}

	// adjust rotidx
	for rotidx < 0 {
		rotidx += (params.N() / 2)
	}

	galEl := params.GaloisElementForColumnRotationBy(rotidx)

	// permute ctIn and put it to ctOut
	for id := range ctIn.Value {
		ringQ.Permute(ctIn.Value[id], galEl, ctOut.Value[id])
	}

	// c0 <- c0 + IP(c_i, rk_i)
	for id := range idset.Value {
		rk := rkSet.GetRotationKey(id, uint(rotidx))
		ks.InternalProduct(level, ctOut.Value[id], rk.Value, ks.polyQPool[0])
		ringQ.AddLvl(level, ctOut.Value["0"], ks.polyQPool[0], ctOut.Value["0"])
	}

	// c_i <- IP(c_i, a)
	a := params.CRS[rotidx]
	for id := range idset.Value {
		ks.InternalProduct(level, ctOut.Value[id], a, ks.polyQPool[0])
		ctOut.Value[id].Copy(ks.polyQPool[0])
	}
}

// Conjugate conjugate ctIn with ctOut with ConjugationKeySet and returns the result in ctOut.
// Input ciphertext should be in NTT form
func (ks *KeySwitcher) Conjugate(ctIn *Ciphertext, ckSet *ConjugationKeySet, ctOut *Ciphertext) {
	level := ctOut.Level()
	idset := ctIn.IDSet()
	params := ks.Parameters
	ringQ := params.RingQ()
	galEl := params.GaloisElementForRowRotation()

	// check ctIn level
	if ctIn.Level() < level {
		panic("Cannot Conjugate: ctIn and ctOut have different levels")
	}

	// permute ctIn and put it to ctOut
	for id := range ctIn.Value {
		ringQ.Permute(ctIn.Value[id], galEl, ctOut.Value[id])
	}

	// c0 <- c0 + IP(c_i, rk_i)
	for id := range idset.Value {
		ck := ckSet.GetConjugationKey(id)
		ks.InternalProduct(level, ctOut.Value[id], ck.Value, ks.polyQPool[0])
		ringQ.AddLvl(level, ctOut.Value["0"], ks.polyQPool[0], ctOut.Value["0"])
	}

	// c_i <- IP(c_i, a)
	a := params.CRS[-2]
	for id := range idset.Value {
		ks.InternalProduct(level, ctOut.Value[id], a, ks.polyQPool[0])
		ctOut.Value[id].Copy(ks.polyQPool[0])
	}
}
