package mkbfv

import "mk-lattigo/mkrlwe"
import "github.com/ldsec/lattigo/v2/ring"

type RelinearizationKey struct {
	Value [2]*mkrlwe.RelinearizationKey
	ID    string
}

// RelinearizationKeySet is a type for a set of multikey BFV relinearization keys.
type RelinearizationKeySet struct {
	Value map[string]*RelinearizationKey

	PolyRPool1 map[string]*ring.Poly
	PolyRPool2 map[string]*ring.Poly

	HoistPool1 [2]*mkrlwe.HoistedCiphertext
	HoistPool2 [2]*mkrlwe.HoistedCiphertext

	params Parameters
}

// NewRelinearizationKey returns a new RelinearizationKey with zero values.
func NewRelinearizationKey(params Parameters, id string) *RelinearizationKey {
	rlk := new(RelinearizationKey)
	rlk.Value[0] = mkrlwe.NewRelinearizationKey(params.Parameters, id)
	rlk.Value[1] = mkrlwe.NewRelinearizationKey(params.Parameters, id)

	rlk.ID = id

	return rlk
}

// NewRelinearizationKeySet returns a new empty RelinearizationKeySet
func NewRelinearizationKeyKeySet(params Parameters) *RelinearizationKeySet {
	rlkSet := new(RelinearizationKeySet)
	rlkSet.Value = make(map[string]*RelinearizationKey)
	rlkSet.params = params

	rlkSet.PolyRPool1 = make(map[string]*ring.Poly)
	rlkSet.PolyRPool2 = make(map[string]*ring.Poly)

	rlkSet.PolyRPool1["0"] = params.RingR().NewPoly()
	rlkSet.PolyRPool2["0"] = params.RingR().NewPoly()

	rlkSet.HoistPool1[0] = mkrlwe.NewHoistedCiphertext()
	rlkSet.HoistPool1[1] = mkrlwe.NewHoistedCiphertext()

	rlkSet.HoistPool2[0] = mkrlwe.NewHoistedCiphertext()
	rlkSet.HoistPool2[1] = mkrlwe.NewHoistedCiphertext()

	return rlkSet
}

// AddRelinearizationKey insert new publickey into RelinearizationKeySet with its id
func (rlkSet *RelinearizationKeySet) AddRelinearizationKey(rlk *RelinearizationKey) {
	rlkSet.Value[rlk.ID] = rlk
	rlkSet.PolyRPool1[rlk.ID] = rlkSet.params.RingR().NewPoly()
	rlkSet.PolyRPool2[rlk.ID] = rlkSet.params.RingR().NewPoly()

	rlkSet.HoistPool1[0].Value[rlk.ID] = mkrlwe.NewSwitchingKey(rlkSet.params.Parameters)
	rlkSet.HoistPool1[1].Value[rlk.ID] = mkrlwe.NewSwitchingKey(rlkSet.params.Parameters)

	rlkSet.HoistPool2[0].Value[rlk.ID] = mkrlwe.NewSwitchingKey(rlkSet.params.Parameters)
	rlkSet.HoistPool2[1].Value[rlk.ID] = mkrlwe.NewSwitchingKey(rlkSet.params.Parameters)
}

// DelRelinearizationKey delete publickey of given id from SecretKeySet
func (rlkSet *RelinearizationKeySet) DelRelinearizationKey(id string) {
	delete(rlkSet.Value, id)
}

// GetRelinearizationKey returns a publickey of given id from RelinearizationKeySet
func (rlkSet *RelinearizationKeySet) GetRelinearizationKey(id string) *RelinearizationKey {
	ret, in := rlkSet.Value[id]

	if !in {
		panic("cannot GetRelinearizationKey: there is no relinearization key with given id")
	}

	return ret
}
