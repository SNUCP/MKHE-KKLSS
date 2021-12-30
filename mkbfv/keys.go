package mkbfv

import "mk-lattigo/mkrlwe"

type RelinearizationKey struct {
	Value [2]*mkrlwe.RelinearizationKey
	ID    string
}

// RelinearizationKeySet is a type for a set of multikey BFV relinearization keys.
type RelinearizationKeySet struct {
	Value map[string]*RelinearizationKey
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
func NewRelinearizationKeyKeySet() *RelinearizationKeySet {
	rlkSet := new(RelinearizationKeySet)
	rlkSet.Value = make(map[string]*RelinearizationKey)
	return rlkSet
}

// AddRelinearizationKey insert new publickey into RelinearizationKeySet with its id
func (rlkSet *RelinearizationKeySet) AddRelinearizationKey(rlk *RelinearizationKey) {
	rlkSet.Value[rlk.ID] = rlk
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
