package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"

// SecretKeySet is a type for generic Multikey RLWE secret keys.
type SecretKey struct {
	rlwe.SecretKey
	ID string
}

// SecretKeySet is a type for a set of multikey RLWE secret keys.
type SecretKeySet struct {
	Value map[string]*SecretKey
}

// PublicKey is a type for generic RLWE public keys.
type PublicKey struct {
	rlwe.PublicKey
	ID string
}

// SwitchingKey is a type for generic RLWE switching keys.
type SwitchingKey struct {
	Value []rlwe.PolyQP
}

// PublicKeySet is a type for a set of multikey RLWE public keys.
type PublicKeySet struct {
	Value map[string]*PublicKey
}

// RelinearizationKey is a type for generic RLWE public relinearization keys.
// It consists of three polynomial vectors
type RelinearizationKey struct {
	Value [3]*SwitchingKey
	ID    string
}

// RotationKey is a type for storing generic RLWE public rotation keys.
type RotationKey struct {
	Value  *SwitchingKey
	ID     string
	RotIdx int
}

// RelinearizationKeySet is a type for a set of multikey RLWE relinearization keys.
type RelinearizationKeySet struct {
	Value map[string]*RelinearizationKey
}

//RotationKeysSet is a type for a set of multikey RLWE rotation keys.
type RotationKeySet struct {
	Value map[string]map[int]*RotationKey
}

// NewSecretKeySet returns a new empty SecretKeySet
func NewSecretKeySet() *SecretKeySet {
	skSet := new(SecretKeySet)
	skSet.Value = make(map[string]*SecretKey)
	return skSet
}

// AddSecretKey insert new secretkey into SecretKeySet with its id
func (skSet *SecretKeySet) AddSecretKey(sk *SecretKey) {
	skSet.Value[sk.ID] = sk
}

// DelSecretKey delete secretkey of given id from SecretKeySet
func (skSet *SecretKeySet) DelSecretKey(id string) {
	delete(skSet.Value, id)
}

// GetSecretKey returns a secretkey of given id from SecretKeySet
func (skSet *SecretKeySet) GetSecretKey(id string) *SecretKey {
	ret, in := skSet.Value[id]

	if !in {
		panic("cannot GetPublicKey: there is no public key with given id")
	}
	return ret
}

// NewPublicKeySet returns a new empty PublicKeySet
func NewPublicKeyKeySet() *PublicKeySet {
	pkSet := new(PublicKeySet)
	pkSet.Value = make(map[string]*PublicKey)
	return pkSet
}

// AddPublicKey insert new publickey into PublicKeySet with its id
func (pkSet *PublicKeySet) AddPublicKey(pk *PublicKey) {
	pkSet.Value[pk.ID] = pk
}

// DelPublicKey delete publickey of given id from SecretKeySet
func (pkSet *PublicKeySet) DelPublicKey(id string) {
	delete(pkSet.Value, id)
}

// GetPublicKey returns a publickey of given id from PublicKeySet
func (pkSet *PublicKeySet) GetPublicKey(id string) *PublicKey {
	ret, in := pkSet.Value[id]

	if !in {
		panic("cannot GetPublicKey: there is no public key with given id")
	}

	return ret
}

// NewRotationKeysSet returns a new empty RotationKeysSet
func NewRotationKeysSet() *RotationKeySet {
	rotSet := new(RotationKeySet)
	rotSet.Value = make(map[string]map[int]*RotationKey)

	return rotSet
}

// AddRotationKeys insert new rotation keys into RotationKeysSet with its id
func (rkSet *RotationKeySet) AddRotationKey(rk *RotationKey) {

	_, ok := rkSet.Value[rk.ID]

	if !ok {
		rkSet.Value[rk.ID] = make(map[int]*RotationKey)
	}

	rkSet.Value[rk.ID][rk.RotIdx] = rk

}

// DelRotationKeys delete rotation keys of given id from RotationKeysSet
func (rkSet *RotationKeySet) DelRotationKey(id string, rotidx int) {
	delete(rkSet.Value[id], rotidx)
}

// GetRotationKeys returns a rotation keys of given id from RotationKeysSet
func (rkSet *RotationKeySet) GetRotationKey(id string, rotidx int) *RotationKey {
	_, in := rkSet.Value[id]
	if !in {
		panic("cannot GetRotationKeys: there is no rotation key with given id")
	}

	_, in = rkSet.Value[id][rotidx]
	if !in {
		panic("cannot GetRotationKeys: there is no rotation key with given id")
	}
	return rkSet.Value[id][rotidx]
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

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params Parameters, id string) *SecretKey {
	sk := new(SecretKey)
	sk.Value = params.RingQP().NewPoly()
	sk.ID = id
	return sk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params Parameters, id string) *PublicKey {
	pk := new(PublicKey)
	pk.Value[0] = params.RingQP().NewPoly()
	pk.Value[1] = params.RingQP().NewPoly()
	pk.ID = id
	return pk
}

//
func NewSwitchingKey(params Parameters) *SwitchingKey {
	beta := params.Beta(params.QCount() - 1)
	ringQP := params.RingQP()
	swk := new(SwitchingKey)
	swk.Value = make([]rlwe.PolyQP, beta)
	for i := 0; i < beta; i++ {
		swk.Value[i] = ringQP.NewPoly()
	}

	return swk
}

// NewRelinearizationKey returns a new RelinearizationKey with zero values.
func NewRelinearizationKey(params Parameters, id string) *RelinearizationKey {
	rlk := new(RelinearizationKey)
	rlk.Value[0] = NewSwitchingKey(params)
	rlk.Value[1] = NewSwitchingKey(params)
	rlk.Value[2] = NewSwitchingKey(params)

	rlk.ID = id

	return rlk
}

func NewRotationKey(params Parameters, rotidx int, id string) *RotationKey {
	rk := new(RotationKey)
	rk.ID = id
	rk.RotIdx = rotidx
	rk.Value = NewSwitchingKey(params)

	return rk
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk *SecretKey) CopyNew() *SecretKey {
	if sk == nil {
		return nil
	}

	ret := new(SecretKey)
	ret.Value = sk.Value.CopyNew()
	ret.ID = sk.ID

	return ret
}

// CopyNew creates a deep copy of the receiver PublicKey and returns it.
func (pk *PublicKey) CopyNew() *PublicKey {
	if pk == nil {
		return nil
	}

	ret := new(PublicKey)
	ret.Value[0] = pk.Value[0].CopyNew()
	ret.Value[1] = pk.Value[1].CopyNew()
	ret.ID = pk.ID

	return ret
}
