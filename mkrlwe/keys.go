package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "math"

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

// PublicKeySet is a type for a set of multikey RLWE public keys.
type PublicKeySet struct {
	Value map[string]*PublicKey
}

// RelinearizationKey is a type for generic RLWE public relinearization keys.
// It consists of three polynomial vectors
type RelinearizationKey struct {
	Value [3][]rlwe.PolyQP
	ID    string
}

// RelinearizationKeySet is a type for a set of multikey RLWE relinearization keys.
type RelinearizationKeySet struct {
	Value map[string]*RelinearizationKey
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

// NewRelinearizationKey returns a new RelinearizationKey with zero values.
func NewRelinearizationKey(params Parameters, levelQ, levelP int, id string) *RelinearizationKey {
	beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))
	rlk := new(RelinearizationKey)
	rlk.Value[0] = make([]rlwe.PolyQP, beta)
	rlk.Value[1] = make([]rlwe.PolyQP, beta)
	rlk.Value[2] = make([]rlwe.PolyQP, beta)

	ringQP := params.RingQP()

	for i := 0; i < beta; i++ {
		rlk.Value[0][i] = ringQP.NewPoly()
		rlk.Value[1][i] = ringQP.NewPoly()
		rlk.Value[2][i] = ringQP.NewPoly()

	}

	rlk.ID = id

	return rlk
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

// CopyNew creates a deep copy of the receiver RelinearizationKey and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	if rlk == nil {
		return nil
	}

	beta := len(rlk.Value[0])
	ret := new(RelinearizationKey)

	ret.Value[0] = make([]rlwe.PolyQP, beta)
	ret.Value[1] = make([]rlwe.PolyQP, beta)
	ret.Value[2] = make([]rlwe.PolyQP, beta)

	for i := 0; i < beta; i++ {
		ret.Value[0][i] = rlk.Value[0][i].CopyNew()
		ret.Value[1][i] = rlk.Value[1][i].CopyNew()
		ret.Value[2][i] = rlk.Value[2][i].CopyNew()
	}
	ret.ID = rlk.ID

	return ret
}
