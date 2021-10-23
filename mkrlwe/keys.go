package mkrlwe

import "github.com/ldsec/lattigo/v2/rlwe"
import "math"

// SecretKeySet is a type for generic Multikey RLWE secret keys.
type SecretKey struct {
	rlwe.SecretKey
}

// SecretKeySet is a type for a set of multikey RLWE secret keys.
type SecretKeySet struct {
	Value map[string]*SecretKey
}

// PublicKey is a type for generic RLWE public keys.
type PublicKey struct {
	rlwe.PublicKey
}

// PublicKeySet is a type for a set of multikey RLWE public keys.
type PublicKeySet struct {
	Value map[string]*PublicKey
}

// RelinearizationKey is a type for generic RLWE public relinearization keys.
// It consists of three polynomial vectors
type RelinearizationKey struct {
	Value [3]*PolyQPVector
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
func (skSet *SecretKeySet) AddSecretKey(id string, sk *SecretKey) {
	skSet.Value[id] = sk
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
func (pkSet *PublicKeySet) AddPublicKey(id string, pk *PublicKey) {
	pkSet.Value[id] = pk
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
func (rlkSet *RelinearizationKeySet) AddRelinearizationKey(id string, rlk *RelinearizationKey) {
	rlkSet.Value[id] = rlk
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
func NewSecretKey(params Parameters) *SecretKey {
	sk := new(SecretKey)
	sk.Value = params.RingQP().NewPoly()
	return sk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params Parameters) *PublicKey {
	pk := new(PublicKey)
	pk.Value[0] = params.RingQP().NewPoly()
	pk.Value[1] = params.RingQP().NewPoly()
	return pk
}

// NewRelinearizationKey returns a new RelinearizationKey with zero values.
func NewRelinearizationKey(params Parameters, levelQ, levelP int) *RelinearizationKey {
	decompSize := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))
	r := RingQP{*params.RingQP()}
	rlk := new(RelinearizationKey)

	rlk.Value[0] = r.NewPolyVector(decompSize)
	rlk.Value[1] = r.NewPolyVector(decompSize)
	rlk.Value[2] = r.NewPolyVector(decompSize)

	return rlk
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk *SecretKey) CopyNew() *SecretKey {
	if sk == nil {
		return nil
	}

	ret := new(SecretKey)
	ret.Value = sk.Value.CopyNew()

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

	return ret
}

// CopyNew creates a deep copy of the receiver RelinearizationKey and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	if rlk == nil {
		return nil
	}

	ret := new(RelinearizationKey)
	ret.Value[0] = rlk.Value[0].CopyNew()
	ret.Value[1] = rlk.Value[1].CopyNew()
	ret.Value[2] = rlk.Value[2].CopyNew()

	return ret
}
