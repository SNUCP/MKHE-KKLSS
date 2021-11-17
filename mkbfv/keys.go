package mkbfv

import "mk-lattigo/mkrlwe"
import "github.com/ldsec/lattigo/v2/rlwe"

// SecretKeySet is a type for generic Multikey RLWE secret keys.
type SecretKey struct {
	*mkrlwe.SecretKey
	ValueQP    rlwe.PolyQP
	ValueQMulP rlwe.PolyQP
	ValueRP    rlwe.PolyQP
	ID         string
}

// SecretKeySet is a type for a set of multikey RLWE secret keys.
type SecretKeySet struct {
	Value map[string]*SecretKey
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
