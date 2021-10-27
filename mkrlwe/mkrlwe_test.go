package mkrlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	//"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
var TestParams = []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218, rlwe.TestPN14QP438, rlwe.TestPN15QP880}

func testString(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestMKRLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	/*
		  if testing.Short() {
				defaultParams = TestParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
			}
	*/

	if *flagParamString != "" {
		var jsonParams rlwe.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []rlwe.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[:] {
		params, err := rlwe.NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		mkparams := NewParameters(params)
		kgen := NewKeyGenerator(mkparams)

		testGenKeyPair(kgen, t)
		testSwitchKeyGen(kgen, t)
		testRelinKeyGen(kgen, t)

		testEncryptor(kgen, t)
		testDecryptor(kgen, t)

		testInternalProduct(kgen, t)

	}

}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients
func log2OfInnerSum(level int, ringQ *ring.Ring, poly *ring.Poly) (logSum int) {
	sumRNS := make([]uint64, level+1)
	var sum uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		qiHalf := qi >> 1
		coeffs := poly.Coeffs[i]
		sum = 0

		for j := 0; j < ringQ.N; j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = ring.CRed(sum+qi-v, qi)
			} else {
				sum = ring.CRed(sum+v, qi)
			}
		}

		sumRNS[i] = sum
	}

	var smallNorm = true
	for i := 1; i < level+1; i++ {
		smallNorm = smallNorm && (sumRNS[0] == sumRNS[i])
	}

	if !smallNorm {
		var qi uint64
		var crtReconstruction *big.Int

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ring.NewUint(1)

		for i := 0; i < level+1; i++ {

			qi = ringQ.Modulus[i]
			QiB.SetUint64(qi)

			modulusBigint.Mul(modulusBigint, QiB)

			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(ringQ.ModulusBigint, QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)

			sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNS[i]), crtReconstruction))
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNS[0])
	}

	return
}

func testGenKeyPair(kgen *KeyGenerator, t *testing.T) {

	params := kgen.params

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PKGen/"), func(t *testing.T) {
		id := "user"
		sk, pk := kgen.GenKeyPair(id)

		ringQ := params.RingQ()
		ringP := params.RingP()
		ringQP := params.RingQP()

		levelQ, levelP := params.QCount()-1, params.PCount()-1

		// [-as + e] + [as]

		ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, sk.Value, pk.Value[1], pk.Value[0])
		ringQP.InvNTTLvl(levelQ, levelP, pk.Value[0], pk.Value[0])

		log2Bound := bits.Len64(uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelQ, ringQ, pk.Value[0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelP, ringP, pk.Value[0].P))
	})
}

func testEncryptor(kgen *KeyGenerator, t *testing.T) {

	params := kgen.params

	t.Run(testString(params, "Encrypt/Pk/MaxLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}

		var user1 string = "tetsUser1"
		var users []string = []string{user1}

		sk, pk := kgen.GenKeyPair(user1)
		ringQ := params.RingQ()

		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value0)
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value0, ciphertext.Value0)
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value0))
	})

	t.Run(testString(params, "Encrypt/Pk/MinLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}

		var user1 string = "tetsUser1"
		var users []string = []string{user1}

		sk, pk := kgen.GenKeyPair(user1)
		ringQ := params.RingQ()

		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value0)
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value0, ciphertext.Value0)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value0))
	})
}

func testSwitchKeyGen(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2

	params := kgen.params
	t.Run(testString(params, "SWKGen/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}
		id := "User"
		ringQ := params.RingQ()
		ringQP := params.RingQP()
		sk := kgen.GenSecretKey(id)

		levelQ, levelP := params.QCount()-1, params.PCount()-1
		beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

		//generate sg
		sg := NewSwitchingKey(params, id)
		kgen.GenSwitchingKey(sk, sg)

		//tmp = P*s
		var pBigInt *big.Int
		if levelP == kgen.params.PCount()-1 {
			pBigInt = kgen.params.RingP().ModulusBigint
		} else {
			P := kgen.params.RingP().Modulus
			pBigInt = new(big.Int).SetUint64(P[0])
			for i := 1; i < levelP+1; i++ {
				pBigInt.Mul(pBigInt, ring.NewUint(P[i]))
			}
		}
		tmp := ringQP.NewPoly()
		ringQ.MulScalarBigint(sk.Value.Q, pBigInt, tmp.Q)
		ringQP.InvMFormLvl(levelQ, levelP, tmp, tmp)

		for i := 0; i < beta; i++ {
			ringQP.SubLvl(levelQ, levelP, tmp, sg.Value[i], tmp)
		}

		ringQP.InvNTTLvl(levelQ, levelP, tmp, tmp)

		//check errors
		log2Bound := bits.Len64(uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()*len(sg.Value)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(tmp.Q.Level(), params.RingQ(), tmp.Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(tmp.P.Level(), params.RingP(), tmp.P))
	})
}

func testRelinKeyGen(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2
	params := kgen.params
	t.Run(testString(params, "RLKGen/"), func(t *testing.T) {
		id := "User"
		ringQ := params.RingQ()
		ringP := params.RingP()
		ringQP := params.RingQP()
		s := kgen.GenSecretKey(id)
		r := kgen.GenSecretKey(id)

		levelQ, levelP := params.QCount()-1, params.PCount()-1
		beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))

		a := params.CRS[0]
		//u := params.CRS[1]

		//generate sg, rg
		sg := NewSwitchingKey(params, id)
		kgen.GenSwitchingKey(s, sg)

		rg := NewSwitchingKey(params, id)
		kgen.GenSwitchingKey(r, rg)

		//gen rlk
		rlk := kgen.GenRelinearizationKey(s, r)

		//check b=sa+e
		b := rlk.Value[0]
		for i := 0; i < beta; i++ {
			ringQP.InvMFormLvl(levelQ, levelP, b.Value[i], b.Value[i])

			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a[i], s.Value, b.Value[i])
			ringQP.InvNTTLvl(levelQ, levelP, b.Value[i], b.Value[i])

			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelQ, ringQ, b.Value[i].Q))
			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelP, ringP, b.Value[i].P))
		}

		//check d = ra + sg + e
		d := rlk.Value[1]
		for i := 0; i < beta; i++ {
			ringQP.InvMFormLvl(levelQ, levelP, d.Value[i], d.Value[i])

			ringQP.MulCoeffsMontgomeryAndSubLvl(levelQ, levelP, a[i], r.Value, d.Value[i])
			ringQP.SubLvl(levelQ, levelP, d.Value[i], sg.Value[i], d.Value[i])
			ringQP.InvNTTLvl(levelQ, levelP, d.Value[i], d.Value[i])

			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelQ, ringQ, d.Value[i].Q))
			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelP, ringP, d.Value[i].P))
		}

	})
}

func testDecryptor(kgen *KeyGenerator, t *testing.T) {
	params := kgen.params

	var user1 string = "tetsUser1"
	var users []string = []string{user1}

	sk, pk := kgen.GenKeyPair(user1)
	ringQ := params.RingQ()
	encryptor := NewEncryptor(params, pk)
	decryptor := NewDecryptor(params)

	skSet := NewSecretKeySet()
	skSet.AddSecretKey(sk)

	t.Run(testString(params, "Decrypt/MaxLevel/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})
}

func testInternalProduct(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2

	params := kgen.params

	t.Run(testString(params, "InternalProductMaxLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		var id string = "tetsUser1"
		var idList []string = []string{id}
		ringQ := params.RingQ()
		ringQP := params.RingQP()
		levelQ := params.QCount() - 1
		levelP := params.PCount() - 1
		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		plaintext := rlwe.NewPlaintext(params.Parameters, levelQ)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, idList, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)

		//generate sg
		sg := NewSwitchingKey(params, id)
		kgen.GenSwitchingKey(sk, sg)

		for i := 0; i < len(sg.Value); i++ {
			ringQP.MFormLvl(levelQ, levelP, sg.Value[i], sg.Value[i])
		}

		//tmp = Inter(c, sg)
		ks := NewKeySwitcher(params)
		tmp := ringQ.NewPolyLvl(ciphertext.Level())
		ks.InternalProduct(ciphertext.Level(), ciphertext.Value0, sg, tmp)

		//tmp2 = c*s
		ringQ.MulCoeffsMontgomeryAndSubLvl(ciphertext.Level(), ciphertext.Value0, sk.Value.Q, tmp)
		ringQ.InvNTTLvl(ciphertext.Level(), tmp, tmp)

		//check errors
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(tmp.Level(), params.RingQ(), tmp))
	})

	t.Run(testString(params, "InternalProductMinLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		var id string = "tetsUser1"
		var idList []string = []string{id}
		ringQP := params.RingQP()
		ringQ := params.RingQ()
		levelQ := params.QCount() - 1
		levelP := params.PCount() - 1
		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		plaintext := rlwe.NewPlaintext(params.Parameters, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, idList, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)

		//generate sg
		sg := NewSwitchingKey(params, id)
		kgen.GenSwitchingKey(sk, sg)

		for i := 0; i < len(sg.Value); i++ {
			ringQP.MFormLvl(levelQ, levelP, sg.Value[i], sg.Value[i])
		}

		//tmp = P*s
		ks := NewKeySwitcher(params)
		tmp := ringQ.NewPolyLvl(ciphertext.Level())
		ks.InternalProduct(ciphertext.Level(), ciphertext.Value0, sg, tmp)

		//tmp2 = c*s
		ringQ.MulCoeffsMontgomeryAndSubLvl(ciphertext.Level(), ciphertext.Value0, sk.Value.Q, tmp)
		ringQ.InvNTTLvl(ciphertext.Level(), tmp, tmp)

		//check errors
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(tmp.Level(), params.RingQ(), tmp))
	})
}
