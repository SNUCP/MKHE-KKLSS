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
	if testing.Short() {
		defaultParams = TestParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}

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
		testRelinKeyGen(kgen, t)
		testEncryptor(*kgen, t)
		testDecryptor(*kgen, t)
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

		// [-as + e] + [as]
		params.RingQP().MulCoeffsMontgomeryAndAddLvl(sk.Value.Q.Level(), sk.Value.P.Level(), sk.Value, pk.Value[1], pk.Value[0])
		params.RingQP().InvNTTLvl(sk.Value.Q.Level(), sk.Value.P.Level(), pk.Value[0], pk.Value[0])

		log2Bound := bits.Len64(uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), params.RingQ(), pk.Value[0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].P.Level(), params.RingP(), pk.Value[0].P))
	})
}

func testRelinKeyGen(kgen *KeyGenerator, t *testing.T) {

	// Checks that switching keys are en encryption under the output key
	// of the RNS decomposition of the input key by
	// 1) Decrypting the RNS decomposed input key
	// 2) Reconstructing the key
	// 3) Checking that the difference with the input key has a small norm

	params := kgen.params
	t.Run(testString(params, "RLKGen/"), func(t *testing.T) {

		id := "user"

		ringQP := RingQP{*params.RingQP()}
		sk := kgen.GenSecretKey(id)
		levelQ, levelP := params.QCount()-1, params.PCount()-1
		beta := int(math.Ceil(float64(levelQ+1) / float64(levelP+1)))
		// Generates RelinearizationKey
		rlk, r := kgen.GenRelinearizationKey(sk)
		tmp := ringQP.NewPolyVector(beta)

		// Decrypts
		// b - sa
		ringQP.MulPolyMontgomeryLvl(levelQ, levelP, params.CRS[0], &sk.Value, tmp)
		ringQP.SubLvl(levelQ, levelP, rlk.Value[0], tmp, rlk.Value[0])

		ringQP.InvNTTLvl(levelQ, levelP, rlk.Value[0], rlk.Value[0])
		ringQP.InvMFormLvl(levelQ, levelP, rlk.Value[0], rlk.Value[0])

		// d - ra - sg
		ringQP.MulPolyMontgomeryLvl(levelQ, levelP, params.CRS[0], &r.Value, tmp)
		ringQP.SubLvl(levelQ, levelP, rlk.Value[1], tmp, rlk.Value[1])
		ringQP.MulPolyMontgomeryLvl(levelQ, levelP, kgen.gadgetVector, &sk.Value, tmp)
		ringQP.SubLvl(levelQ, levelP, rlk.Value[1], tmp, rlk.Value[1])

		ringQP.InvNTTLvl(levelQ, levelP, rlk.Value[1], rlk.Value[1])
		ringQP.InvMFormLvl(levelQ, levelP, rlk.Value[1], rlk.Value[1])

		// v - su - rg
		ringQP.MulPolyMontgomeryLvl(levelQ, levelP, params.CRS[1], &sk.Value, tmp)
		ringQP.SubLvl(levelQ, levelP, rlk.Value[2], tmp, rlk.Value[2])
		ringQP.MulPolyMontgomeryLvl(levelQ, levelP, kgen.gadgetVector, &r.Value, tmp)
		ringQP.SubLvl(levelQ, levelP, rlk.Value[2], tmp, rlk.Value[2])

		ringQP.InvNTTLvl(levelQ, levelP, rlk.Value[2], rlk.Value[2])
		ringQP.InvMFormLvl(levelQ, levelP, rlk.Value[2], rlk.Value[2])

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys

		for j := 0; j < beta; j++ {

			log2Bound := bits.Len64(uint64(math.Floor(rlwe.DefaultSigma*6)) * uint64(params.N()))

			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[0].Value[j].Q.Level(), params.RingQ(), rlk.Value[0].Value[j].Q))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[0].Value[j].P.Level(), params.RingP(), rlk.Value[0].Value[j].P))

			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[1].Value[j].Q.Level(), params.RingQ(), rlk.Value[1].Value[j].Q))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[1].Value[j].P.Level(), params.RingP(), rlk.Value[1].Value[j].P))

			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[2].Value[j].Q.Level(), params.RingQ(), rlk.Value[2].Value[j].Q))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(rlk.Value[2].Value[j].P.Level(), params.RingP(), rlk.Value[2].Value[j].P))

		}

	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.params

	var user1 string = "tetsUser1"
	var users []string = []string{user1}

	sk, pk := kgen.GenKeyPair(user1)

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/Slow/MaxLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := Plaintext{*rlwe.NewPlaintext(params.Parameters, params.MaxLevel())}
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params.Parameters, users, 1, plaintext.Level())
		encryptor.Encrypt(&plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value0)
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value0, ciphertext.Value0)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value0))
	})

	t.Run(testString(params, "Encrypt/Pk/Slow/MinLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := Plaintext{*rlwe.NewPlaintext(params.Parameters, params.MaxLevel())}
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params.Parameters, users, 1, plaintext.Level())
		encryptor.Encrypt(&plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value0)
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value0, ciphertext.Value0)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value0))
	})
}

func testDecryptor(kgen KeyGenerator, t *testing.T) {
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
		plaintext := Plaintext{*rlwe.NewPlaintext(params.Parameters, params.MaxLevel())}
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params.Parameters, users, 1, plaintext.Level())
		encryptor.Encrypt(&plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, &plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel/"), func(t *testing.T) {
		plaintext := Plaintext{*rlwe.NewPlaintext(params.Parameters, params.MaxLevel())}
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params.Parameters, users, 1, plaintext.Level())
		encryptor.Encrypt(&plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, &plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})
}
