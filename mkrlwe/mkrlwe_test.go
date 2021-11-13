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

var PN16QP1761 = rlwe.ParametersLiteral{
	LogN: 16,
	Q: []uint64{0x80000000080001, 0x2000000a0001, 0x2000000e0001, 0x1fffffc20001, // 55 + 33 x 45
		0x200000440001, 0x200000500001, 0x200000620001, 0x1fffff980001,
		0x2000006a0001, 0x1fffff7e0001, 0x200000860001, 0x200000a60001,
		0x200000aa0001, 0x200000b20001, 0x200000c80001, 0x1fffff360001,
		0x200000e20001, 0x1fffff060001, 0x200000fe0001, 0x1ffffede0001,
		0x1ffffeca0001, 0x1ffffeb40001, 0x200001520001, 0x1ffffe760001,
		0x2000019a0001, 0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
		0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001, 0x200002480001,
		0x1ffffdb60001, 0x200002560001},
	P:     []uint64{0x80000000440001, 0x7fffffffba0001, 0x80000000500001, 0x7fffffffaa0001}, // 4 x 55
	Sigma: rlwe.DefaultSigma,
}

//var TestParams = []rlwe.ParametersLiteral{PN16QP1761}

var TestParams = []rlwe.ParametersLiteral{rlwe.TestPN12QP109, rlwe.TestPN13QP218, rlwe.TestPN14QP438, rlwe.TestPN15QP880}

//PN16QP1761}

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

		gamma := 2

		if params.PCount() < gamma {
			continue
		}

		mkparams := NewParameters(params, gamma)
		kgen := NewKeyGenerator(mkparams)

		testGenKeyPair(kgen, t)
		testSwitchKeyGen(kgen, t)
		testRelinKeyGen(kgen, t)

		testEncryptor(kgen, t)
		testDecryptor(kgen, t)

		testDecompose(kgen, t)
		testInternalProduct(kgen, t)
		testHadamardProduct(kgen, t)
		testRelinearize(kgen, t)

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

	smallNorm = false
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
		users := NewIDSet()
		users.Add(user1)

		sk, pk := kgen.GenKeyPair(user1)
		ringQ := params.RingQ()

		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value["0"])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value["0"], ciphertext.Value["0"])
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value["0"]))
	})

	t.Run(testString(params, "Encrypt/Pk/MinLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}

		var user1 string = "tetsUser1"
		users := NewIDSet()
		users.Add(user1)

		sk, pk := kgen.GenKeyPair(user1)
		ringQ := params.RingQ()

		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[user1], sk.Value.Q, ciphertext.Value["0"])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value["0"], ciphertext.Value["0"])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value["0"]))
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
		beta := params.Beta(levelQ)

		//generate sg
		sg := NewSwitchingKey(params)
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

		for i := 0; i < beta; i++ {
			ringQP.SubLvl(levelQ, levelP, tmp, sg.Value[i], tmp)
		}

		ringQP.InvMFormLvl(levelQ, levelP, tmp, tmp)
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
		beta := params.Beta(levelQ)

		a := params.CRS[0]
		u := params.CRS[1]

		//generate sg, rg
		sg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(s, sg)

		rg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(r, rg)

		//gen rlk
		rlk := kgen.GenRelinearizationKey(s, r)

		//check b=-sa+e
		b := rlk.Value[0]
		for i := 0; i < beta; i++ {

			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, a.Value[i], s.Value, b.Value[i])

			ringQP.InvMFormLvl(levelQ, levelP, b.Value[i], b.Value[i])
			ringQP.InvNTTLvl(levelQ, levelP, b.Value[i], b.Value[i])

			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelQ, ringQ, b.Value[i].Q))
			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelP, ringP, b.Value[i].P))
		}

		//check d = -ra + sg + e
		d := rlk.Value[1]
		for i := 0; i < beta; i++ {

			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, a.Value[i], r.Value, d.Value[i])
			ringQP.SubLvl(levelQ, levelP, d.Value[i], sg.Value[i], d.Value[i])

			ringQP.InvMFormLvl(levelQ, levelP, d.Value[i], d.Value[i])
			ringQP.InvNTTLvl(levelQ, levelP, d.Value[i], d.Value[i])

			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelQ, ringQ, d.Value[i].Q))
			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelP, ringP, d.Value[i].P))
		}

		//check v = -su -rg + e
		v := rlk.Value[2]
		for i := 0; i < beta; i++ {

			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, u.Value[i], s.Value, v.Value[i])
			ringQP.AddLvl(levelQ, levelP, v.Value[i], rg.Value[i], v.Value[i])

			ringQP.InvMFormLvl(levelQ, levelP, v.Value[i], v.Value[i])
			ringQP.InvNTTLvl(levelQ, levelP, v.Value[i], v.Value[i])

			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelQ, ringQ, v.Value[i].Q))
			require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(levelP, ringP, v.Value[i].P))
		}

	})
}

func testDecryptor(kgen *KeyGenerator, t *testing.T) {
	params := kgen.params

	var user1 string = "tetsUser1"
	users := NewIDSet()
	users.Add(user1)

	sk, pk := kgen.GenKeyPair(user1)
	ringQ := params.RingQ()
	encryptor := NewEncryptor(params)
	decryptor := NewDecryptor(params)

	skSet := NewSecretKeySet()
	skSet.AddSecretKey(sk)

	t.Run(testString(params, "Decrypt/MaxLevel/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		ciphertext := NewCiphertext(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "DecryptNTT/MaxLevel/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)
		decryptor.Decrypt(ciphertext, skSet, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Decrypt/Multikey/"), func(t *testing.T) {
		plaintext := rlwe.NewPlaintext(params.Parameters, params.MaxLevel())
		plaintext.Value.IsNTT = true

		user1 := "user1"
		user2 := "user2"
		idset1 := NewIDSet()
		idset2 := NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)
		idset := idset1.Union(idset2)

		sk1, pk1 := kgen.GenKeyPair(user1)
		sk2, pk2 := kgen.GenKeyPair(user2)

		skSet.AddSecretKey(sk1)
		skSet.AddSecretKey(sk2)

		level := plaintext.Level()

		ct1 := NewCiphertextNTT(params, idset1, level)
		ct2 := NewCiphertextNTT(params, idset2, level)
		ctOut := NewCiphertextNTT(params, idset, level)

		encryptor.Encrypt(plaintext, pk1, ct1)
		encryptor.Encrypt(plaintext, pk2, ct2)

		ringQ.AddLvl(level, ct1.Value["0"], ct2.Value["0"], ctOut.Value["0"])
		ctOut.Value[user1].Copy(ct1.Value[user1])
		ctOut.Value[user2].Copy(ct2.Value[user2])

		decryptor.Decrypt(ct1, skSet, plaintext)
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ct1.Level(), ringQ, plaintext.Value))

		decryptor.Decrypt(ct2, skSet, plaintext)
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ct2.Level(), ringQ, plaintext.Value))

		decryptor.Decrypt(ctOut, skSet, plaintext)
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(ctOut.Level(), ringQ, plaintext.Value))

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
		users := NewIDSet()
		users.Add(id)

		levelQ := params.QCount() - 1

		ringQ := params.RingQ()
		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		plaintext := rlwe.NewPlaintext(params.Parameters, levelQ)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)

		//generate sg
		sg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(sk, sg)

		//tmp = P*s
		ks := NewKeySwitcher(params)
		tmp := ringQ.NewPolyLvl(ciphertext.Level())
		ks.InternalProduct(ciphertext.Level(), ciphertext.Value["0"], sg, tmp)

		//tmp2 = c*s
		ringQ.MulCoeffsMontgomeryAndSubLvl(ciphertext.Level(), ciphertext.Value["0"], sk.Value.Q, tmp)
		ringQ.InvNTTLvl(ciphertext.Level(), tmp, tmp)

		//check errors
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(tmp.Level(), params.RingQ(), tmp))
	})
}

func testDecompose(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2

	params := kgen.params

	t.Run(testString(params, "DecomposeMaxLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		var id string = "tetsUser1"
		users := NewIDSet()
		users.Add(id)

		ringQ := params.RingQ()
		ringQP := params.RingQP()
		levelQ := params.QCount() - 1
		levelP := params.PCount() - 1
		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		plaintext := rlwe.NewPlaintext(params.Parameters, levelQ-1)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)

		beta := params.Beta(levelQ)

		//generate sg
		sg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(sk, sg)

		//generate cd
		ks := NewKeySwitcher(params)
		cd := NewSwitchingKey(params)
		ks.Decompose(ciphertext.Level(), ciphertext.Value["0"], cd)

		//internal product on tmp
		tmp := ringQP.NewPolyLvl(ciphertext.Level(), levelP)
		for i := 0; i < beta; i++ {
			ringQP.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), levelP, cd.Value[i], sg.Value[i], tmp)
		}
		ks.Baseconverter.ModDownQPtoQNTT(ciphertext.Level(), levelP, tmp.Q, tmp.P, tmp.Q)

		//tmp2 = c*s
		ringQ.MulCoeffsMontgomeryAndSubLvl(ciphertext.Level(), ciphertext.Value["0"], sk.Value.Q, tmp.Q)
		ringQ.InvNTTLvl(ciphertext.Level(), tmp.Q, tmp.Q)

		//check errors
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(tmp.Q.Level(), params.RingQ(), tmp.Q))
	})

	t.Run(testString(params, "DecomposeMinLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		var id string = "tetsUser1"
		users := NewIDSet()
		users.Add(id)

		ringQ := params.RingQ()
		ringQP := params.RingQP()
		levelP := params.PCount() - 1
		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		plaintext := rlwe.NewPlaintext(params.Parameters, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params)
		ciphertext := NewCiphertextNTT(params, users, plaintext.Level())
		encryptor.Encrypt(plaintext, pk, ciphertext)

		beta := params.Beta(0)

		//generate sg
		sg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(sk, sg)

		//generate cd
		ks := NewKeySwitcher(params)
		cd := NewSwitchingKey(params)
		ks.Decompose(ciphertext.Level(), ciphertext.Value["0"], cd)

		//internal product on tmp
		tmp := ringQP.NewPolyLvl(ciphertext.Level(), levelP)
		for i := 0; i < beta; i++ {
			ringQP.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), levelP, cd.Value[i], sg.Value[i], tmp)
		}
		ks.Baseconverter.ModDownQPtoQNTT(ciphertext.Level(), levelP, tmp.Q, tmp.P, tmp.Q)

		//tmp2 = c*s
		ringQ.MulCoeffsMontgomeryAndSubLvl(ciphertext.Level(), ciphertext.Value["0"], sk.Value.Q, tmp.Q)
		ringQ.InvNTTLvl(ciphertext.Level(), tmp.Q, tmp.Q)

		//check errors
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(tmp.Q.Level(), params.RingQ(), tmp.Q))
	})
}

func testRelinearize(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2

	params := kgen.params

	t.Run(testString(params, "RelinMaxLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		ringQ := params.RingQ()
		levelQ := params.QCount() - 1

		level := levelQ

		pt := rlwe.NewPlaintext(params.Parameters, level)
		pt.Value.IsNTT = true

		user1 := "user1"
		user2 := "user2"
		idset1 := NewIDSet()
		idset2 := NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)
		idset := idset1.Union(idset2)

		sk1, pk1 := kgen.GenKeyPair(user1)
		r1 := kgen.GenSecretKey(user1)
		rlk1 := kgen.GenRelinearizationKey(sk1, r1)

		sk2, pk2 := kgen.GenKeyPair(user2)
		r2 := kgen.GenSecretKey(user2)
		rlk2 := kgen.GenRelinearizationKey(sk2, r2)

		skSet := NewSecretKeySet()
		skSet.AddSecretKey(sk1)
		skSet.AddSecretKey(sk2)

		rlkSet := NewRelinearizationKeyKeySet()
		rlkSet.AddRelinearizationKey(rlk1)
		rlkSet.AddRelinearizationKey(rlk2)

		enc := NewEncryptor(params)
		dec := NewDecryptor(params)

		ct1 := NewCiphertextNTT(params, idset1, level)
		ct2 := NewCiphertextNTT(params, idset2, level)

		enc.Encrypt(pt, pk1, ct1)
		enc.Encrypt(pt, pk2, ct2)

		ct3 := NewCiphertextNTT(params, idset, level)
		ks := NewKeySwitcher(params)
		ks.MulAndRelin(ct1, ct2, rlkSet, ct3)

		dec.Decrypt(ct1, skSet, pt)
		ringQ.InvNTTLvl(level, pt.Value, pt.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(level, ringQ, pt.Value))

		dec.Decrypt(ct2, skSet, pt)
		ringQ.InvNTTLvl(level, pt.Value, pt.Value)
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(level, ringQ, pt.Value))

		dec.Decrypt(ct3, skSet, pt)
		ringQ.InvNTTLvl(level, pt.Value, pt.Value)
		require.GreaterOrEqual(t, 2*(9+params.LogN()), log2OfInnerSum(level, ringQ, pt.Value))

	})

}

func testHadamardProduct(kgen *KeyGenerator, t *testing.T) {

	// Checks that internal product works properly
	// 1) generate two pk (-as+e, a)
	// 2) generate sg
	// 3) check Interal(-as+e, sg) similar to -as^2

	params := kgen.params

	t.Run(testString(params, "HadamardProductMaxLevel/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip()
		}

		var id string = "tetsUser1"
		users := NewIDSet()
		users.Add(id)

		ringQ := params.RingQ()
		levelQ := params.QCount() - 1

		level := levelQ

		sk := kgen.GenSecretKey(id)
		pk := kgen.GenPublicKey(sk)
		pt := rlwe.NewPlaintext(params.Parameters, level)
		pt.Value.IsNTT = true
		enc := NewEncryptor(params)
		ct := NewCiphertextNTT(params, users, level)
		enc.Encrypt(pt, pk, ct)

		//generate sg
		sg := NewSwitchingKey(params)
		kgen.GenSwitchingKey(sk, sg)

		ks := NewKeySwitcher(params)

		cd := NewSwitchingKey(params)
		ks.Decompose(level, ct.Value["0"], cd)

		beta := params.Beta(level)
		ringQP := params.RingQP()
		levelP := params.PCount() - 1

		for i := 0; i < beta; i++ {
			ringQP.MulCoeffsMontgomeryLvl(level, levelP, cd.Value[i], sg.Value[i], sg.Value[i])
			ringQP.MFormLvl(level, levelP, sg.Value[i], sg.Value[i])
		}

		//tmp = Inter(c, csg)
		tmp := ringQ.NewPolyLvl(level)
		tmp.IsNTT = true
		ks.InternalProduct(level, ct.Value["0"], sg, tmp)

		//tmp2 = c*s
		tmp2 := ringQ.NewPolyLvl(level)
		ringQ.MulCoeffsMontgomeryLvl(level, ct.Value["0"], sk.Value.Q, tmp2)
		ringQ.MFormLvl(level, tmp2, tmp2)
		ringQ.MulCoeffsMontgomeryAndSubLvl(level, ct.Value["0"], tmp2, tmp)
		ringQ.InvNTTLvl(level, tmp, tmp)

		//check errors
		require.GreaterOrEqual(t, 10+2*params.LogN(), log2OfInnerSum(tmp.Level(), params.RingQ(), tmp))
	})

}
