package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/stretchr/testify/require"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"

import "fmt"
import "testing"
import "strconv"
import "math/big"
import "math/bits"

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQP=%d/logRP=%d/",
		opname,
		params.LogN(),
		params.LogQP(),
		params.paramsRP.LogQP(),
	)
}

var PN15QP873 = ParametersLiteral{
	LogN: 15,

	Q: []uint64{
		// 12 x 45
		0x200000b20001, 0x200000c80001, 0x1fffff360001,
		0x200000e20001, 0x1fffff060001, 0x200000fe0001,
		0x1ffffede0001, 0x1ffffeca0001, 0x1ffffeb40001,
		0x200001520001, 0x1ffffe760001, 0x2000019a0001,
	},

	QMul: []uint64{
		// 12 x 45
		0x2000000a0001, 0x2000000e0001, 0x1fffffc20001,
		0x200000440001, 0x200000500001, 0x200000620001,
		0x1fffff980001, 0x2000006a0001, 0x1fffff7e0001,
		0x200000860001, 0x200000a60001, 0x200000aa0001,
	},

	P: []uint64{
		// 6 x 55
		0x80000000080001, 0x80000000130001, 0x80000000190001,
		0x800000001d0001, 0x80000000440001, 0x80000000490001,
	},
	T:     65537,
	Sigma: rlwe.DefaultSigma,
}

var PN14QP435 = ParametersLiteral{
	LogN: 14,

	Q: []uint64{
		0x10000140001,
		0x100003e0001,
		0x100004b0001,
		0x10000650001,
		0x10000960001,
		0x10000ab0001,
		// 6 x 40
	},

	QMul: []uint64{
		// 6 x 40
		0x10000290001,
		0x10000470001,
		0x10000500001,
		0x10000890001,
		0x10000a40001,
		0x10000b60001,
	},

	P: []uint64{
		// 50 x 4
		0x4000000120001,
		0x40000001b0001,
		0x4000000270001,
		0x4000000350001,
	},
	T:     65537,
	Sigma: rlwe.DefaultSigma,
}

var PN13QP220 = ParametersLiteral{
	LogN: 13,

	Q: []uint64{
		// 4 x 35
		0x800008001,
		0x800250001,
		0x800280001,
		0x8002f0001,
	},

	QMul: []uint64{
		// 4 x 35
		0x8000f8001,
		0x800260001,
		0x8002a8001,
		0x800398001,
	},

	P: []uint64{
		// 40 x 2
		0x10000290001,
		0x10000470001,
	},
	T:     65537,
	Sigma: rlwe.DefaultSigma,
}

type testParams struct {
	params    Parameters
	ringQ     *ring.Ring
	ringP     *ring.Ring
	prng      utils.PRNG
	kgen      *KeyGenerator
	skSet     *mkrlwe.SecretKeySet
	pkSet     *mkrlwe.PublicKeySet
	rlkSet    *mkrlwe.RelinearizationKeySet
	rtkSet    *mkrlwe.RotationKeySet
	cjkSet    *mkrlwe.ConjugationKeySet
	encryptor *Encryptor
	decryptor *Decryptor
	evaluator *Evaluator
	idset     *mkrlwe.IDSet
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

func genTestParams(defaultParam Parameters, idset *mkrlwe.IDSet) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam

	testContext.kgen = NewKeyGenerator(testContext.params)
	testContext.evaluator = NewEvaluator(testContext.params)

	testContext.skSet = mkrlwe.NewSecretKeySet()
	testContext.pkSet = mkrlwe.NewPublicKeyKeySet()
	testContext.rlkSet = mkrlwe.NewRelinearizationKeyKeySet()
	testContext.rtkSet = mkrlwe.NewRotationKeySet()
	testContext.cjkSet = mkrlwe.NewConjugationKeySet()

	testContext.rlkSet.AddRelinearizationKey(testContext.kgen.GenConstRelinearizationKey())

	for id := range idset.Value {
		sk, pk := testContext.kgen.GenKeyPair(id)
		r := testContext.kgen.GenSecretKey(id)
		rlk := testContext.kgen.GenRelinearizationKey(sk, r)
		cjk := testContext.kgen.GenConjugationKey(sk)

		testContext.kgen.GenDefaultRotationKeys(sk, testContext.rtkSet)

		testContext.skSet.AddSecretKey(sk)
		testContext.pkSet.AddPublicKey(pk)
		testContext.rlkSet.AddRelinearizationKey(rlk)
		testContext.cjkSet.AddConjugationKey(cjk)
	}

	testContext.ringQ = defaultParam.RingQ()

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encryptor = NewEncryptor(testContext.params)
	testContext.decryptor = NewDecryptor(testContext.params)

	return testContext, nil

}

func TestMKBFV(t *testing.T) {
	defaultParams := []ParametersLiteral{PN13QP220, PN14QP435, PN15QP873}
	for _, defaultParam := range defaultParams {
		params := NewParametersFromLiteral(defaultParam)

		maxUsers := 4
		userList := make([]string, maxUsers)
		idset := mkrlwe.NewIDSet()

		for i := range userList {
			userList[i] = "user" + strconv.Itoa(i)
			idset.Add(userList[i])
		}

		testContext, _ := genTestParams(params, idset)

		testEncAndDec(testContext, userList, t)

		for numUsers := 2; numUsers <= maxUsers; numUsers *= 2 {
			testEvaluatorAdd(testContext, userList[:numUsers], t)
			testEvaluatorSub(testContext, userList[:numUsers], t)
			testEvaluatorMul(testContext, userList[:numUsers], t)
			testEvaluatorRot(testContext, userList[:numUsers], t)
			testEvaluatorConj(testContext, userList[:numUsers], t)
		}
	}
}

func newTestVectors(testContext *testParams, id string, a, b int64) (msg *Message, ciphertext *Ciphertext) {

	params := testContext.params
	msg = NewMessage(params)

	for i := 0; i < params.N(); i++ {
		msg.Value[i] = int64(utils.RandUint64()/2)%(b-a) + a
	}

	if testContext.encryptor != nil {
		ciphertext = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot newTestVectors: encryptor is not initialized!")
	}

	return msg, ciphertext
}

func testEncAndDec(testContext *testParams, userList []string, t *testing.T) {
	params := testContext.params
	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	skSet := testContext.skSet
	dec := testContext.decryptor

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], -int64(params.T())/4, int64(params.T())/4)
	}

	t.Run(GetTestName(testContext.params, "MKBFVEncAndDec: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {

		for i := range userList {
			msgOut := dec.Decrypt(ctList[i], skSet)
			for j := range msgList[i].Value {
				delta := msgList[i].Value[j] - msgOut.Value[j]
				require.Equal(t, int64(0), delta, fmt.Sprintf("%v vs %v", msgList[i].Value[j], msgOut.Value[j]))
			}
		}
	})

}

func testEvaluatorAdd(testContext *testParams, userList []string, t *testing.T) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], -100, -20)
	}

	ct := ctList[0]
	msg := msgList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])

		for j := range msg.Value {
			msg.Value[j] += msgList[i].Value[j]
		}
	}

	t.Run(GetTestName(testContext.params, "MKBFVAdd: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := ct
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := range msgRes.Value {
			delta := msgRes.Value[i] - msg.Value[i]
			require.Equal(t, int64(0), delta, fmt.Sprintf("%v vs %v", msgRes.Value[i], msg.Value[i]))
		}
	})

}

func testEvaluatorSub(testContext *testParams, userList []string, t *testing.T) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], -2, 2)
	}

	ct := ctList[0]
	msg := msgList[0]

	for i := range userList {
		ct = eval.SubNew(ct, ctList[i])

		for j := range msg.Value {
			msg.Value[j] -= msgList[i].Value[j]
		}
	}

	t.Run(GetTestName(testContext.params, "MKBFVSub: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := ct
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := range msgRes.Value {
			delta := msgRes.Value[i] - msg.Value[i]
			require.Equal(t, int64(0), delta, fmt.Sprintf("%v vs %v", msgRes.Value[i], msg.Value[i]))
		}
	})

}

func testEvaluatorMul(testContext *testParams, userList []string, t *testing.T) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], 0, 2)
	}

	ct := ctList[0]
	msg := msgList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])

		for j := range msg.Value {
			msg.Value[j] += msgList[i].Value[j]
		}
	}

	for j := range msg.Value {
		msg.Value[j] *= msg.Value[j]
	}

	t.Run(GetTestName(testContext.params, "MKBFVMulAndRelin: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := eval.MulRelinNew(ct, ct, rlkSet)
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := range msgRes.Value {
			delta := msgRes.Value[i] - msg.Value[i]
			require.Equal(t, int64(0), delta, fmt.Sprintf("%v: %v vs %v", i, msgRes.Value[i], msg.Value[i]))
		}
	})

}

func testEvaluatorRot(testContext *testParams, userList []string, t *testing.T) {

	params := testContext.params
	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rtkSet := testContext.rtkSet
	eval := testContext.evaluator
	slots := eval.params.N() / 2

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], 0, 2)
	}

	ct := ctList[0]
	msg := msgList[0]
	rot := int(utils.RandUint64()%uint64(2*params.N())) - params.N()

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])

		for j := range msg.Value {
			msg.Value[j] += msgList[i].Value[j]
		}
	}

	t.Run(GetTestName(testContext.params, "MKRotate: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := eval.RotateNew(ct, rot, rtkSet)
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := 0; i < slots; i++ {
			var delta int64
			if rot > 0 {
				delta = msgRes.Value[i] - msg.Value[(i+rot)%slots]
			} else {
				delta = msg.Value[i] - msgRes.Value[(i-rot)%slots]
			}
			require.Equal(t, int64(0), delta)
		}

		for i := 0; i < slots; i++ {
			var delta int64
			if rot > 0 {
				delta = msgRes.Value[i+slots] - msg.Value[(i+rot)%slots+slots]
			} else {
				delta = msg.Value[i+slots] - msgRes.Value[(i-rot)%slots+slots]
			}
			require.Equal(t, int64(0), delta)
		}

	})

}

func testEvaluatorConj(testContext *testParams, userList []string, t *testing.T) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	cjkSet := testContext.cjkSet
	eval := testContext.evaluator
	slots := eval.params.N() / 2

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], 0, 2)
	}

	ct := ctList[0]
	msg := msgList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])

		for j := range msg.Value {
			msg.Value[j] += msgList[i].Value[j]
		}
	}

	t.Run(GetTestName(testContext.params, "MKConjugate: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := eval.ConjugateNew(ct, cjkSet)
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := 0; i < slots; i++ {
			delta := msgRes.Value[i] - msg.Value[(i+slots)]
			require.Equal(t, int64(0), delta)
		}

		for i := 0; i < slots; i++ {
			delta := msgRes.Value[i+slots] - msg.Value[i]
			require.Equal(t, int64(0), delta)
		}

	})

}
