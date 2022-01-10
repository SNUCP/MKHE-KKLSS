package mkckks

import (
	"flag"
	"fmt"
	"strconv"
	"testing"

	"mk-lattigo/mkrlwe"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"

	"github.com/stretchr/testify/require"
	"math"
	"math/cmplx"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var minPrec float64 = 15.0

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1)
}

type testParams struct {
	params Parameters
	ringQ  *ring.Ring
	ringP  *ring.Ring
	prng   utils.PRNG
	kgen   *mkrlwe.KeyGenerator
	skSet  *mkrlwe.SecretKeySet
	pkSet  *mkrlwe.PublicKeySet
	rlkSet *mkrlwe.RelinearizationKeySet
	rtkSet *mkrlwe.RotationKeySet
	cjkSet *mkrlwe.ConjugationKeySet

	encryptor *Encryptor
	decryptor *Decryptor
	evaluator *Evaluator
	idset     *mkrlwe.IDSet
}

var (
	PN15QP878 = ckks.ParametersLiteral{
		LogN:     15,
		LogSlots: 14,
		//56 + 13*46
		Q: []uint64{
			0x100000000060001,

			0x400000060001,
			0x400000080001,
			0x4000000c0001,
			0x400000180001,
			0x400000290001,
			0x400000420001,
			0x3ffffff70001,
			0x3fffffe50001,
			0x3fffffcc0001,
			0x3fffffc70001,
			0x3fffffb20001,
			0x3fffffa80001,
			0x3fffff960001,
		},
		P: []uint64{
			//56 * 4
			0x1000000002a0001,
			0x100000000450001,
			0x100000000480001,
			0x1000000005f0001,
		},
		Scale: 1 << 46,
		Sigma: rlwe.DefaultSigma,
	}
	// PN14QP438 is a default parameter set for logN=14 and logQP=438
	PN14QP441 = ckks.ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Q: []uint64{
			// 49 + 5*39
			0x20000000b0001,

			0x80001d0001,
			0x8000410001,
			0x8000430001,
			0x7ffffb0001,
			0x7fffe60001,
		},
		P: []uint64{
			//49*4
			0x20000001a0001,
			0x20000003b0001,
			0x20000005e0001,
			0x20000006d0001,
		},
		Scale: 1 << 39,
		Sigma: rlwe.DefaultSigma,
	}
)

func TestCKKS(t *testing.T) {

	defaultParams := []ckks.ParametersLiteral{PN15QP878, PN14QP441}

	for _, defaultParam := range defaultParams {
		ckksParams, err := ckks.NewParametersFromLiteral(defaultParam)

		if ckksParams.PCount() < 2 {
			continue
		}

		if err != nil {
			panic(err)
		}

		params := NewParameters(ckksParams)
		maxUsers := 4
		userList := make([]string, maxUsers)
		idset := mkrlwe.NewIDSet()

		for i := range userList {
			userList[i] = "user" + strconv.Itoa(i)
			idset.Add(userList[i])
		}

		var testContext *testParams
		if testContext, err = genTestParams(params, idset); err != nil {
			panic(err)
		}

		testEncAndDec(testContext, userList, t)

		for numUsers := 2; numUsers <= maxUsers; numUsers *= 2 {
			testEvaluatorMul(testContext, userList[:numUsers], t)
			testEvaluatorRot(testContext, userList[:numUsers], t)
			testEvaluatorConj(testContext, userList[:numUsers], t)
		}
	}
}

func genTestParams(defaultParam Parameters, idset *mkrlwe.IDSet) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam

	testContext.kgen = NewKeyGenerator(testContext.params)

	testContext.skSet = mkrlwe.NewSecretKeySet()
	testContext.pkSet = mkrlwe.NewPublicKeyKeySet()
	testContext.rlkSet = mkrlwe.NewRelinearizationKeyKeySet()
	testContext.rtkSet = mkrlwe.NewRotationKeySet()
	testContext.cjkSet = mkrlwe.NewConjugationKeySet()

	// gen sk, pk, rlk, rk

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

	testContext.evaluator = NewEvaluator(testContext.params)

	return testContext, nil

}

func newTestVectors(testContext *testParams, id string, a, b complex128) (msg *Message, ciphertext *Ciphertext) {

	params := testContext.params
	logSlots := testContext.params.LogSlots()

	msg = NewMessage(params)

	for i := 0; i < 1<<logSlots; i++ {
		msg.Value[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
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
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, -1), complex(1, 1))
	}

	t.Run(GetTestName(testContext.params, "MKCKKSEncAndDec: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {

		for i := range userList {
			msgOut := dec.Decrypt(ctList[i], skSet)
			for j := range msgList[i].Value {
				delta := msgList[i].Value[j] - msgOut.Value[j]
				require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
			}
		}
	})

}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Add/CtCt/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1))
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1))
		msg3 := NewMessage(params)

		for i := range msg3.Value {
			msg3.Value[i] = msg1.Value[i] + msg2.Value[i]
		}

		user1 := "user1"
		user2 := "user2"
		idset1 := mkrlwe.NewIDSet()
		idset2 := mkrlwe.NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)

		ct3 := testContext.evaluator.AddNew(ct1, ct2)

		msg1Out := testContext.decryptor.Decrypt(ct1, testContext.skSet)
		msg2Out := testContext.decryptor.Decrypt(ct2, testContext.skSet)
		msg3Out := testContext.decryptor.Decrypt(ct3, testContext.skSet)

		for i := range msg1Out.Value {
			delta := msg1.Value[i] - msg1Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg2Out.Value {
			delta := msg2.Value[i] - msg2Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg3Out.Value {
			delta := msg3.Value[i] - msg3Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

	})

}

func testEvaluatorSub(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Sub/CtCt/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1))
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1))
		msg3 := NewMessage(params)

		for i := range msg3.Value {
			msg3.Value[i] = msg1.Value[i] - msg2.Value[i]
		}

		user1 := "user1"
		user2 := "user2"
		idset1 := mkrlwe.NewIDSet()
		idset2 := mkrlwe.NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)

		ct3 := testContext.evaluator.SubNew(ct1, ct2)

		msg1Out := testContext.decryptor.Decrypt(ct1, testContext.skSet)
		msg2Out := testContext.decryptor.Decrypt(ct2, testContext.skSet)
		msg3Out := testContext.decryptor.Decrypt(ct3, testContext.skSet)

		for i := range msg1Out.Value {
			delta := msg1.Value[i] - msg1Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg2Out.Value {
			delta := msg2.Value[i] - msg2Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg3Out.Value {
			delta := msg3.Value[i] - msg3Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

	})

}

func testEvaluatorMul(testContext *testParams, userList []string, t *testing.T) {

	params := testContext.params
	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i],
			complex(0.1/float64(numUsers), 1.0/float64(numUsers)),
			complex(0.1/float64(numUsers), 1.0/float64(numUsers)))
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

	t.Run(GetTestName(testContext.params, "MKMulAndRelin: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := eval.MulRelinNew(ct, ct, rlkSet)
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := range msgRes.Value {
			delta := msgRes.Value[i] - msg.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(real(delta))))
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(imag(delta))))
		}
	})

}

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Rescale/Single/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1))
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1))
		msg3 := NewMessage(params)

		constant := testContext.ringQ.Modulus[ct1.Level()]

		for i := range msg3.Value {
			msg3.Value[i] = msg1.Value[i] + msg2.Value[i]
		}

		user1 := "user1"
		user2 := "user2"
		idset1 := mkrlwe.NewIDSet()
		idset2 := mkrlwe.NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)

		ct3 := testContext.evaluator.AddNew(ct1, ct2)

		testContext.evaluator.MultByConst(ct3, constant, ct3)
		ct3.Scale *= float64(constant)
		testContext.evaluator.Rescale(ct3, params.Scale(), ct3)

		msg1Out := testContext.decryptor.Decrypt(ct1, testContext.skSet)
		msg2Out := testContext.decryptor.Decrypt(ct2, testContext.skSet)
		msg3Out := testContext.decryptor.Decrypt(ct3, testContext.skSet)

		for i := range msg1Out.Value {
			delta := msg1.Value[i] - msg1Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg2Out.Value {
			delta := msg2.Value[i] - msg2Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
		}

		for i := range msg3Out.Value {
			delta := msg3.Value[i] - msg3Out.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+8, math.Log2(math.Abs(real(delta))))
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

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i],
			complex(0.5/float64(numUsers), 1.0/float64(numUsers)),
			complex(0.5/float64(numUsers), 1.0/float64(numUsers)))
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

		for i := range msgRes.Value {
			var delta complex128
			if rot > 0 {
				delta = msgRes.Value[i] - msg.Value[(i+rot)%len(msg.Value)]
			} else {
				delta = msg.Value[i] - msgRes.Value[(i-rot)%len(msg.Value)]
			}
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(real(delta))))
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(imag(delta))))
		}
	})

}

func testEvaluatorConj(testContext *testParams, userList []string, t *testing.T) {

	params := testContext.params
	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	cjkSet := testContext.cjkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i],
			complex(0.5/float64(numUsers), 1.0/float64(numUsers)),
			complex(0.5/float64(numUsers), 1.0/float64(numUsers)))
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
		msg.Value[j] = cmplx.Conj(msg.Value[j])
	}

	t.Run(GetTestName(testContext.params, "MKConjugate: "+strconv.Itoa(numUsers)+"/ "), func(t *testing.T) {
		ctRes := eval.ConjugateNew(ct, cjkSet)
		msgRes := testContext.decryptor.Decrypt(ctRes, testContext.skSet)

		for i := range msgRes.Value {
			delta := msgRes.Value[i] - msg.Value[i]
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(real(delta))))
			require.GreaterOrEqual(t, -math.Log2(params.Scale())+float64(params.LogSlots())+11, math.Log2(math.Abs(imag(delta))))
		}
	})

}
