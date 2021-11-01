package mkckks

import (
	"flag"
	"fmt"
	"testing"

	"mk-lattigo/mkrlwe"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"

	"github.com/stretchr/testify/require"
	"math"
	/*

		"github.com/ldsec/lattigo/v2/rlwe"
			"encoding/json"
			"runtime"

			"github.com/stretchr/testify/assert"
	*/)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var minPrec float64 = 15.0

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.Beta())
}

type testParams struct {
	params    Parameters
	ringQ     *ring.Ring
	ringP     *ring.Ring
	prng      utils.PRNG
	kgen      *mkrlwe.KeyGenerator
	skSet     *mkrlwe.SecretKeySet
	pkSet     *mkrlwe.PublicKeySet
	rlkSet    *mkrlwe.RelinearizationKeySet
	encryptor *Encryptor
	decryptor *Decryptor
	evaluator *Evaluator
	idset     *mkrlwe.IDSet
}

func TestCKKS(t *testing.T) {

	defaultParams := ckks.DefaultParams[:4]                                      // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	defaultParams = append(ckks.DefaultParams, ckks.DefaultPostQuantumParams...) // the long test suite runs for all default parameters

	for _, defaultParam := range defaultParams {
		ckksParams, err := ckks.NewParametersFromLiteral(defaultParam)
		params := NewParameters(ckksParams)
		if err != nil {
			panic(err)
		}
		var testContext *testParams
		idset := mkrlwe.NewIDSet()
		idset.Add("user1")
		idset.Add("user2")

		if testContext, err = genTestParams(params, idset); err != nil {
			panic(err)
		}

		testEvaluatorAdd(testContext, t)
		testEvaluatorSub(testContext, t)
		testEvaluatorRescale(testContext, t)
		testEvaluatorMul(testContext, t)
	}
}

func genTestParams(defaultParam Parameters, idset *mkrlwe.IDSet) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam

	testContext.kgen = NewKeyGenerator(testContext.params)

	testContext.skSet = mkrlwe.NewSecretKeySet()
	testContext.pkSet = mkrlwe.NewPublicKeyKeySet()
	testContext.rlkSet = mkrlwe.NewRelinearizationKeyKeySet()

	for id := range idset.Value {
		sk, pk := testContext.kgen.GenKeyPair(id)
		r := testContext.kgen.GenSecretKey(id)
		rlk := testContext.kgen.GenRelinearizationKey(sk, r)

		testContext.skSet.AddSecretKey(sk)
		testContext.pkSet.AddPublicKey(pk)
		testContext.rlkSet.AddRelinearizationKey(rlk)

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

func newTestVectors(testContext *testParams, id string, a, b complex128, t *testing.T) (msg *Message, ciphertext *Ciphertext) {

	params := testContext.params
	logSlots := testContext.params.LogSlots()

	msg = NewMessage(params)

	for i := 0; i < 1<<logSlots; i++ {
		msg.Value[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	msg.Value[0] = complex(0.607538, 0)

	if testContext.encryptor != nil {
		ciphertext = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot newTestVectors: encryptor is not initialized!")
	}

	return msg, ciphertext
}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Add/CtCt/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1), t)
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1), t)
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
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1), t)
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1), t)
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

func testEvaluatorMul(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Mul/CtCt/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1), t)
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1), t)
		msg3 := NewMessage(params)

		for i := range msg3.Value {
			msg3.Value[i] = msg1.Value[i] * msg2.Value[i]
		}

		user1 := "user1"
		user2 := "user2"
		idset1 := mkrlwe.NewIDSet()
		idset2 := mkrlwe.NewIDSet()

		idset1.Add(user1)
		idset2.Add(user2)

		ct3 := testContext.evaluator.MulRelinNew(ct1, ct2, testContext.rlkSet)

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

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(GetTestName(testContext.params, "Evaluator/Rescale/Single/"), func(t *testing.T) {
		params := testContext.params
		msg1, ct1 := newTestVectors(testContext, "user1", complex(-1, -1), complex(1, 1), t)
		msg2, ct2 := newTestVectors(testContext, "user2", complex(-1, -1), complex(1, 1), t)
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
