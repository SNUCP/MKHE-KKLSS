package mkbfv

import "mk-lattigo/mkrlwe"

import "github.com/stretchr/testify/require"
import "github.com/ldsec/lattigo/v2/rlwe"
import "github.com/ldsec/lattigo/v2/ring"
import "github.com/ldsec/lattigo/v2/utils"

import "fmt"
import "testing"
import "strconv"

//import "math"

func GetTestName(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQP=%d/",
		opname,
		params.paramsQP.LogN(),
		params.paramsQP.LogQP(),
	)
}

var PN15QP840 = ParametersLiteral{
	LogN: 15,
	Q: []uint64{
		0x200000440001, 0x200000500001, 0x200000620001,
		0x1fffff980001, 0x2000006a0001, 0x1fffff7e0001,
		0x200000860001, 0x200000a60001, 0x200000aa0001,
		0x200000b20001, 0x200000c80001, 0x1fffff360001,
		0x200000e20001, 0x1fffff060001, 0x200000fe0001,
	}, // 15x45
	QMul: []uint64{
		0x1ffffede0001, 0x1ffffeca0001, 0x1ffffeb40001,
		0x200001520001, 0x1ffffe760001, 0x2000019a0001,
		0x1ffffe640001, 0x200001a00001, 0x1ffffe520001,
		0x200001e80001, 0x1ffffe0c0001, 0x1ffffdee0001,
		0x200002480001, 0x2000000a0001, 0x2000000e0001,
	}, // 15x45
	P: []uint64{
		0x80000000440001, 0x7fffffffba0001, 0x80000000500001,
	}, // 3x55
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
	encryptor *Encryptor
	decryptor *Decryptor
	//evaluator *Evaluator
	idset *mkrlwe.IDSet
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

	return testContext, nil

}

func TestCKKS(t *testing.T) {
	defaultParams := []ParametersLiteral{PN15QP840}
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
