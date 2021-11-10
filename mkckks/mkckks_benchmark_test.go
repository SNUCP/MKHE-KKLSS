package mkckks

import (
	"strconv"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"mk-lattigo/mkrlwe"
)

func BenchmarkMKCKKS(b *testing.B) {

	defaultParams := []ckks.ParametersLiteral{PN15QP830}

	for _, defaultParam := range defaultParams {
		ckksParams, err := ckks.NewParametersFromLiteral(defaultParam)

		if ckksParams.PCount() < 2 {
			continue
		}

		if err != nil {
			panic(err)
		}

		params := NewParameters(ckksParams)
		maxUsers := 16
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

		benchRelinAndMul(testContext, userList[:2], b)
		benchRelinAndMul(testContext, userList[:4], b)
		benchRelinAndMul(testContext, userList[:8], b)
		benchRelinAndMul(testContext, userList[:16], b)
	}
}

func benchRelinAndMul(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, 1), complex(1, 1))
	}

	ctRes := ctList[0]

	for i := range userList {
		ctRes = eval.AddNew(ctRes, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "Evaluator/Mul/CtCt/"+strconv.Itoa(numUsers)), func(b *testing.B) {
		ctRes = eval.MulRelinNew(ctRes, ctRes, rlkSet)
	})
}
