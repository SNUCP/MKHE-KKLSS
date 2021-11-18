package mkbfv

import (
	"strconv"
	"testing"

	"mk-lattigo/mkrlwe"
)

func BenchmarkMKCKKS(b *testing.B) {

	defaultParams := []ParametersLiteral{PN15QP873}

	for _, defaultParam := range defaultParams {
		params := NewParametersFromLiteral(defaultParam)

		maxUsers := 64
		userList := make([]string, maxUsers)
		idset := mkrlwe.NewIDSet()

		for i := range userList {
			userList[i] = "user" + strconv.Itoa(i)
			idset.Add(userList[i])
		}

		testContext, _ := genTestParams(params, idset)

		benchMulAndRelin(testContext, userList[:2], b)
		benchMulAndRelin(testContext, userList[:4], b)
		benchMulAndRelin(testContext, userList[:8], b)
		benchMulAndRelin(testContext, userList[:16], b)
		benchMulAndRelin(testContext, userList[:32], b)
		benchMulAndRelin(testContext, userList[:64], b)
	}
}

func benchMulAndRelin(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], 0, 2)
	}

	ct := ctList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "MKMulAndRelin: "+strconv.Itoa(numUsers)+"/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct, ct, rlkSet)
		}

	})
}
