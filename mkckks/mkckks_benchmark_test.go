package mkckks

import (
	"strconv"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"mk-lattigo/mkrlwe"
)

func BenchmarkMKCKKS(b *testing.B) {

	defaultParams := []ckks.ParametersLiteral{PN14QP439, PN15QP880}

	for _, defaultParam := range defaultParams {
		ckksParams, err := ckks.NewParametersFromLiteral(defaultParam)

		if ckksParams.PCount() < 2 {
			continue
		}

		if err != nil {
			panic(err)
		}

		params := NewParameters(ckksParams)
		userList := make([]string, *maxUsers)
		idset := mkrlwe.NewIDSet()

		for i := range userList {
			userList[i] = "user" + strconv.Itoa(i)
			idset.Add(userList[i])
		}

		var testContext *testParams
		if testContext, err = genTestParams(params, idset); err != nil {
			panic(err)
		}

		for numUsers := 2; numUsers <= *maxUsers; numUsers *= 2 {
			benchMulAndRelin(testContext, userList[:numUsers], b)
			//benchMulAndRelinHoisted(testContext, userList[:numUsers], b)
			//benchSquareHoisted(testContext, userList[:numUsers], b)
		}

		/*

			for numUsers := 2; numUsers <= maxUsers; numUsers *= 2 {
				benchRotate(testContext, userList[:numUsers], b)
			}

		*/

	}
}

func benchMulAndRelin(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, 1), complex(-1, 1))
	}

	ct0 := ctList[0]
	ct1 := ctList[0]

	for i := range userList {
		ct0 = eval.AddNew(ct0, ctList[i])
		ct1 = eval.SubNew(ct1, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "MKMulAndRelin: "+strconv.Itoa(numUsers)+"/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MulRelinNew(ct0, ct1, rlkSet)
		}

	})
}

func benchRotate(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rtkSet := testContext.rtkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, 1), complex(-1, 1))
	}

	ct := ctList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "MKRotate: "+strconv.Itoa(numUsers)+"/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.RotateNew(ct, 2, rtkSet)
		}

	})
}

func benchMulAndRelinHoisted(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, 1), complex(-1, 1))
	}

	ct0 := ctList[0]
	ct1 := ctList[1]

	for i := range userList {
		ct0 = eval.AddNew(ct0, ctList[i])
		ct1 = eval.SubNew(ct1, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "MKMulAndRelinHoisted: "+strconv.Itoa(numUsers)+"/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ct0Hoisted := eval.HoistedForm(ct0)
			ct1Hoisted := eval.HoistedForm(ct1)
			eval.MulRelinHoistedNew(ct0, ct1, ct0Hoisted, ct1Hoisted, rlkSet)
		}

	})
}

func benchSquareHoisted(testContext *testParams, userList []string, b *testing.B) {

	numUsers := len(userList)
	msgList := make([]*Message, numUsers)
	ctList := make([]*Ciphertext, numUsers)

	rlkSet := testContext.rlkSet
	eval := testContext.evaluator

	for i := range userList {
		msgList[i], ctList[i] = newTestVectors(testContext, userList[i], complex(-1, 1), complex(-1, 1))
	}

	ct := ctList[0]

	for i := range userList {
		ct = eval.AddNew(ct, ctList[i])
	}

	b.Run(GetTestName(testContext.params, "SquareHoisted: "+strconv.Itoa(numUsers)+"/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ctHoisted := eval.HoistedForm(ct)
			eval.MulRelinHoistedNew(ct, ct, ctHoisted, ctHoisted, rlkSet)
		}

	})
}
