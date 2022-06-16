package cnn

import (
	//"strconv"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"mk-lattigo/mkckks"
	"mk-lattigo/mkrlwe"
)

func BenchmarkCNN(b *testing.B) {
	defaultParam := PN14QP433
	ckksParams, err := ckks.NewParametersFromLiteral(defaultParam)
	params := mkckks.NewParameters(ckksParams)

	if err != nil {
		panic(err)
	}

	var testContext *testParams
	idset := mkrlwe.NewIDSet()
	idset.Add(modelOwner)
	idset.Add(dataOwner)

	if testContext, err = genTestParams(params, idset); err != nil {
		panic(err)
	}

	readData()

	// Precomputation
	imageIndex := 2
	ctImage := encryptImage(testContext, dataOwner, imageIndex)
	ctKernels := encryptKernels(testContext, modelOwner)
	ctFC1 := encryptFC1(testContext, modelOwner)
	ctFC2 := encryptFC2(testContext, modelOwner)
	ctB1 := encryptB1(testContext, modelOwner)
	ctB2 := encryptB2(testContext, modelOwner)

	eval := testContext.evaluator
	rlkSet := testContext.rlkSet

	ctImageHoisted := eval.HoistedForm(ctImage)
	ctKernelsHoisted := make([]*mkrlwe.HoistedCiphertext, len(ctKernels))
	for i := 0; i < len(ctKernelsHoisted); i++ {
		ctKernelsHoisted[i] = eval.HoistedForm(ctKernels[i])
	}
	ctFC1Hoisted := make([]*mkrlwe.HoistedCiphertext, len(ctFC1))
	for i := 0; i < len(ctFC1Hoisted); i++ {
		ctFC1Hoisted[i] = eval.HoistedForm(ctFC1[i])
	}

	numSlots := testContext.params.Slots()
	mask := make([]complex128, numSlots)
	for i := 0; i < numSlots; i += 128 {
		mask[i] = 1
	}
	msg := mkckks.NewMessage(testContext.params)
	msg.Value = mask
	ptMask := testContext.encryptor.EncodeMsgNew(msg)

	// Evaluation
	convOut := benchConvolution(testContext, b, ctImage, ctImageHoisted, ctKernels, ctKernelsHoisted)
	convOutHoisted := eval.HoistedForm(convOut)
	square1Out := eval.MulRelinHoistedNew(convOut, convOut, convOutHoisted, convOutHoisted, rlkSet)
	square1OutHoisted := eval.HoistedForm(square1Out)

	fc1Out := benchFC1(testContext, b, square1Out, square1OutHoisted, ctFC1, ctFC1Hoisted, ctB1)
	fc1OutHoisted := eval.HoistedForm(fc1Out)
	square2Out := eval.MulRelinHoistedNew(fc1Out, fc1Out, fc1OutHoisted, fc1OutHoisted, rlkSet)

	benchFC2(testContext, b, square2Out, ctFC2, ctB2, ptMask)

}

func benchConvolution(testContext *testParams, b *testing.B,
	ctImage *mkckks.Ciphertext, ctImageHoisted *mkrlwe.HoistedCiphertext,
	ctKernels []*mkckks.Ciphertext, ctKernelsHoisted []*mkrlwe.HoistedCiphertext) (convOut *mkckks.Ciphertext) {

	b.Run(GetTestName(testContext.params, "Convolution/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval := testContext.evaluator
			rlkSet := testContext.rlkSet
			rtkSet := testContext.rtkSet
			convOut = Convolution(eval, rlkSet, rtkSet, ctImage, ctImageHoisted, ctKernels, ctKernelsHoisted)
		}

	})

	return
}

func benchFC1(testContext *testParams, b *testing.B,
	ctVec *mkckks.Ciphertext, ctVecHoisted *mkrlwe.HoistedCiphertext,
	ctMat []*mkckks.Ciphertext, ctMatHoisted []*mkrlwe.HoistedCiphertext, ctBias *mkckks.Ciphertext) (fc1Out *mkckks.Ciphertext) {

	b.Run(GetTestName(testContext.params, "FC1/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval := testContext.evaluator
			rlkSet := testContext.rlkSet
			rtkSet := testContext.rtkSet
			fc1Out = FC1Layer(eval, rlkSet, rtkSet, ctVec, ctVecHoisted, ctMat, ctMatHoisted, ctBias)
		}

	})

	return
}

func benchFC2(testContext *testParams, b *testing.B,
	ctVec *mkckks.Ciphertext, ctMat *mkckks.Ciphertext, ctBias *mkckks.Ciphertext, ptMask *ckks.Plaintext) (fc2Out *mkckks.Ciphertext) {

	b.Run(GetTestName(testContext.params, "FC2/ "), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval := testContext.evaluator
			rlkSet := testContext.rlkSet
			rtkSet := testContext.rtkSet
			fc2Out = FC2Layer(eval, rlkSet, rtkSet, ctVec, ctMat, ctBias, ptMask)
		}

	})

	return
}
