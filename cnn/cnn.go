package cnn

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"math"
	"mk-lattigo/mkckks"
	"mk-lattigo/mkrlwe"
)

func Convolution(eval *mkckks.Evaluator, rlkSet *mkrlwe.RelinearizationKeySet, rtkSet *mkrlwe.RotationKeySet,
	ctImage *mkckks.Ciphertext, ctImageHoisted *mkrlwe.HoistedCiphertext,
	ctKernels []*mkckks.Ciphertext, ctKernelsHoisted []*mkrlwe.HoistedCiphertext) (convOut *mkckks.Ciphertext) {
	// OPTIMIZE: reuse encryptedImage three times - RotateNew & MulRelinNew
	// Note : encryptedKernels rotated by 0,1,14,15  and hoisting is applied

	convOut = eval.MulRelinHoistedNew(ctImage, ctKernels[0], ctImageHoisted, ctKernelsHoisted[0], rlkSet)

	temp := eval.RotateHoistedNew(ctImage, 1, ctImageHoisted, rtkSet)
	tempHoisted := eval.HoistedForm(temp)
	temp = eval.MulRelinHoistedNew(temp, ctKernels[1], tempHoisted, ctKernelsHoisted[1], rlkSet)
	convOut = eval.AddNew(convOut, temp)

	temp = eval.RotateHoistedNew(ctImage, 14, ctImageHoisted, rtkSet)
	tempHoisted = eval.HoistedForm(temp)
	temp = eval.MulRelinHoistedNew(temp, ctKernels[2], tempHoisted, ctKernelsHoisted[2], rlkSet)
	convOut = eval.AddNew(convOut, temp)

	temp = eval.RotateHoistedNew(ctImage, 15, ctImageHoisted, rtkSet)
	tempHoisted = eval.HoistedForm(temp)
	temp = eval.MulRelinHoistedNew(temp, ctKernels[3], tempHoisted, ctKernelsHoisted[3], rlkSet)
	convOut = eval.AddNew(convOut, temp)

	temp = eval.RotateNew(convOut, 2048, rtkSet)
	convOut = eval.AddNew(convOut, temp)

	temp = eval.RotateNew(convOut, 1024, rtkSet)
	convOut = eval.AddNew(convOut, temp)

	return convOut
}

func FC1Layer(eval *mkckks.Evaluator, rlkSet *mkrlwe.RelinearizationKeySet, rtkSet *mkrlwe.RotationKeySet,
	ctVec *mkckks.Ciphertext, ctVecHoisted *mkrlwe.HoistedCiphertext,
	ctMat []*mkckks.Ciphertext, ctMatHoisted []*mkrlwe.HoistedCiphertext,
	ctBias *mkckks.Ciphertext) (fc1Out *mkckks.Ciphertext) {

	var temp *mkckks.Ciphertext
	var tempHoisted *mkrlwe.HoistedCiphertext

	// OPTIMIZE: reuse ctVec - RotateNew & MulRelinNew
	for i := 0; i < len(ctMat); i++ {
		temp = ctVec.CopyNew()
		temp = eval.RotateHoistedNew(ctVec, i*128, ctVecHoisted, rtkSet)
		tempHoisted = eval.HoistedForm(temp)
		temp = eval.MulRelinHoistedNew(temp, ctMat[i], tempHoisted, ctMatHoisted[i], rlkSet)
		if i == 0 {
			fc1Out = temp.CopyNew()
		} else {
			fc1Out = eval.AddNew(fc1Out, temp)
		}
	}

	logn := int(math.Log2(float64(128)))
	for i := 0; i < logn; i++ {
		temp = eval.RotateNew(fc1Out, (1 << i), rtkSet)
		fc1Out = eval.AddNew(fc1Out, temp)
	}
	fc1Out = eval.AddNew(fc1Out, ctBias)

	return
}

func FC2Layer(eval *mkckks.Evaluator, rlkSet *mkrlwe.RelinearizationKeySet, rtkSet *mkrlwe.RotationKeySet,
	ctVec *mkckks.Ciphertext, ctMat *mkckks.Ciphertext,
	ctBias *mkckks.Ciphertext, ptMask *ckks.Plaintext) (fc2Out *mkckks.Ciphertext) {

	fc2Out = eval.MulPtxtNew(ctVec, ptMask)

	var temp *mkckks.Ciphertext
	logn := int(math.Log2(float64(16))) //10
	for i := 0; i < logn; i++ {
		temp = eval.RotateNew(fc2Out, -1*(1<<i), rtkSet)
		fc2Out = eval.AddNew(fc2Out, temp)
	}

	fc2Out = eval.MulRelinNew(fc2Out, ctMat, rlkSet)

	logn = int(math.Log2(float64(64)))
	for i := 0; i < logn; i++ {
		temp = eval.RotateNew(fc2Out, 128*(1<<i), rtkSet)
		fc2Out = eval.AddNew(fc2Out, temp)
	}

	fc2Out = eval.AddNew(fc2Out, ctBias)
	return
}
