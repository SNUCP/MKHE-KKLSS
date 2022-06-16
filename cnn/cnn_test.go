package cnn

import (
	"fmt"
	"io/ioutil"
	"strconv"
	"strings"

	"github.com/stretchr/testify/require"
	"testing"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"

	"mk-lattigo/mkckks"
	"mk-lattigo/mkrlwe"
)

var (
	imageSize = 28
	numClass  = 10
	numImage  = 10000

	numKernels  = 5
	kernelSize  = 4
	stride      = 2
	blockSize   = 14
	convOutSize = 13
	numFCUnit   = 64

	gap = 128

	dataOwner  = "dataOwner"
	modelOwner = "modelOwner"
)

func GetTestName(params mkckks.Parameters, opname string) string {
	return fmt.Sprintf("%sLogN=%d/LogSlots=%d/LogQP=%d/Levels=%d/",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1)
}

type testParams struct {
	params mkckks.Parameters
	ringQ  *ring.Ring
	ringP  *ring.Ring
	prng   utils.PRNG
	kgen   *mkrlwe.KeyGenerator
	skSet  *mkrlwe.SecretKeySet
	pkSet  *mkrlwe.PublicKeySet
	rlkSet *mkrlwe.RelinearizationKeySet
	rtkSet *mkrlwe.RotationKeySet

	encryptor *mkckks.Encryptor
	decryptor *mkckks.Decryptor
	evaluator *mkckks.Evaluator
	idset     *mkrlwe.IDSet
}

const imagefile = "./data/mnist_test.csv"
const k0file = "./data/k1.txt"
const FC1file = "./data/FC1.txt"
const FC2file = "./data/FC2.txt"
const B1file = "./data/B1.txt"
const B2file = "./data/B2.txt"

var (
	classes   []complex128     // 10000
	images    [][][]complex128 // 10000 x 28 x 28
	kernels   [][][]complex128 // 5 x 4 x 4
	FC1       [][]complex128   // 845 x 64
	FC2       [][]complex128   // 64 x 10
	B1        []complex128     // 64
	B2        []complex128     // 10
	PN14QP433 = ckks.ParametersLiteral{
		LogN:     14,
		LogSlots: 13,
		Q: []uint64{
			//57 + 47 x 6
			0x2000000002b0001,
			0x800000020001, 0x800000280001,
			0x800000520001, 0x800000770001,
			0x800000aa0001, 0x800000ad0001,
		},
		P: []uint64{
			//47 x 2
			0x800000df0001, 0x800000f80001,
		},
		Scale: 1 << 47,
		Sigma: rlwe.DefaultSigma,
	}
)

func TestCNN(t *testing.T) {
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
	rtkSet := testContext.rtkSet

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
	convOut := Convolution(eval, rlkSet, rtkSet, ctImage, ctImageHoisted, ctKernels, ctKernelsHoisted)

	convOutHoisted := eval.HoistedForm(convOut)
	square1Out := eval.MulRelinHoistedNew(convOut, convOut, convOutHoisted, convOutHoisted, rlkSet)
	square1OutHoisted := eval.HoistedForm(square1Out)

	fc1Out := FC1Layer(eval, rlkSet, rtkSet, square1Out, square1OutHoisted, ctFC1, ctFC1Hoisted, ctB1)
	fc1OutHoisted := eval.HoistedForm(fc1Out)
	square2Out := eval.MulRelinHoistedNew(fc1Out, fc1Out, fc1OutHoisted, fc1OutHoisted, rlkSet)

	fc2Out := FC2Layer(eval, rlkSet, rtkSet, square2Out, ctFC2, ctB2, ptMask)

	// Decrypt
	ptResult := testContext.decryptor.Decrypt(fc2Out, testContext.skSet)
	value := ptResult.Value
	maxIndex := -1
	maxValue := -100.0
	for i := 0; i < 10; i++ {
		if real(value[i]) > maxValue {
			maxValue = real(value[i])
			maxIndex = i
		}
	}

	answer := int(real(classes[imageIndex]))
	require.Equal(t, answer, maxIndex)
}

func genTestParams(defaultParam mkckks.Parameters, idset *mkrlwe.IDSet) (testContext *testParams, err error) {
	testContext = new(testParams)

	testContext.params = defaultParam

	rots := []int{14, 15, 384, 512, 640, 768, 896, 8191, 8190, 8188, 8184}

	for _, rot := range rots {
		testContext.params.AddCRS(rot)
	}

	testContext.kgen = mkckks.NewKeyGenerator(testContext.params)

	testContext.skSet = mkrlwe.NewSecretKeySet()
	testContext.pkSet = mkrlwe.NewPublicKeyKeySet()
	testContext.rlkSet = mkrlwe.NewRelinearizationKeyKeySet()
	testContext.rtkSet = mkrlwe.NewRotationKeySet()

	for i := 0; i < testContext.params.LogN()-1; i++ {
		rots = append(rots, 1<<i)
	}

	for id := range idset.Value {
		sk, pk := testContext.kgen.GenKeyPair(id)
		r := testContext.kgen.GenSecretKey(id)
		rlk := testContext.kgen.GenRelinearizationKey(sk, r)

		for _, rot := range rots {
			rk := testContext.kgen.GenRotationKey(rot, sk)
			testContext.rtkSet.AddRotationKey(rk)
		}

		testContext.skSet.AddSecretKey(sk)
		testContext.pkSet.AddPublicKey(pk)
		testContext.rlkSet.AddRelinearizationKey(rlk)
	}

	testContext.ringQ = defaultParam.RingQ()

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encryptor = mkckks.NewEncryptor(testContext.params)
	testContext.decryptor = mkckks.NewDecryptor(testContext.params)
	testContext.evaluator = mkckks.NewEvaluator(testContext.params)

	return testContext, nil
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data Reading Functions //////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

func readData() {
	readTestData(imagefile)
	readKData(k0file)
	B1 = readBData(B1file, numFCUnit)
	B2 = readBData(B2file, numClass)
	FC1 = readFCData(FC1file, convOutSize*convOutSize*numKernels, numFCUnit)
	FC2 = readFCData(FC2file, numFCUnit, numClass)
}

func readTestData(filename string) {
	classes = make([]complex128, numImage)
	images = make([][][]complex128, numImage)
	for i := 0; i < numImage; i++ {
		images[i] = make([][]complex128, imageSize)
		for j := 0; j < imageSize; j++ {
			images[i][j] = make([]complex128, imageSize)
		}
	}

	f, error := ioutil.ReadFile(filename)
	if error != nil {
		panic(error)
	}

	lines := strings.Split(string(f), "\n")
	for i, l := range lines {
		values := strings.Split(l, ",")
		for j, v := range values {
			value, error := strconv.ParseFloat(v, 64)

			if error != nil {
				panic(error)
			}
			if j == 0 {
				classes[i] = complex(value, 0)
			} else {
				images[i][(j-1)/imageSize][(j-1)%imageSize] = complex(value/255, 0)
			}
		}
	}
}

func readKData(filename string) {
	kernels = make([][][]complex128, numKernels)
	for i := 0; i < numKernels; i++ {
		kernels[i] = make([][]complex128, kernelSize)
		for j := 0; j < kernelSize; j++ {
			kernels[i][j] = make([]complex128, kernelSize)
		}
	}

	f, error := ioutil.ReadFile(filename)
	if error != nil {
		panic(error)
	}

	lines := strings.Split(string(f), "\n")
	for i, l := range lines {
		values := strings.Split(l, " ")
		for j, v := range values {
			value, error := strconv.ParseFloat(v, 64)
			if error != nil {
				panic(error)
			}
			kernels[i][j/kernelSize][j%kernelSize] = complex(float64(value), 0)
		}
	}
}

func readFCData(filename string, insize int, outsize int) (FCData [][]complex128) {
	FCData = make([][]complex128, insize)
	for i := 0; i < insize; i++ {
		FCData[i] = make([]complex128, outsize)
	}

	f, error := ioutil.ReadFile(filename)
	if error != nil {
		panic(error)
	}

	lines := strings.Split(string(f), "\n")
	for i, l := range lines {
		values := strings.Split(l, " ")
		for j, v := range values {
			value, error := strconv.ParseFloat(v, 64)
			if error != nil {
				panic(error)
			}
			FCData[i][j] = complex(float64(value), 0)
		}
	}

	return FCData
}

func readBData(filename string, size int) (BData []complex128) {
	BData = make([]complex128, size)

	f, error := ioutil.ReadFile(filename)
	if error != nil {
		panic(error)
	}

	line := strings.Split(string(f), " ")
	for i, v := range line {
		value, error := strconv.ParseFloat(v, 64)
		if error != nil {
			panic(error)
		}
		BData[i] = complex(float64(value), 0)
	}

	return BData
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Encryption Functions ////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

func encryptImage(testContext *testParams, id string, imageIndex int) (ctOut *mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		encodedImage := make([]complex128, 8192)

		for k := 0; k < numKernels; k++ {
			for i := 0; i < blockSize; i++ {
				for j := 0; j < blockSize; j++ {
					index := blockSize*blockSize*k + blockSize*i + j
					encodedImage[index] = images[imageIndex][2*i][2*j]

					index += 1024
					encodedImage[index] = images[imageIndex][2*i][2*j+1]

					index += 1024
					encodedImage[index] = images[imageIndex][2*i+1][2*j]

					index += 1024
					encodedImage[index] = images[imageIndex][2*i+1][2*j+1]
				}
			}
		}

		for i := 0; i < 4096; i++ {
			encodedImage[i+4096] = encodedImage[i]
		}

		msg := mkckks.NewMessage(testContext.params)
		msg.Value = encodedImage
		ctOut = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}

func encryptKernels(testContext *testParams, id string) (ctOut []*mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		encodedKernels := make([][]complex128, kernelSize)
		for i := 0; i < kernelSize; i++ {
			encodedKernels[i] = make([]complex128, 8192)
		}

		for i := 0; i < numKernels; i++ {
			for j := 0; j < convOutSize; j++ {
				for k := 0; k < convOutSize; k++ {
					index := blockSize*blockSize*i + blockSize*j + k
					encodedKernels[0][index] = kernels[i][0][0]
					encodedKernels[1][index] = kernels[i][0][2]
					encodedKernels[2][index] = kernels[i][2][0]
					encodedKernels[3][index] = kernels[i][2][2]

					index += 1024
					encodedKernels[0][index] = kernels[i][0][1]
					encodedKernels[1][index] = kernels[i][0][3]
					encodedKernels[2][index] = kernels[i][2][1]
					encodedKernels[3][index] = kernels[i][2][3]

					index += 1024
					encodedKernels[0][index] = kernels[i][1][0]
					encodedKernels[1][index] = kernels[i][1][2]
					encodedKernels[2][index] = kernels[i][3][0]
					encodedKernels[3][index] = kernels[i][3][2]

					index += 1024
					encodedKernels[0][index] = kernels[i][1][1]
					encodedKernels[1][index] = kernels[i][1][3]
					encodedKernels[2][index] = kernels[i][3][1]
					encodedKernels[3][index] = kernels[i][3][3]
				}
			}
		}

		for i := 0; i < 4; i++ {
			for j := 0; j < 4096; j++ {
				encodedKernels[i][j+4096] = encodedKernels[i][j]
			}
		}

		ctOut = make([]*mkckks.Ciphertext, kernelSize)
		for i := 0; i < kernelSize; i++ {
			msg := mkckks.NewMessage(testContext.params)
			msg.Value = encodedKernels[i]
			ctOut[i] = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
		}
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}

func encryptFC1(testContext *testParams, id string) (ctOut []*mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		tempFC1 := make([][]complex128, numFCUnit)
		for i := 0; i < numFCUnit; i++ {
			tempFC1[i] = make([]complex128, 1024)
		}

		encodedFC1 := make([][]complex128, 8)
		for i := 0; i < 8; i++ {
			encodedFC1[i] = make([]complex128, 8192)
		}

		// FC1: 5*13*13 x 64
		// tempFC1: 64 x 1024
		for i := 0; i < numKernels; i++ {
			for j := 0; j < convOutSize; j++ {
				for k := 0; k < convOutSize; k++ {
					for l := 0; l < numFCUnit; l++ {
						tempFC1[l][blockSize*blockSize*i+blockSize*j+k] = FC1[i+numKernels*(j*convOutSize+k)][l]
					}
				}
			}
		}

		// encodedFC1: 8 x 8192(64*128)
		for i := 0; i < 8; i++ {
			for j := 0; j < 64; j++ {
				for k := 0; k < 128; k++ {
					encodedFC1[i][128*j+k] = tempFC1[j][128*((i+j)%8)+k]
				}
			}
		}

		ctOut = make([]*mkckks.Ciphertext, 8)
		for i := 0; i < 8; i++ {
			msg := mkckks.NewMessage(testContext.params)
			msg.Value = encodedFC1[i]
			ctOut[i] = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
		}
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}

func encryptFC2(testContext *testParams, id string) (ctOut *mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		numSlots := testContext.params.Slots()
		encodedFC2 := make([]complex128, numSlots)

		numRows := 10
		numColumns := 64
		for i := 0; i < numSlots; i++ {
			x := i / gap
			y := i % gap
			if y < numRows && x < numColumns {
				encodedFC2[i] = FC2[x][y]
			}
		}

		msg := mkckks.NewMessage(testContext.params)
		msg.Value = encodedFC2
		ctOut = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}

func encryptB1(testContext *testParams, id string) (ctOut *mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		numSlots := testContext.params.Slots()
		encodedB1 := make([]complex128, numSlots)
		for i := 0; i < numFCUnit; i++ {
			encodedB1[i*gap] = B1[i]
		}

		msg := mkckks.NewMessage(testContext.params)
		msg.Value = encodedB1
		ctOut = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}

func encryptB2(testContext *testParams, id string) (ctOut *mkckks.Ciphertext) {
	if testContext.encryptor != nil {
		numSlots := testContext.params.Slots()
		encodedB2 := make([]complex128, numSlots)
		for i := 0; i < numClass; i++ {
			encodedB2[i] = B2[i]
		}

		msg := mkckks.NewMessage(testContext.params)
		msg.Value = encodedB2
		ctOut = testContext.encryptor.EncryptMsgNew(msg, testContext.pkSet.GetPublicKey(id))
	} else {
		panic("cannot encrypt image: encryptor is not initialized")
	}
	return ctOut
}
