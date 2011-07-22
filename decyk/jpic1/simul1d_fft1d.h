/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class simul1d_fft1d */

#ifndef _Included_simul1d_fft1d
#define _Included_simul1d_fft1d
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     simul1d_fft1d
 * Method:    fft_init
 * Signature: ([I[DI)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fft_1init
  (JNIEnv *, jclass, jintArray, jdoubleArray, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fft1
 * Signature: ([DI[I[D[DII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fft1
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fft2
 * Signature: ([DI[I[D[DIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fft2
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fst_init
 * Signature: ([I[DI)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fst_1init
  (JNIEnv *, jclass, jintArray, jdoubleArray, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fst1
 * Signature: ([DI[I[D[DII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fst1
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fct1
 * Signature: ([DI[I[D[DII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fct1
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fst2
 * Signature: ([DI[I[D[DIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fst2
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fct2
 * Signature: ([DI[I[D[DIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fct2
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint, jint);

/*
 * Class:     simul1d_fft1d
 * Method:    fcst
 * Signature: ([DI[I[D[DIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_fft1d_fcst
  (JNIEnv *, jclass, jdoubleArray, jint, jintArray, jdoubleArray, jdoubleArray, jint, jint, jint);

#ifdef __cplusplus
}
#endif
#endif