/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class simul1d_bpush1d */

#ifndef _Included_simul1d_bpush1d
#define _Included_simul1d_bpush1d
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     simul1d_bpush1d
 * Method:    djpost
 * Signature: ([D[DIDD[DIIIIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_bpush1d_djpost
  (JNIEnv *, jclass, jdoubleArray, jdoubleArray, jint, jdouble, jdouble, jdoubleArray, jint, jint, jint, jint, jint, jint);

/*
 * Class:     simul1d_bpush1d
 * Method:    push3
 * Signature: ([D[D[DDIDDD[D[DIIIIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_bpush1d_push3
  (JNIEnv *, jclass, jdoubleArray, jdoubleArray, jdoubleArray, jdouble, jint, jdouble, jdouble, jdouble, jdoubleArray, jdoubleArray, jint, jint, jint, jint, jint, jint);

/*
 * Class:     simul1d_bpush1d
 * Method:    retard
 * Signature: ([DIDIII)V
 */
JNIEXPORT void JNICALL Java_simul1d_bpush1d_retard
  (JNIEnv *, jclass, jdoubleArray, jint, jdouble, jint, jint, jint);

#ifdef __cplusplus
}
#endif
#endif
