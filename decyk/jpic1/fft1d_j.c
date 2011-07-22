//----------------------------------------------------------------------
#include <jni.h>
#include <sys/time.h>
#include "simul1d_fft1d.h"
#include "fft1lib_f.h"
#include "globals.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

// Java wrappers for 1d PIC Fortran77 library fft1lib.f

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fft_1init(JNIEnv *env, jclass cls, jintArray mixup,
                             jdoubleArray sct, jint indx) {
   jdouble *psct = NULL, *f = NULL, *t = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   int isign = 0, nxd = 1;
   int nxhd = (*env)->GetArrayLength(env,mixup);
// Extract arrays from Java
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psct = (*env)->GetDoubleArrayElements(env,sct,&isc);
//
   fft1rx(f,(double complex*)t,isign,(int *)pmixup,(double complex*)psct,
          indx,nxd,nxhd);
//
// Release the arrays back to Java
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sct,psct,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fft1(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sct,
                        jdoubleArray tfft, jint indx, jint order) {
   jdouble *pf = NULL, *psct = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f);
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psct = (*env)->GetDoubleArrayElements(env,sct,&isc);
   if (order==LINEAR) {
      fft1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psct,indx,
             nxd,nxhd);
   }
   else {
      fft1rx(&pf[1],t,isign,(int *)pmixup,(double complex*)psct,indx,
             nxd,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sct,psct,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fft2(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sct,
                        jdoubleArray tfft, jint indx, jint ndim,
                        jint order) {
   jdouble *pf = NULL, *psct = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f)/ndim;
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[ndim*nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psct = (*env)->GetDoubleArrayElements(env,sct,&isc);
   switch (ndim) {
   case 1:
      if (order==LINEAR) {
          fft1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psct,
                 indx,nxd,nxhd);
      }
      else {
          fft1rx(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psct,
                 indx,nxd,nxhd);
      }
      break;
   case 2:
      if (order==LINEAR) {
         fft1r2(&pf[0],t,isign,(int *)pmixup,(double complex*)psct,
                indx,nxd,nxhd);
      }
      else {
         fft1r2(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psct,
                indx,nxd,nxhd);
      }
      break;
   case 3:
      if (order==LINEAR) {
         fft1r3(&pf[0],t,isign,(int *)pmixup,(double complex*)psct,
                indx,nxd,nxhd);
      }
      else {
         fft1r3(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psct,
                indx,nxd,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sct,psct,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fst_1init(JNIEnv *env, jclass cls, jintArray mixup,
                             jdoubleArray sctd, jint indx) {
   jdouble *psctd = NULL, *f = NULL, *t = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   int isign = 0;
   int nxd = (*env)->GetArrayLength(env,sctd);
   int nxhd = (*env)->GetArrayLength(env,mixup);
// Extract arrays from Java
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
//
   fst1rx(f,(double complex*)t,isign,(int *)pmixup,
          (double complex*)psctd,indx,nxd,nxhd);
//
// Release the arrays back to Java
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fst1(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sctd,
                        jdoubleArray tfft, jint indx, jint order) {
   jdouble *pf = NULL, *psctd = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f);
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
   if (order==LINEAR) {
      fst1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,indx,
             nxd,nxhd);
   }
   else {
      fst1rx(&pf[1],t,isign,(int *)pmixup,(double complex*)psctd,indx,
             nxd,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fct1(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sctd,
                        jdoubleArray tfft, jint indx, jint order) {
   jdouble *pf = NULL, *psctd = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f);
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
   if (order==LINEAR) {
      fct1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,indx,
             nxd,nxhd);
   }
   else {
      fct1rx(&pf[1],t,isign,(int *)pmixup,(double complex*)psctd,indx,
             nxd,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fst2(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sctd,
                        jdoubleArray tfft, jint indx, jint ndim,
                        jint order) {
   jdouble *pf = NULL, *psctd = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f)/ndim;
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[ndim*nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
   switch (ndim) {
   case 1:
      if (order==LINEAR) {
         fst1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      else {
         fst1rx(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      break;
   case 2:
      if (order==LINEAR) {
         fst1r2(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      else {
         fst1r2(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fct2(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sctd,
                        jdoubleArray tfft, jint indx, jint ndim,
                        jint order) {
   jdouble *pf = NULL, *psctd = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f)/ndim;
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[ndim*nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
   switch (ndim) {
   case 1:
      if (order==LINEAR) {
         fct1rx(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      else {
         fct1rx(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      break;
   case 2:
      if (order==LINEAR) {
         fct1r2(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,
                indx,nxd,nxhd);
      }
      else {
         fct1r2(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psctd,
               indx,nxd,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_fft1d_fcst(JNIEnv *env, jclass cls, jdoubleArray f,
                        jint isign, jintArray mixup, jdoubleArray sctd,
                        jdoubleArray tfft, jint indx, jint ndim,
                        jint order) {
   jdouble *pf = NULL, *psctd = NULL, *ptfft = NULL;
   jint *pmixup = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxd = (*env)->GetArrayLength(env,f)/ndim;
   int nxhd = (*env)->GetArrayLength(env,mixup);
   double complex t[ndim*nxhd];
   ptfft = (*env)->GetDoubleArrayElements(env,tfft,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pmixup = (*env)->GetIntArrayElements(env,mixup,&isc);
   psctd = (*env)->GetDoubleArrayElements(env,sctd,&isc);
   switch (ndim) {
   case 3:
      if (order==LINEAR) {
         fcst1r3(&pf[0],t,isign,(int *)pmixup,(double complex*)psctd,
                 indx,nxd,nxhd);
      }
      else {
         fcst1r3(&pf[ndim],t,isign,(int *)pmixup,(double complex*)psctd,
                 indx,nxd,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,mixup,pmixup,0);
   (*env)->ReleaseDoubleArrayElements(env,sctd,psctd,0);
// record time
   dtimer(&dtime,&itime,1);
   ptfft[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tfft,ptfft,0);
   return;
}

