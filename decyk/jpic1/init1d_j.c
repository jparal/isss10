//----------------------------------------------------------------------
#include <jni.h>
#include "simul1d_init1d.h"
#include "init1lib_f.h"
#include "globals.h"

// Java wrappers for 1d PIC Fortran77 library init1lib.f

JNIEXPORT void JNICALL
Java_simul1d_init1d_distr(JNIEnv *env, jclass cls, jdoubleArray part,
                          jint nstart, jint nop, jdouble vtx,
                          jdouble vdx, jint npx, jint nx, jint ipbc,
                          jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
   int ns = idimp*(nstart - 1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
//
   distr1(&ppart[ns],vtx,vdx,npx,idimp,nop,nx,ipbc);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_init1d_distrh(JNIEnv *env, jclass cls, jdoubleArray part,
                           jint nstart, jint nop, jdouble vtx,
                           jdouble vty, jdouble vtz, jdouble vdx,
                           jdouble vdy, jdouble vdz, jint npx, jint nx,
                           jint ipbc, jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
   int ns = idimp*(nstart - 1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
//
   distr1h(&ppart[ns],vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop,nx,ipbc);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_init1d_fdistr(JNIEnv *env, jclass cls, jdoubleArray part,
                           jint nstart, jint nop, jdouble ampx,
                           jdouble scalex, jdouble shiftx, jint npx,
                           jint nx, jint ndpro, jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
   int ns = idimp*(nstart - 1);
   int ierr = 0;
   double sxi = 0.0, zero = 0.0;
   if (scalex != 0.0)
      sxi = 1.0/scalex;
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
// uniform density
   if (ndpro==0) {
      fdistr1(&ppart[ns],fldistr1_,zero,zero,zero,npx,idimp,nop,nx,
              &ierr);
   }
// linear density
   else if (ndpro==1) {
      fdistr1(&ppart[ns],fldistr1_,ampx,sxi,shiftx,npx,idimp,nop,nx,
              &ierr);
   }
// sinusoidal density
   else if (ndpro==2) {
      fdistr1(&ppart[ns],fsdistr1_,ampx,sxi,shiftx,npx,idimp,nop,nx,
              &ierr);
   }
// gaussian density
   else if (ndpro==3) {
      fdistr1(&ppart[ns],fgdistr1_,ampx,sxi,shiftx,npx,idimp,nop,nx,
              &ierr);
   }
// hyperbolic secant squared density
   else if (ndpro==4) {
      fdistr1(&ppart[ns],fhdistr1_,ampx,sxi,shiftx,npx,idimp,nop,nx,
              &ierr);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_init1d_vdistr(JNIEnv *env, jclass cls, jdoubleArray part,
                           jint nstart, jint nop, jdouble vtx,
                           jdouble vdx, jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
   int ns = idimp*(nstart - 1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
//
   vdistr1(&ppart[ns],vtx,vdx,idimp,nop);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_init1d_vdistrh(JNIEnv *env, jclass cls, jdoubleArray part,
                            jint nstart, jint nop, jdouble vtx,
                            jdouble vty, jdouble vtz, jdouble vdx,
                            jdouble vdy, jdouble vdz, jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
   int ns = idimp*(nstart - 1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
//
   vdistr1h(&ppart[ns],vtx,vty,vtz,vdx,vdy,vdz,idimp,nop);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_init1d_bdistr(JNIEnv *env, jclass cls, jdoubleArray part,
                           jdoubleArray byz, jint nop, jdouble qbm,
                           jint nx, jint ipbc, jint idimp, jint ndim,
                           jint inorder) {
   jdouble *ppart = NULL, *pbyz = NULL;
   jboolean isc;
   int nxv = (*env)->GetArrayLength(env,byz)/ndim;
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         gbdistr1l(ppart,&pbyz[0],qbm,idimp,nop,nx,nxv,ipbc);
      }
      else {
         gbdistr1l(ppart,&pbyz[0],qbm,idimp,nop,nx,nxv,ipbc);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   return;
}


