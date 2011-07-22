//----------------------------------------------------------------------
#include <jni.h>
#include <sys/time.h>
#include "simul1d_bpush1d.h"
#include "bpush1lib_f.h"
#include "globals.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

// Java wrappers for 1d PIC Fortran77 library bpush1lib.f

JNIEXPORT void JNICALL
Java_simul1d_bpush1d_djpost(JNIEnv *env, jclass cls, jdoubleArray part,
                            jdoubleArray cu, jint nop, jdouble qm,
                            jdouble dt, jdoubleArray tdjpost, jint nx,
                            jint ipbc, jint idimp, jint ndim, 
                            jint inorder, jint djopt) {
   jdouble *ppart = NULL, *pcu = NULL, *ptdjpost = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,cu)/ndim;
   ptdjpost = (*env)->GetDoubleArrayElements(env,tdjpost,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         if (djopt==LOOKAHEAD) {
            gsjpost1l(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
         else if (djopt==VECTOR) {
            gsjpost1xl(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
         else {
            gjpost1l(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
      }
      else {
         if (djopt==LOOKAHEAD) {
            gsjpost1(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
         else if (djopt==VECTOR) {
            gsjpost1x(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
         else {
            gjpost1(ppart,pcu,qm,dt,nop,idimp,nx,nxv,ipbc);
         }
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
// record time
   dtimer(&dtime,&itime,1);
   ptdjpost[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tdjpost,ptdjpost,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_bpush1d_push3(JNIEnv *env, jclass cls, jdoubleArray part,
                           jdoubleArray fxyz, jdoubleArray byz,
                           jdouble omx, jint nop, jdouble qbm,
                           jdouble dt, jdouble dtc, jdoubleArray ek,
                           jdoubleArray tpush, jint nx, jint ipbc,
                           jint idimp, jint ndim, jint inorder,
                           jint popt) {
   jdouble *ppart = NULL, *pfxyz = NULL, *pbyz = NULL, *ptpush = NULL;
   jdouble *pek = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,byz)/ndim;
   ptpush = (*env)->GetDoubleArrayElements(env,tpush,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pfxyz = (*env)->GetDoubleArrayElements(env,fxyz,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pek = (*env)->GetDoubleArrayElements(env,ek,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         if (popt==LOOKAHEAD) {
            gsbpush13l(ppart,pfxyz,pbyz,omx,qbm,dt,dtc,pek,idimp,nop,nx,
                       nxv,ipbc);
         }
         else {
            gbpush13l(ppart,pfxyz,pbyz,omx,qbm,dt,dtc,pek,idimp,nop,nx,
                      nxv,ipbc);
         }
      }
      else {
         if (popt==LOOKAHEAD) {
            gsbpush13(ppart,pfxyz,pbyz,omx,qbm,dt,dtc,pek,idimp,nop,nx,
                      nxv,ipbc);
         }
         else {
            gbpush13(ppart,pfxyz,pbyz,omx,qbm,dt,dtc,pek,idimp,nop,nx,
                     nxv,ipbc);
         }
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,fxyz,pfxyz,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ek,pek,0);
// record time
   dtimer(&dtime,&itime,1);
   ptpush[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tpush,ptpush,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_bpush1d_retard(JNIEnv *env, jclass cls, jdoubleArray part,
                            jint nop, jdouble dtc, jint nx, jint ipbc,
                            jint idimp) {
   jdouble *ppart = NULL;
   jboolean isc;
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
//
   retard1(ppart,dtc,idimp,nop,nx,ipbc);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   return;
}
