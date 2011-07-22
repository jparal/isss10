//----------------------------------------------------------------------
#include <jni.h>
#include <sys/time.h>
#include "simul1d_bpush1d.h"
#include "dpush1lib_f.h"
#include "globals.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

JNIEXPORT void JNICALL
Java_simul1d_dpush1d_dmjpost(JNIEnv *env, jclass cls, jdoubleArray part,
                             jdoubleArray amu, jint nop, jdouble qm,
                             jdoubleArray tdcjpost, jint idimp,
                             jint ndim, jint inorder, jint djopt) {
   jdouble *ppart = NULL, *pamu = NULL, *ptdcjpost = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,amu)/ndim;
   ptdcjpost = (*env)->GetDoubleArrayElements(env,tdcjpost,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pamu = (*env)->GetDoubleArrayElements(env,amu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         if (djopt==LOOKAHEAD) {
            gsmjpost1l(ppart,pamu,qm,nop,idimp,nxv);
         }
         else {
            gmjpost1l(ppart,pamu,qm,nop,idimp,nxv);
         }
      }
      else {
         if (djopt==LOOKAHEAD) {
            gsmjpost1(ppart,pamu,qm,nop,idimp,nxv);
         }
         else {
            gmjpost1(ppart,pamu,qm,nop,idimp,nxv);
         }
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,amu,pamu,0);
// record time
   dtimer(&dtime,&itime,1);
   ptdcjpost[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tdcjpost,ptdcjpost,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_dpush1d_dcjpost(JNIEnv *env, jclass cls, jdoubleArray part,
                             jdoubleArray fxyz, jdoubleArray byz,
                             jdoubleArray cu, jdoubleArray dcu,
                             jdoubleArray amu, jdouble omx, jint nop,
                             jdouble qm, jdouble qbm, jdouble dt,
                             jdoubleArray tdcjpost, jint idimp,
                             jint ndim, jint inorder, int djopt) {
   jdouble *ppart = NULL, *pfxyz = NULL, *pbyz = NULL, *pcu = NULL;
   jdouble *pdcu = NULL, *pamu = NULL, *ptdcjpost = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,byz)/ndim;
   ptdcjpost = (*env)->GetDoubleArrayElements(env,tdcjpost,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pfxyz = (*env)->GetDoubleArrayElements(env,fxyz,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   pdcu = (*env)->GetDoubleArrayElements(env,dcu,&isc);
   pamu = (*env)->GetDoubleArrayElements(env,amu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         if (djopt==LOOKAHEAD) {
            gsdcjpost1l(ppart,pfxyz,pbyz,pcu,pdcu,pamu,omx,qm,qbm,dt,
                        idimp,nop,nxv);
         }
         else {
           gdcjpost1l(ppart,pfxyz,pbyz,pcu,pdcu,pamu,omx,qm,qbm,dt,
                      idimp,nop,nxv);
         }
      }
      else {
         if (djopt==LOOKAHEAD) {
            gsdcjpost1(ppart,pfxyz,pbyz,pcu,pdcu,pamu,omx,qm,qbm,dt,
                       idimp,nop,nxv);
         }
         else {
            gdcjpost1(ppart,pfxyz,pbyz,pcu,pdcu,pamu,omx,qm,qbm,dt,
                      idimp,nop,nxv);
         }
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,fxyz,pfxyz,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   (*env)->ReleaseDoubleArrayElements(env,dcu,pdcu,0);
   (*env)->ReleaseDoubleArrayElements(env,amu,pamu,0);
// record time
   dtimer(&dtime,&itime,1);
   ptdcjpost[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tdcjpost,ptdcjpost,0);
   return;
}

JNIEXPORT void JNICALL
//----------------------------------------------------------------------
Java_simul1d_diag1d_rwnml1(JNIEnv *env, jclass cls, jstring fname,
                           jstring name, jobjectArray cn,
                           jdoubleArray ddata, jint iunit, jint isign,
                           jint nl, jintArray ml, jintArray ierr) {
   return;
}
