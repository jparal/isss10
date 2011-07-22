//----------------------------------------------------------------------
#include <jni.h>
#include <sys/time.h>
#include "simul1d_push1d.h"
#include "push1lib_f.h"
#include "globals.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

// Java wrappers for 1d PIC Fortran77 library push1lib.f

JNIEXPORT void JNICALL
Java_simul1d_push1d_dpost(JNIEnv *env, jclass cls, jdoubleArray part,
                          jdoubleArray q, jint nop, jdouble qm,
                          jdoubleArray tdpost, jint idimp, jint inorder,
                          jint dopt) {
   jdouble *ppart = NULL, *pq = NULL, *ptdpost = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,q);
   ptdpost = (*env)->GetDoubleArrayElements(env,tdpost,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pq = (*env)->GetDoubleArrayElements(env,q,&isc);
   if (inorder==LINEAR) {
      if (dopt==LOOKAHEAD) {
         gspost1l(ppart,pq,qm,nop,idimp,nxv);
      }
      else if (dopt==VECTOR) {
         gspost1xl(ppart,pq,qm,nop,idimp,nxv);
      }
      else {
         gpost1l(ppart,pq,qm,nop,idimp,nxv);
      }
   }
   else {
      if (dopt==LOOKAHEAD) {
         gspost1(ppart,pq,qm,nop,idimp,nxv);
      }
      else if (dopt==VECTOR) {
         gspost1x(ppart,pq,qm,nop,idimp,nxv);
      }
      else {
         gpost1(ppart,pq,qm,nop,idimp,nxv);
      }
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,q,pq,0);
// record time
   dtimer(&dtime,&itime,1);
   ptdpost[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tdpost,ptdpost,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_push1d_push(JNIEnv *env, jclass cls, jdoubleArray part,
                         jdoubleArray fx, jint nop, jdouble qbm,
                         jdouble dt, jdoubleArray ek,
                         jdoubleArray tpush, jint nx, jint ipbc,
                         jint idimp, jint inorder, jint popt) {
   jdouble *ppart = NULL, *pfx = NULL, *ptpush = NULL, *pek = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxv = (*env)->GetArrayLength(env,fx);
   ptpush = (*env)->GetDoubleArrayElements(env,tpush,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   pek = (*env)->GetDoubleArrayElements(env,ek,&isc);
   if (inorder==LINEAR) {
      if (popt==LOOKAHEAD) {
         gspush1l(ppart,pfx,qbm,dt,pek,idimp,nop,nx,nxv,ipbc);
      }
      else {
         gpush1l(ppart,pfx,qbm,dt,pek,idimp,nop,nx,nxv,ipbc);
      }
   }
   else {
      if (popt==LOOKAHEAD) {
         gspush1(ppart,pfx,qbm,dt,pek,idimp,nop,nx,nxv,ipbc);
      }
      else {
         gpush1(ppart,pfx,qbm,dt,pek,idimp,nop,nx,nxv,ipbc);
      }
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   (*env)->ReleaseDoubleArrayElements(env,ek,pek,0);
// record time
   dtimer(&dtime,&itime,1);
   ptpush[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tpush,ptpush,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_push1d_sortp(JNIEnv *env, jclass cls, jdoubleArray part,
                          jdoubleArray pt, jintArray ip, jint nop,
                          jintArray npic, jdoubleArray tsort,
                          jint idimp, jint inorder) {
   jdouble *ppart = NULL, *ppt = NULL, *ptsort = NULL;
   jint *pip = NULL, *pnpic = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nx1 = (*env)->GetArrayLength(env,npic);
   ptsort = (*env)->GetDoubleArrayElements(env,tsort,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   ppt = (*env)->GetDoubleArrayElements(env,pt,&isc);
   pip = (*env)->GetIntArrayElements(env,ip,&isc);
   pnpic = (*env)->GetIntArrayElements(env,npic,&isc);
   if (inorder==LINEAR) {
      sortp1xl(ppart,ppt,(int *)pip,(int *)pnpic,idimp,nop,nx1);
   }
   else {
      sortp1x(ppart,ppt,(int *)pip,(int *)pnpic,idimp,nop,nx1);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,pt,ppt,0);
   (*env)->ReleaseIntArrayElements(env,ip,pip,0);
   (*env)->ReleaseIntArrayElements(env,npic,pnpic,0);
// record time
   dtimer(&dtime,&itime,1);
   ptsort[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tsort,ptsort,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_push1d_dpostgl(JNIEnv *env, jclass cls, jdoubleArray part,
                            jdoubleArray q, jint nop, jdouble qm,
                            jint nx, jint nxh, jdoubleArray tdpost,
                            jint idimp) {
   jdouble *ppart = NULL, *pq = NULL, *ptdpost = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxvh = (*env)->GetArrayLength(env,q)/2;
   double complex sctx[nxvh];
   ptdpost = (*env)->GetDoubleArrayElements(env,tdpost,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pq = (*env)->GetDoubleArrayElements(env,q,&isc);
//
   dpost1gl(ppart,(double complex *)q,sctx,qm,nop,idimp,nx,nxh,nxvh);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,q,pq,0);
// record time
   dtimer(&dtime,&itime,1);
   ptdpost[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tdpost,ptdpost,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_push1d_pushgl(JNIEnv *env, jclass cls, jdoubleArray part,
                           jdoubleArray fx, jint nop, jdouble qbm,
                           jdouble dt, jdoubleArray ek, jint nx,
                           jint nxh, jdoubleArray tpush, jint idimp) {
   jdouble *ppart = NULL, *pfx = NULL, *ptpush = NULL, *pek = NULL;
   jboolean isc;
   struct timeval itime;
   double dtime;
   int nxvh = (*env)->GetArrayLength(env,fx)/2;
   double complex sctx[nxvh];
   ptpush = (*env)->GetDoubleArrayElements(env,tpush,&isc);
// initialize timer
   dtimer(&dtime,&itime,-1);
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   pek = (*env)->GetDoubleArrayElements(env,ek,&isc);
//
   push1gl(ppart,(double complex *)pfx,sctx,qbm,dt,pek,idimp,nop,nx,nxh,
           nxvh);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   (*env)->ReleaseDoubleArrayElements(env,ek,pek,0);
// record time
   dtimer(&dtime,&itime,1);
   ptpush[0] += dtime;
   (*env)->ReleaseDoubleArrayElements(env,tpush,ptpush,0);
   return;
}
