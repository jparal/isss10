//----------------------------------------------------------------------
#include <jni.h>
#include "simul1d_field1d.h"
#include "field1lib_f.h"
#include "globals.h"

// Java wrappers for 1d PIC Fortran77 library field1lib.f

JNIEXPORT void JNICALL
Java_simul1d_field1d_cguard(JNIEnv *env, jclass cls, jdoubleArray fxy,
                            jint nx, jint ndim, jint inorder) {
   jdouble *pfxy = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,fxy)/ndim;
// Extract array from Java
   pfxy = (*env)->GetDoubleArrayElements(env,fxy,&isc);
   switch (ndim) {
   case 1:
      if (inorder==LINEAR) {
         dguard1l(pfxy,nx,nxe);
      }
      else {
         dguard1(pfxy,nx,nxe);
      }
      break;
   case 2:
      if (inorder==LINEAR) {
         cguard1l(pfxy,nx,nxe);
      }
      else {
         cguard1(pfxy,nx,nxe);
      }
      break;
   case 3:
      if (inorder==LINEAR) {
         bguard1l(pfxy,nx,nxe);
      }
      else {
         bguard1(pfxy,nx,nxe);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fxy,pfxy,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_dguard(JNIEnv *env, jclass cls, jdoubleArray fx,
                            jint nx, jint inorder) {
   jdouble *pfx = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,fx);
// Extract array from Java
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   if (inorder==LINEAR) {
      dguard1l(pfx,nx,nxe);
   }
   else {
      dguard1(pfx,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_scguard(JNIEnv *env, jclass cls, jdoubleArray cu,
                             jdouble yj0, jdouble zj0, jint nx,
                             jint inorder) {
   jdouble *pcu = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,cu)/2;
// Extract array from Java
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   if (inorder==LINEAR) {
      scguard1l(pcu,yj0,zj0,nx,nxe);
   }
   else {
      scguard1(pcu,yj0,zj0,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_sguard(JNIEnv *env, jclass cls, jdoubleArray q,
                            jdouble qi0, jint nx, jint inorder) {
   jdouble *pq = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,q);
// Extract array from Java
   pq = (*env)->GetDoubleArrayElements(env,q,&isc);
   if (inorder==LINEAR) {
      sguard1l(pq,qi0,nx,nxe);
   }
   else {
      sguard1(pq,qi0,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,q,pq,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_acguard(JNIEnv *env, jclass cls, jdoubleArray cu,
                             jint nx, jint inorder) {
   jdouble *pcu = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,cu)/2;
// Extract array from Java
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   if (inorder==LINEAR) {
      acguard1l(pcu,nx,nxe);
   }
   else {
      acguard1(pcu,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_aguard(JNIEnv *env, jclass cls, jdoubleArray q,
                            jint nx, jint inorder) {
   jdouble *pq = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,q);
// Extract array from Java
   pq = (*env)->GetDoubleArrayElements(env,q,&isc);
   if (inorder==LINEAR) {
      aguard1l(pq,nx,nxe);
   }
   else {
      aguard1(pq,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,q,pq,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_pois_1init(JNIEnv *env, jclass cls,
                                jdoubleArray ffc, jdouble ax,
                                jdouble affp, jint nx) {
   jdouble *pffc = NULL;
   jboolean isc;
   double *q = NULL, *fx = NULL, *we = NULL;
   int isign = 0;
// Extract array from Java
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
//
   poisp1(q,fx,isign,pffc,ax,affp,we,nx);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_pois(JNIEnv *env, jclass cls, jdoubleArray q,
                          jdoubleArray fx, jint isign, jdoubleArray ffc,
                          jdoubleArray we, jint nx, jint inorder) {
   jdouble *pq = NULL, *pfx = NULL, *pffc = NULL, *pwe = NULL;
   jboolean isc;
   double ax = 0.0, affp = 0.0;
// Extract array from Java
   pq = (*env)->GetDoubleArrayElements(env,q,&isc);
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   pwe = (*env)->GetDoubleArrayElements(env,we,&isc);
   if (inorder==LINEAR) {
      poisp1(&pq[0],&pfx[0],isign,pffc,ax,affp,pwe,nx);
   }
   else {
      poisp1(&pq[1],&pfx[1],isign,pffc,ax,affp,pwe,nx);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,q,pq,0);
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   (*env)->ReleaseDoubleArrayElements(env,we,pwe,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_bpois(JNIEnv *env, jclass cls, jdoubleArray cu,
                           jdoubleArray byz, jdoubleArray ffc,
                           jdouble ci, jdoubleArray wm, jint nx,
                           jint ndim, jint inorder) {
   jdouble *pcu = NULL, *pbyz = NULL, *pffc = NULL, *pwm = NULL;
   jboolean isc;
   int isign = -1;
   double ax = 0.0, affp = 0.0;
   int nxvh = (*env)->GetArrayLength(env,cu)/(2*ndim);
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   pwm = (*env)->GetDoubleArrayElements(env,wm,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         bpois13((double complex*)&pcu[0],(double complex*)&pbyz[0],
                  isign,(double complex*)pffc,ax,affp,ci,pwm,nx,nxvh,
                  nxhd);
     }
     else {
         bpois13((double complex*)&pcu[ndim],
                 (double complex*)&pbyz[ndim],isign,
                 (double complex*)pffc,ax,affp,ci,pwm,nx,nxvh,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   (*env)->ReleaseDoubleArrayElements(env,wm,pwm,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_apois(JNIEnv *env, jclass cls, jdoubleArray cu,
                           jdoubleArray ayz, jdoubleArray ffc,
                           jdouble ci, jdoubleArray wm, jint nx,
                           jint ndim, jint inorder) {
   jdouble *pcu = NULL, *payz = NULL, *pffc = NULL, *pwm = NULL;
   jboolean isc;
   int isign = 1;
   double ax = 0.0, affp = 0.0;
   int nxvh = (*env)->GetArrayLength(env,cu)/(2*ndim);
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   payz = (*env)->GetDoubleArrayElements(env,ayz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   pwm = (*env)->GetDoubleArrayElements(env,wm,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         bpois13((double complex*)&pcu[0],(double complex*)&payz[0],
                 isign,(double complex*)pffc,ax,affp,ci,pwm,nx,nxvh,
                 nxhd);
      }
      else {
         bpois13((double complex*)&pcu[ndim],
                 (double complex*)&payz[ndim],isign,
                 (double complex*)pffc,ax,affp,ci,pwm,nx,nxvh,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   (*env)->ReleaseDoubleArrayElements(env,ayz,payz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   (*env)->ReleaseDoubleArrayElements(env,wm,pwm,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_ibpois(JNIEnv *env, jclass cls, jdoubleArray cu,
                            jdoubleArray byz, jdoubleArray ffc,
                            jdouble ci, jdoubleArray wm, jint nx,
                            jint inorder) {
   jdouble *pcu = NULL, *pbyz = NULL, *pffc = NULL, *pwm = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,cu)/4;
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   pwm = (*env)->GetDoubleArrayElements(env,wm,&isc);
   if (inorder==LINEAR) {
      ibpois13((double complex*)&pcu[0],(double complex*)&pbyz[0],
               (double complex*)pffc,ci,pwm,nx,nxvh,nxhd);
   }
   else {
      ibpois13((double complex*)&pcu[2],(double complex*)&pbyz[2],
               (double complex*)pffc,ci,pwm,nx,nxvh,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   (*env)->ReleaseDoubleArrayElements(env,wm,pwm,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_maxwel(JNIEnv *env, jclass cls, jdoubleArray eyz,
                            jdoubleArray byz, jdoubleArray cu,
                            jdoubleArray ffc, jdouble ci, jdouble dt,
                            jdoubleArray wf, jdoubleArray wm, jint nx,
                            jint inorder) {
   jdouble *peyz = NULL, *pbyz = NULL, *pcu = NULL, *pffc = NULL;
   jdouble *pwf = NULL, *pwm = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,cu)/4;
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   pwf = (*env)->GetDoubleArrayElements(env,wf,&isc);
   pwm = (*env)->GetDoubleArrayElements(env,wm,&isc);
   if (inorder==LINEAR) {
      maxwel1((double complex*)peyz,(double complex*)pbyz,
              (double complex*)&pcu[0],(double complex*)pffc,ci,dt,pwf,
              pwm,nx,nxvh,nxhd);
   }
   else {
      maxwel1((double complex*)peyz,(double complex*)pbyz,
              (double complex*)&pcu[2],(double complex*)pffc,ci,dt,pwf,
              pwm,nx,nxvh,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   (*env)->ReleaseDoubleArrayElements(env,wf,pwf,0);
   (*env)->ReleaseDoubleArrayElements(env,wm,pwm,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_emfield(JNIEnv *env, jclass cls, jdoubleArray fxyz,
                             jdoubleArray fx, jdoubleArray eyz,
                             jdoubleArray ffc, jint nx, jint inorder) {
   jdouble *pfxyz = NULL, *pfx = NULL, *peyz = NULL, *pffc = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,fxyz)/6;
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   pfxyz = (*env)->GetDoubleArrayElements(env,fxyz,&isc);
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   if (inorder==LINEAR) {
      emfield1((double complex*)&pfxyz[0],(double complex*)&pfx[0],
               (double complex*)peyz,(double complex*)pffc,nx,nxvh,
               nxhd);
   }
   else {
      emfield1((double complex*)&pfxyz[3],(double complex*)&pfx[1],
               (double complex*)peyz,(double complex*)pffc,nx,nxvh,
               nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fxyz,pfxyz,0);
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_bmfield(JNIEnv *env, jclass cls, jdoubleArray fyz,
                             jdoubleArray eyz, jdoubleArray ffc,
                             jint nx, jint inorder) {
   jdouble *pfyz = NULL, *peyz = NULL, *pffc = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,fyz)/4;
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   pfyz = (*env)->GetDoubleArrayElements(env,fyz,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   if (inorder==LINEAR) {
      bmfield1((double complex*)&pfyz[0],(double complex*)peyz,
               (double complex*)pffc,nx,nxvh,nxhd);
   }
   else {
      bmfield1((double complex*)&pfyz[2],(double complex*)peyz,
               (double complex*)pffc,nx,nxvh,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fyz,pfyz,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_emfieldr(JNIEnv *env, jclass cls,
                              jdoubleArray fxyz, jdoubleArray fx,
                              jdoubleArray eyz, jdoubleArray ffc,
                              jint nx, jint inorder) {
   jdouble *pfxyz = NULL, *pfx = NULL, *peyz = NULL, *pffc = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,fxyz)/3;
   int nxd = (*env)->GetArrayLength(env,ffc);
// Extract array from Java
   pfxyz = (*env)->GetDoubleArrayElements(env,fxyz,&isc);
   pfx = (*env)->GetDoubleArrayElements(env,fx,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
    if (inorder==LINEAR) {
      emfieldr1(&pfxyz[0],&pfx[0],peyz,(double complex*)pffc,nx,nxe,
                nxd);
   }
   else {
      emfieldr1(&pfxyz[3],&pfx[1],peyz,(double complex*)pffc,nx,nxe,
                nxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fxyz,pfxyz,0);
   (*env)->ReleaseDoubleArrayElements(env,fx,pfx,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_bmfieldr(JNIEnv *env, jclass cls, jdoubleArray fyz,
                              jdoubleArray eyz, jdoubleArray ffc,
                              jint nx, jint inorder) {
   jdouble *pfyz = NULL, *peyz = NULL, *pffc = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,fyz)/2;
   int nxd = (*env)->GetArrayLength(env,ffc);
// Extract array from Java
   pfyz = (*env)->GetDoubleArrayElements(env,fyz,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
    if (inorder==LINEAR) {
      bmfieldr1(&pfyz[0],peyz,(double complex*)pffc,nx,nxe,nxd);
   }
   else {
      bmfieldr1(&pfyz[2],peyz,(double complex*)pffc,nx,nxe,nxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,fyz,pfyz,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_avpot(JNIEnv *env, jclass cls, jdoubleArray byz,
                            jdoubleArray ayz, jint nx, jint inorder) {
   jdouble *pbyz = NULL, *payz = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,byz)/4;
// Extract array from Java
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   payz = (*env)->GetDoubleArrayElements(env,ayz,&isc);
   if (inorder==LINEAR) {
      avpot13((double complex*)pbyz,(double complex*)&payz[0],nx,nxvh);
   }
   else {
      avpot13((double complex*)pbyz,(double complex*)&payz[2],nx,nxvh);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ayz,payz,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_avrpot(JNIEnv *env, jclass cls, jdoubleArray ayz,
                            jdoubleArray byz, jdoubleArray ffc,
                            jdouble ci, jint nx, jint inorder) {
   jdouble *payz = NULL, *pbyz = NULL, *pffc = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,ayz)/4;
   int nxhd = (*env)->GetArrayLength(env,ffc)/2;
// Extract array from Java
   payz = (*env)->GetDoubleArrayElements(env,ayz,&isc);
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   pffc = (*env)->GetDoubleArrayElements(env,ffc,&isc);
   if (inorder==LINEAR) {
      avrpot13((double complex*)&payz[0],(double complex*)pbyz,
               (double complex*)pffc,ci,nx,nxvh,nxhd);
   }
   else {
      avrpot13((double complex*)&payz[2],(double complex*)pbyz,
               (double complex*)pffc,ci,nx,nxvh,nxhd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,ayz,payz,0);
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffc,pffc,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_gtmodes(JNIEnv *env, jclass cls, jdoubleArray pot,
                             jdoubleArray pott, jint nx, jint modesx,
                             jint order) {
   jdouble *ppot = NULL, *ppott = NULL;
   jboolean isc;
   int it = 1, nt2 = 2;
   int nxe = (*env)->GetArrayLength(env,pot);
   int modesxd = (*env)->GetArrayLength(env,pott)/2;
// Extract array from Java
   ppot = (*env)->GetDoubleArrayElements(env,pot,&isc);
   ppott = (*env)->GetDoubleArrayElements(env,pott,&isc);
   if (order==LINEAR) {
      gtmodes1(&ppot[0],ppott,nx,it,modesx,nxe,nt2,
               modesxd);
   }
   else {
       gtmodes1(&ppot[1],ppott,nx,it,modesx,nxe,nt2,
                modesxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,pot,ppot,0);
   (*env)->ReleaseDoubleArrayElements(env,pott,ppott,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_ptmodes(JNIEnv *env, jclass cls, jdoubleArray pot,
                             jdoubleArray pott, jint nx, jint modesx,
                             jint order) {
   jdouble *ppot = NULL, *ppott = NULL;
   jboolean isc;
   int it = 1, nt2 = 2;
   int nxe = (*env)->GetArrayLength(env,pot);
   int modesxd = (*env)->GetArrayLength(env,pott)/2;
// Extract array from Java
   ppot = (*env)->GetDoubleArrayElements(env,pot,&isc);
   ppott = (*env)->GetDoubleArrayElements(env,pott,&isc);
   if (order==LINEAR) {
      ptmodes1(&ppot[0],ppott,nx,it,modesx,nxe,nt2,modesxd);
   }
   else {
      ptmodes1(&ppot[1],ppott,nx,it,modesx,nxe,nt2,modesxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,pot,ppot,0);
   (*env)->ReleaseDoubleArrayElements(env,pott,ppott,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_gtvmodes(JNIEnv *env, jclass cls, 
                              jdoubleArray vpot, jdoubleArray vpott,
                              jint nx, jint modesx, jint ndim,
                              jint order) {
   jdouble *pvpot = NULL, *pvpott = NULL;
   jboolean isc;
   int it = 1, nt = 1;
   int nxvh = (*env)->GetArrayLength(env,vpot)/(2*ndim);
   int modesxd = (*env)->GetArrayLength(env,vpott)/(2*ndim);
// Extract array from Java
   pvpot = (*env)->GetDoubleArrayElements(env,vpot,&isc);
   pvpott = (*env)->GetDoubleArrayElements(env,vpott,&isc);
   if (order==LINEAR) {
      gtvmodes1((double complex*)&pvpot[0],(double complex*)pvpott,
                nx,it,modesx,ndim,nxvh,nt,modesxd);
   }
   else {
      gtvmodes1((double complex*)&pvpot[ndim],(double complex*)pvpott,
                nx,it,modesx,ndim,nxvh,nt,modesxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,vpot,pvpot,0);
   (*env)->ReleaseDoubleArrayElements(env,vpott,pvpott,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_ptvmodes(JNIEnv *env, jclass cls,
                              jdoubleArray vpot, jdoubleArray vpott,
                              jint nx, jint modesx, jint ndim,
                              jint order) {
   jdouble *pvpot = NULL, *pvpott = NULL;
   jboolean isc;
   int it = 1, nt = 1;
   int nxvh = (*env)->GetArrayLength(env,vpot)/(2*ndim);
   int modesxd = (*env)->GetArrayLength(env,vpott)/(2*ndim);
// Extract array from Java
   pvpot = (*env)->GetDoubleArrayElements(env,vpot,&isc);
   pvpott = (*env)->GetDoubleArrayElements(env,vpott,&isc);
   if (order==LINEAR) {
      ptvmodes1((double complex*)&pvpot[0],(double complex*)pvpott,nx,
                 it,modesx,ndim,nxvh,nt,modesxd);
   }
   else {
      ptvmodes1((double complex*)&pvpot[ndim],(double complex*)pvpott,
                nx,it,modesx,ndim,nxvh,nt,modesxd);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,vpot,pvpot,0);
   (*env)->ReleaseDoubleArrayElements(env,vpott,pvpott,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_sfguard(JNIEnv *env, jclass cls, jdoubleArray cus,
                             jdoubleArray cu, jdouble q2m0, jint nx,
                             jint ndim, jint inorder) {
   jdouble *pcus = NULL, *pcu = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,cus)/ndim;
// Extract array from Java
   pcus = (*env)->GetDoubleArrayElements(env,cus,&isc);
   pcu = (*env)->GetDoubleArrayElements(env,cu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         scfguard1l(pcus,pcu,q2m0,nx,nxe);
      }
      else {
         scfguard1(pcus,pcu,q2m0,nx,nxe);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,cus,pcus,0);
   (*env)->ReleaseDoubleArrayElements(env,cu,pcu,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_dcuperp(JNIEnv *env, jclass cls, jdoubleArray dcu,
                             jdoubleArray amu, jint nx, jint ndim,
                             jint inorder) {
   jdouble *pdcu = NULL, *pamu = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,dcu)/(2*ndim);
// Extract array from Java
   pdcu = (*env)->GetDoubleArrayElements(env,dcu,&isc);
   pamu = (*env)->GetDoubleArrayElements(env,amu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         dcuperp13((double complex*)&pdcu[0],(double complex*)&pamu[0],
                   nx,nxvh);
      }
      else {
         dcuperp13((double complex*)&pdcu[ndim],
                   (double complex*)&pamu[ndim],nx,nxvh);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,dcu,pdcu,0);
   (*env)->ReleaseDoubleArrayElements(env,amu,pamu,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_adcuperp(JNIEnv *env, jclass cls, jdoubleArray dcu,
                              jdoubleArray amu, jint nx, jint ndim,
                              jint inorder) {
   jdouble *pdcu = NULL, *pamu = NULL;
   jboolean isc;
   int nxvh = (*env)->GetArrayLength(env,dcu)/(2*ndim);
// Extract array from Java
   pdcu = (*env)->GetDoubleArrayElements(env,dcu,&isc);
   pamu = (*env)->GetDoubleArrayElements(env,amu,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         adcuperp13((double complex*)&pdcu[0],(double complex*)&pamu[0],
                     nx,nxvh);
      }
      else {
         adcuperp13((double complex*)&pdcu[ndim],
                    (double complex*)&pamu[ndim],nx,nxvh);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,dcu,pdcu,0);
   (*env)->ReleaseDoubleArrayElements(env,amu,pamu,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_epois_1init(JNIEnv *env, jclass cls,
                                 jdoubleArray ffe, jdouble ax,
                                 jdouble affp, jdouble wp0, jdouble ci,
                                 jint nx) {
   jdouble *pffe = NULL;
   jboolean isc;
   double complex *dcu = NULL, *eyz = NULL;
   double *wf = NULL;
   int isign = 0, nxvh = 1;
   int nxhd = (*env)->GetArrayLength(env,ffe)/2;
// Extract array from Java
   pffe = (*env)->GetDoubleArrayElements(env,ffe,&isc);
//
   epois13(dcu,eyz,isign,(double complex*)pffe,ax,affp,wp0,ci,wf,nx,
           nxvh,nxhd);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,ffe,pffe,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_epois(JNIEnv *env, jclass cls, jdoubleArray dcu,
                           jdoubleArray eyz, jdoubleArray ffe,
                           jdouble ci, jdoubleArray wf, jint nx, 
                           jint ndim, jint inorder) {
   jdouble *pdcu = NULL, *peyz = NULL, *pffe = NULL, *pwf = NULL;
   jboolean isc;
   int isign = -1;
   double ax = 0.0, affp = 0.0, wp0 = 0.0;
   int nxvh = (*env)->GetArrayLength(env,dcu)/(2*ndim);
   int nxhd = (*env)->GetArrayLength(env,ffe)/2;
// Extract array from Java
   pdcu = (*env)->GetDoubleArrayElements(env,dcu,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffe = (*env)->GetDoubleArrayElements(env,ffe,&isc);
   pwf = (*env)->GetDoubleArrayElements(env,wf,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         epois13((double complex*)&pdcu[0],(double complex*)&peyz[0],
                  isign,(double complex*)pffe,ax,affp,wp0,ci,pwf,nx,
                  nxvh,nxhd);
      }
      else {
         epois13((double complex*)&pdcu[ndim],
                 (double complex*)&peyz[ndim],isign,
                 (double complex*)pffe,ax,affp,wp0,ci,pwf,nx,nxvh,nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,dcu,pdcu,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffe,pffe,0);
   (*env)->ReleaseDoubleArrayElements(env,wf,pwf,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_iepois(JNIEnv *env, jclass cls, jdoubleArray dcu,
                            jdoubleArray eyz, jdoubleArray ffe,
                            jdouble ci, jdoubleArray wf, jint nx,
                            jint ndim, jint inorder) {
   jdouble *pdcu = NULL, *peyz = NULL, *pffe = NULL, *pwf = NULL;
   jboolean isc;
   int isign = 1;
   double ax = 0.0, affp = 0.0, wp0 = 0.0;
   int nxvh = (*env)->GetArrayLength(env,dcu)/(2*ndim);
   int nxhd = (*env)->GetArrayLength(env,ffe)/2;
// Extract array from Java
   pdcu = (*env)->GetDoubleArrayElements(env,dcu,&isc);
   peyz = (*env)->GetDoubleArrayElements(env,eyz,&isc);
   pffe = (*env)->GetDoubleArrayElements(env,ffe,&isc);
   pwf = (*env)->GetDoubleArrayElements(env,wf,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         epois13((double complex*)&pdcu[0],(double complex*)peyz,isign,
                 (double complex*)pffe,ax,affp,wp0,ci,pwf,nx,nxvh,nxhd);
      }
      else {
         epois13((double complex*)&pdcu[ndim],(double complex*)peyz,
                 isign,(double complex*)pffe,ax,affp,wp0,ci,pwf,nx,nxvh,
                 nxhd);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,dcu,pdcu,0);
   (*env)->ReleaseDoubleArrayElements(env,eyz,peyz,0);
   (*env)->ReleaseDoubleArrayElements(env,ffe,pffe,0);
   (*env)->ReleaseDoubleArrayElements(env,wf,pwf,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_wpmxn(JNIEnv *env, jclass cls, jdoubleArray qe,
                           jdouble qi0, jdouble qbme,
                           jdoubleArray wpmax, jdoubleArray wpmin,
                           jint nx, jint inorder) {
   jdouble *pqe = NULL, *pwpmax = NULL, *pwpmin = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,qe);
// Extract array from Java
   pqe = (*env)->GetDoubleArrayElements(env,qe,&isc);
   pwpmax = (*env)->GetDoubleArrayElements(env,wpmax,&isc);
   pwpmin = (*env)->GetDoubleArrayElements(env,wpmin,&isc);
   if (inorder==LINEAR) {
      wpmxn1(&pqe[0],qi0,qbme,pwpmax,pwpmin,nx,nxe);
   }
   else {
      wpmxn1(&pqe[1],qi0,qbme,pwpmax,pwpmin,nx,nxe);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,qe,pqe,0);
   (*env)->ReleaseDoubleArrayElements(env,wpmax,pwpmax,0);
   (*env)->ReleaseDoubleArrayElements(env,wpmin,pwpmin,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_baddext(JNIEnv *env, jclass cls, jdoubleArray byz,
                             jdouble omy, jdouble omz, jint nx,
                             jint ndim, jint inorder) {
   jdouble *pbyz = NULL;
   jboolean isc;
   int nxe = (*env)->GetArrayLength(env,byz)/ndim;
// Extract array from Java
   pbyz = (*env)->GetDoubleArrayElements(env,byz,&isc);
   switch (ndim) {
   case 2:
      if (inorder==LINEAR) {
         baddext1(&pbyz[0],omy,omz,nx,nxe);
      }
      else {
         baddext1(&pbyz[ndim],omy,omz,nx,nxe);
      }
      break;
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,byz,pbyz,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_ivrcopy(JNIEnv *env, jclass cls, jdoubleArray f,
                             jdoubleArray g, jint nx, jint ndim,
                             jint inorder) {
   jdouble *pf = NULL, *pg = NULL;
   jboolean isc;
   int nxv = (*env)->GetArrayLength(env,g)/ndim;
   if (nx > nxv)
      return;
// Extract array from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pg = (*env)->GetDoubleArrayElements(env,g,&isc);
   if (inorder==LINEAR) {
      vrcopy1(pf,&pg[0],nx,ndim,nxv);
   }
   else {
      vrcopy1(pf,&pg[ndim],nx,ndim,nxv);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseDoubleArrayElements(env,g,pg,0);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_field1d_ivccopy(JNIEnv *env, jclass cls, jdoubleArray f,
                             jdoubleArray g, jint nx, jint ndim,
                             jint inorder) {
   jdouble *pf = NULL, *pg = NULL;
   jboolean isc;
   int nxv = (*env)->GetArrayLength(env,g)/(2*ndim);
   if (nx > nxv)
      return;
// Extract array from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pg = (*env)->GetDoubleArrayElements(env,g,&isc);
   if (inorder==LINEAR) {
      vccopy1((double complex*)&pf[0],(double complex*)pg,nx,ndim,nxv);
   }
   else {
      vccopy1((double complex*)&pf[ndim],(double complex*)pg,nx,ndim,
              nxv);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseDoubleArrayElements(env,g,pg,0);
   return;
}
