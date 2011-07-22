//----------------------------------------------------------------------
#include <jni.h>
#include <string.h>
#include "simul1d_diag1d.h"
#include "diag1lib_f.h"
#include "globals.h"

// Java wrappers for 1d PIC Fortran77 libraries diag1lib.f, diag1xlib.f

JNIEXPORT void JNICALL
Java_simul1d_diag1d_writebf(JNIEnv *env, jclass cls, jdoubleArray f,
                            jint nx, jintArray iunit, jintArray nrec,
                            jstring name, jint order) {
   jdouble *pf = NULL;
   jint *piunit = NULL, *pnrec = NULL;
   jboolean isc;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   piunit = (*env)->GetIntArrayElements(env,iunit,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   if (pnrec[0] <= 0) {
      lrec = fgetlrec(pf,nx);
      if (pnrec[0] < 0)
         piunit[0] = fgetunit(piunit[0]);
   }
   if (order==LINEAR) {
      fwrite1(&pf[0],nx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
   else {
      fwrite1(&pf[1],nx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,iunit,piunit,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_readbf(JNIEnv *env, jclass cls, jdoubleArray f,
                           jint nx, jint iunit, jintArray nrec,
                           jintArray ierr, jstring name, jint order) {
   jdouble *pf = NULL;
   jint *pnrec = NULL, *pierr = NULL;
   jboolean isc;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   pierr = (*env)->GetIntArrayElements(env,ierr,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   if (pnrec[0] <= 0) {
      lrec = fgetlrec(pf,nx);
   }
   if (order==LINEAR) {
      fread1(&pf[0],nx,nxv,iunit,(int *)pnrec,lrec,str,
            (int *)&pierr[0]);
   }
   else {
      fread1(&pf[1],nx,nxv,iunit,(int *)pnrec,lrec,str,
             (int *)&pierr[0]);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseIntArrayElements(env,ierr,pierr,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_writebfv(JNIEnv *env, jclass cls, jdoubleArray f,
                             jint nx, jintArray iunit, jintArray nrec,
                             jstring name, jint ndim, jint order) {
   jdouble *pf = NULL;
   jint *piunit = NULL, *pnrec = NULL;
   jboolean isc;
   int nnx;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   piunit = (*env)->GetIntArrayElements(env,iunit,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   nnx = ndim*nx;
   if (pnrec[0] <= 0) {
      lrec = fgetlrec(pf,nnx);
      if (pnrec[0] < 0)
         piunit[0] = fgetunit(piunit[0]);
   }
   if (order==LINEAR) {
      fwrite1(&pf[0],nnx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
   else {
      fwrite1(&pf[1],nnx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,iunit,piunit,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_readbfv(JNIEnv *env, jclass cls, jdoubleArray f,
                            jint nx, jint iunit, jintArray nrec,
                            jintArray ierr, jstring name, jint ndim,
                            jint order) {
   jdouble *pf = NULL;
   jint *pnrec = NULL, *pierr = NULL;
   jboolean isc;
   int nnx;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   pierr = (*env)->GetIntArrayElements(env,ierr,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   nnx = ndim*nx;
   if (pnrec[0] <= 0) {
      lrec = fgetlrec(pf,nnx);
   }
   if (order==LINEAR) {
      fread1(&pf[0],nnx,nxv,iunit,(int *)pnrec,lrec,str,
             (int *)&pierr[0]);
   }
   else {
      fread1(&pf[1],nnx,nxv,iunit,(int *)pnrec,lrec,str,
             (int *)&pierr[0]);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseIntArrayElements(env,ierr,pierr,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_writebfc(JNIEnv *env, jclass cls, jdoubleArray f,
                             jint nx, jintArray iunit, jintArray nrec,
                             jstring name, jint order) {
   jdouble *pf = NULL;
   jint *piunit = NULL, *pnrec = NULL;
   jboolean isc;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f)/2;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   piunit = (*env)->GetIntArrayElements(env,iunit,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   if (pnrec[0] <= 0) {
      lrec = fgetclrec((double complex*)pf,nx);
      if (pnrec[0] < 0)
         piunit[0] = fgetunit(piunit[0]);
   }
   if (order==LINEAR) {
      fcwrite1(&pf[0],nx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
   else {
      fcwrite1(&pf[2],nx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,iunit,piunit,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_readbfc(JNIEnv *env, jclass cls, jdoubleArray f,
                            jint nx, jint iunit, jintArray nrec,
                            jintArray ierr, jstring name, jint order) {
   jdouble *pf = NULL;
   jint *pnrec = NULL, *pierr = NULL;
   jboolean isc;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f)/2;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   pierr = (*env)->GetIntArrayElements(env,ierr,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   if (pnrec[0] <= 0) {
      lrec = fgetclrec((double complex *)pf,nx);
   }
   if (order==LINEAR) {
      fcread1(&pf[0],nx,nxv,iunit,(int *)pnrec,lrec,str,
             (int *)&pierr[0]);
   }
   else {
      fcread1(&pf[2],nx,nxv,iunit,(int *)pnrec,lrec,str,
              (int *)&pierr[0]);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseIntArrayElements(env,ierr,pierr,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_writebfvc(JNIEnv *env, jclass cls, jdoubleArray f,
                              jint nx, jintArray iunit, jintArray nrec,
                              jstring name, jint ndim, jint order) {
   jdouble *pf = NULL;
   jint *piunit = NULL, *pnrec = NULL;
   jboolean isc;
   int nnx;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f)/2;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   piunit = (*env)->GetIntArrayElements(env,iunit,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   nnx = ndim*nx;
   if (pnrec[0] <= 0) {
      lrec = fgetclrec((double complex*)pf,nnx);
      if (pnrec[0] < 0)
         piunit[0] = fgetunit(piunit[0]);
   }
   if (order==LINEAR) {
      fcwrite1(&pf[0],nnx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
   else {
      fcwrite1(&pf[2],nnx,nxv,(int *)piunit,(int *)pnrec,lrec,str);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,iunit,piunit,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_readbfvc(JNIEnv *env, jclass cls, jdoubleArray f,
                             jint nx, jint iunit, jintArray nrec,
                             jintArray ierr, jstring name, jint ndim,
                             jint order) {
   jdouble *pf = NULL;
   jint *pnrec = NULL, *pierr = NULL;
   jboolean isc;
   int nnx;
   int lrec = 0;
   int nxv = (*env)->GetArrayLength(env,f)/2;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   pierr = (*env)->GetIntArrayElements(env,ierr,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
   nnx = ndim*nx;
   if (pnrec[0] <= 0) {
      lrec = fgetclrec((double complex *)pf,nnx);
   }
   if (order==LINEAR) {
      fcread1(&pf[0],nnx,nxv,iunit,(int *)pnrec,lrec,str,
              (int *)&pierr[0]);
   }
   else {
      fcread1(&pf[2],nnx,nxv,iunit,(int *)pnrec,lrec,str,
              (int *)&pierr[0]);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseIntArrayElements(env,ierr,pierr,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_writef(JNIEnv *env, jclass cls, jdoubleArray f,
                           jint iunit, jintArray nrec, jstring name) {
   jdouble *pf = NULL;
   jint *pnrec = NULL;
   jboolean isc;
   int lrec = 0;
   int nxp = (*env)->GetArrayLength(env,f);
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
//
   fwrite0(pf,nxp,iunit,(int *)pnrec,str);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_vdist(JNIEnv *env, jclass cls, jdoubleArray part,
                          jdoubleArray fv, jdoubleArray  fvm, jint np,
                          jint nmv, jint idimp, jint idimv) {
   jdouble *ppart = NULL, *pfv = NULL, *pfvm = NULL;
   jboolean isc;
   int nmvf = (*env)->GetArrayLength(env,fv)/idimv;
// Extract arrays from Java
   ppart = (*env)->GetDoubleArrayElements(env,part,&isc);
   pfv = (*env)->GetDoubleArrayElements(env,fv,&isc);
   pfvm = (*env)->GetDoubleArrayElements(env,fvm,&isc);
   if (idimv==1) {
      vdist1(ppart,pfv,pfvm,idimp,np,nmv,nmvf);
   }
   else if (idimv==3) {
      vdist13(ppart,pfv,pfvm,idimp,np,nmv,nmvf);
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,part,ppart,0);
   (*env)->ReleaseDoubleArrayElements(env,fv,pfv,0);
   (*env)->ReleaseDoubleArrayElements(env,fvm,pfvm,0);
   return;
}

JNIEXPORT jint JNICALL
Java_simul1d_diag1d_fgetunit(JNIEnv *env, jclass cls, jint iunit) {
   return fgetunit(iunit);
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_fbfopen(JNIEnv *env, jclass cls, jdoubleArray f,
                            jint nx, jint iunit, jintArray nrec,
                            jstring fname) {
   jdouble *pf = NULL;
   jint *pnrec = NULL;
   jboolean isc;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,fname,0);
//
   fbfopen(pf,nx,iunit,(int *)pnrec,str);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,fname,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_fbfvopen(JNIEnv *env, jclass cls, jdoubleArray f,
                             jint nx, jint ndim, jint iunit, 
                             jintArray nrec, jstring fname) {
   jdouble *pf = NULL;
   jint *pnrec = NULL;
   jboolean isc;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,fname,0);
//
   fbfvopen(pf,nx,ndim,iunit,(int *)pnrec,str);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,fname,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_fbfcopen(JNIEnv *env, jclass cls, jdoubleArray f,
                             jint nx, jint iunit, jintArray nrec, 
                             jstring fname) {
   jdouble *pf = NULL;
   jint *pnrec = NULL;
   jboolean isc;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,fname,0);
//
   fbfcopen((double complex *)pf,nx,iunit,(int *)pnrec,str);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,fname,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_fbfvcopen(JNIEnv *env, jclass cls, jdoubleArray f,
                              jint nx, jint ndim, jint iunit, 
                              jintArray nrec, jstring fname) {
   jdouble *pf = NULL;
   jint *pnrec = NULL;
   jboolean isc;
// Extract arrays from Java
   pf = (*env)->GetDoubleArrayElements(env,f,&isc);
   pnrec = (*env)->GetIntArrayElements(env,nrec,&isc);
   const char *str = (*env)->GetStringUTFChars(env,fname,0);
//
   fbfvcopen((double complex *)pf,nx,ndim,iunit,(int *)pnrec,str);
//
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,f,pf,0);
   (*env)->ReleaseIntArrayElements(env,nrec,pnrec,0);
   (*env)->ReleaseStringUTFChars(env,fname,str);
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_frclose(JNIEnv *env, jclass cls, jint iunit) {
//
   frclose(iunit);
//
   return;
}

JNIEXPORT void JNICALL
Java_simul1d_diag1d_frwnml1(JNIEnv *env, jclass cls, jstring fname,
                            jstring name, jobjectArray cn,
                            jdoubleArray ddata, jobjectArray sdata,
                            jint iunit, jint isign, jintArray ml,
                            jintArray ms, jintArray ierr) {
   jdouble *pddata = NULL;
   jint *pml = NULL, *pms = NULL, *pierr = NULL;
   jboolean isc;
   jstring js;
   int j, mcn, msd;
   const char *cs;
   char cnull[] = "";
   int nl = (*env)->GetArrayLength(env,ddata);
   int ns = (*env)->GetArrayLength(env,sdata);
// Extract arrays from Java
// Do we need to do something with cn and sdata?????
   pddata = (*env)->GetDoubleArrayElements(env,ddata,&isc);
   pml = (*env)->GetIntArrayElements(env,ml,&isc);
   pms = (*env)->GetIntArrayElements(env,ms,&isc);
   pierr = (*env)->GetIntArrayElements(env,ierr,&isc);
   const char *fstr = (*env)->GetStringUTFChars(env,fname,0);
   const char *str = (*env)->GetStringUTFChars(env,name,0);
// First determine the number of characters in Fortran strings
   frwnml1(fstr,str,cnull,pddata,cnull,0,0,nl,ns,&mcn,&msd,(int *)pml,
           (int *)pms,(int *)pierr);
// Then allocate appropriate C string memory
   char ccn[mcn*(nl+ns)], cdata[msd*ns];
// Copy Java strings to C
   if (isign > 0) {
// Copy variable names, insert blanks for Fortran
      for (j = 0; j < mcn*(nl+ns); j++) {
         ccn[j] = ' ';
      }
      for (j = 0; j < *pml+*pms; j++) {
         js = (*env)->GetObjectArrayElement(env,cn,j);
         cs = (*env)->GetStringUTFChars(env,js,0);
         strcpy(&ccn[j*mcn],cs);
         (*env)->ReleaseStringUTFChars(env,js,cs);
      }
// Copy String variables, insert blanks for Fortran
      for (j = 0; j < msd*ns; j++) {
         cdata[j] = ' ';
      }
      for (j = 0; j < *pms; j++) {
         js = (*env)->GetObjectArrayElement(env,sdata,j);
         cs = (*env)->GetStringUTFChars(env,js,0);
         strcpy(&cdata[j*msd],cs);
         (*env)->ReleaseStringUTFChars(env,js,cs);
      }
   }
// reads or write Fortran namelist
   frwnml1(fstr,str,ccn,pddata,cdata,iunit,isign,nl,ns,&mcn,&msd,
           (int *)pml,(int *)pms,(int *)pierr);
// Copy C strings to Java
   if (isign < 0) {
// Copy variable names
      for (j = 0; j < *pml+*pms; j++) {
         js = (*env)->NewStringUTF(env,&ccn[j*mcn]);
         (*env)->SetObjectArrayElement(env,cn,j,js);
      }
// Copy String variables
      for (j = 0; j < *pms; j++) {
         js = (*env)->NewStringUTF(env,&cdata[j*msd]);
         (*env)->SetObjectArrayElement(env,sdata,j,js);
      }
   }
// Release the arrays back to Java
   (*env)->ReleaseDoubleArrayElements(env,ddata,pddata,0);
   (*env)->ReleaseIntArrayElements(env,ml,pml,0);
   (*env)->ReleaseIntArrayElements(env,ms,pms,0);
   (*env)->ReleaseIntArrayElements(env,ierr,pierr,0);
   (*env)->ReleaseStringUTFChars(env,fname,fstr);
   (*env)->ReleaseStringUTFChars(env,name,str);
   return;
}


