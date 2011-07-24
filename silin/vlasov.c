#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#ifdef _KAPPA
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.141592653589793
#define TPI 6.28318530717958647692528676655900577
#define TPI3 248.050213442398561403810520536811162

//variables for the Adams-Bashforth scheme
#define B1 2.291666666666667
#define B2 -2.458333333333333
#define B3 1.541666666666667
#define B4 -0.375

// declaration of all generic and auxilary variables
char tag[20], header[200];
int nx, ny, nz, nr, i, ir, ix, iv, iy, iz, k, n, il, i0;
int ivx, ivy, ivz, nw, nwy, nwz, iout, nvmax, nvmax2, nvymax, nmax;
int savestep, pdfstep, step, histstep, endstep, mbratio, setup, trigger;
long double EE, EB, Etot;
double xmax, dx, t, dt, psi0, dv, dv3, fact1, alpha, delta, beta, tres;
double *rho, *jx, *jy, *jz, *Ex, *Bx, *By, *Bz, *Ey, *Ez, B0, Bx0, By0, Bz0;
double *gam, *dft, *th, *fluc, *flucx, *flucy, *flucz, mratio, wce, wpe, wcp, wpp;
double *Ax, *Ay, *Az, *Gx, *Gy, *Gz, *d2Ax, *d2Ay, *d2Az, Gxa, Gya, Gza;
double *Axt, *Ayt, *Azt, *Gxt, *Gyt, *Gzt, *Axn, *Ayn, *Azn, *Gxn, *Gyn, *Gzn;
double *Ax0, *Ay0, *Az0, *Gx0, *Gy0, *Gz0, bet, A, B, C, D;
double *in, *out, *kx1, *dGx, *dGy, *dGz, *dAydz, *dAxdz, *dAxdy;
double save_t, hist_t, pdf_t, end_t, E0, Et, coef3, coef4, dt1, dt2, dt3, norm;
double n0, c, c2, lDp, lDe, L, vA, vTp, vTe, vs, dv, dv2, dv3, Pw, Pxtot, ptot, b;
double complex *df2, *fbk, coef, coef1, coef2, fact2, coef5, dummy1, dummy2, dummy3;
double complex *dfdvy, *daxk, *dayk, *dazk, *kx, *bp, *bm, *dfdv, *G, *rhok;
double complex *jxk, *jyk, *jzk, *ky, *kz, *fink, *foutk;
double complex *axk, *ayk, *azk, *ink, *outk, *f1, *f2;

fftw_plan db1, db2, dxf, dxb, dyf, dyb, dzf, dzb, *fdxf, *fdxb; 
fftw_plan transjx, transjy, transjz, transjxk, transjyk, transjzk;
//fftw_plan transaxk,transazk, transayk, transax, transay, transaz;
fftw_plan transrho, transay0, transaz0;
fftw_plan dfxin, dfxout, dfyin, dfyout, dfzin, dfzout, trans1, trans2;

//output file handles, generic
FILE *sEx, *sEy, *sEz, *sBx, *sBy, *sBz, *save, *saves, *dump, *restart, *ref, *streamw, *sinit;
//FILE *smax, *smaxs, *srho, *sjx, *sjy, *sjz, *sAx, *sAy, *sAz; //useful only for tests

//declaration of all particle-specific variables

//the first particle population is defined without any conditions
double complex *F1, *Fn1, *Ft1, *G0, *G1, *G2;
double n1, ud1, M1, T1, Q1, QM1, A1, k1, vmax1, vT1, dV1, dv1, wc1, wp1, Px1;
double *rho1, *jx1, *jy1, *jz1, *Vx1, *Vy1, *Vz1, *u1x, *t1x, *t1y, *t1z;
int nv1, nv1y, nv1z, n11, nv12, nv13, n12, i10; //array dimensions
long double Et1, Ek1, E1;
FILE *sf1, *srho1, *sjx1,  *sjy1, *sjz1;

#if defined(_NONINERT)
double udave, udold, ud, udnew, uxt, *uyt, *uzt, *duy, *duz;
double *uy, *uz, *uyold, *uzold, *uyave, *uzave, *uynew, *uznew;
FILE *suy, *suz;
	#if defined(_AB4)
	double ud0, dud1, dud2, dud3, dud4;
	double *uy0, *uz0, *duy1, *duy2, *duy3, *duy4, *duz1, *duz2, *duz3, *duz4;	
	#endif
#endif

#if defined(_AB4)
//variables for the Adams-Bashforth scheme
double *dAy1, *dAy2, *dAy3, *dAy4, *dAz1, *dAz2, *dAz3, *dAz4;
double complex *F0, *dF1, *dF2, *dF3, *dF4;
	#if defined(_2POP) || defined(_3POP)
	double complex *F20, *dF21, *dF22, *dF23, *dF24;
	#endif
	#if defined(_3POP)
	double complex *F30, *dF31, *dF32, *dF33, *dF34;
	#endif
#endif

#if !defined(_HYBRID)
//if code isn't hybrid, we'll need the electrostatic potential
double *phi, *dphi;
FILE *sphi; 
#if defined(_2D) || defined(_3D)
double *cx, *cy;
#endif
#if defined(_3D)
double *cz;
#endif
#else
// for hybrid codes we'll need these
double *uex, *rhoe, *ubx, *dubx, *pe, *dpe, ge, *pe0, *pet, *pen;
#endif
//the second particle population/species
#if defined(_2POP) || defined(_3POP)
double complex *F2, *Fn2, *Ft2;
int nv2, nv2y, nv2z, n21, nv22, nv23, n22, i20; //array dimensions
double n2, ud2, M2, T2, Q2, QM2, A2, k2, vmax2, vT2, dV2, dv2, wc2, wp2, Px2;
double *rho2, *jx2, *jy2, *jz2, *Vx2, *Vy2, *Vz2, *u2x, *t2x, *t2y, *t2z, *flucy2, *flucz2;
long double Et2, Ek2, E2;
FILE *sf2, *srho2, *sjx2, *sjy2, *sjz2;
#endif
//the third particle population/species
#if defined(_3POP)
double complex *F3, *Fn3, *Ft3;
int nv3, nv3y, nv3z, n31, nv32, nv33, n32, i30; //array dimensions
double n3, ud3, M3, T3, Q3, QM3, A3, k3, vmax3, vT3, dV3, dv3, wc3, wp3, Px3;
double *rho3, *jx3, *jy3, *jz3, *Vx3, *Vy3, *Vz3, *u3x, *t3x, *t3y, *t3z, *flucy3, *flucz3;
long double Et3, Ek3, E3;
FILE *sf3, *srho3, *sjx3, *sjy3, *sjz3;
#endif



#if defined(_TESTELEC)
//this is an unfinished part for test-particle electrons, to calculate particle orbits
int nphi, nth, nev, netot, nep, ip, iphi, ith, ive;
double dth, dphi, dve, dte, mu;
double *xe, *ue, *xeold, *dxe, *vxe, *vxeold, *dvxe, *vye, *vyeold, *dvye, *vze, *vzeold, *dvze, *fe;
nep=nphi*nth*nev; netot=nep*nx;
xe= (double *) calloc(netot, sizeof(double));
vxe= (double *) calloc(netot, sizeof(double));
vye= (double *) calloc(netot, sizeof(double));
vze= (double *) calloc(netot, sizeof(double));
fe= (double *) calloc(netot, sizeof(double));

	for(ix=0;ix<nx;ix++){
		for(ip=0;ip<nep;ip++){xe[ix*nep+ip]=dx*ix;}
		for(ive=0;ive<nev;ive++){
			ue[ix*nep+ive*nth*nphi]=c*sqrt((.195449+ive*.195449)*(.195449+ive*.195449)+2.*(.195449+ive*.195449));}
		for(ive=0;ive<nev;ive++){	
			for(ith=0;ith<nth;ith++)
			for(iphi=0;iphi<nphi;iphi++){i=ix*nep+ive*nth*nphi+ith*nphi+iphi;
			vxe[i]=
			fe[i]=mu/(4.*Pi*c*c*c*k2mu)*exp(-mu*sqrt(1.+
	}


void test_electrons(void){
	for(i=0;i<netot;i++){
		xe[i]+=
	}
}
#endif




double myrand(double limit){
	double random;
	random=limit*rand()/RAND_MAX;
	return random;
}
double Noise(double limit){
	double Grandom;
	Grandom=(myrand(limit)+myrand(limit)+myrand(limit)+myrand(limit)
   	+myrand(limit))/2.5-limit;
	return Grandom;
}
#if defined(_OPENMP)
int runtime(int t0){ //check total runtime since start
	static int tstart;
	double t;
	if(t0==0){tstart=(int)omp_get_wtime(); t=0.;}
	else{t=(int)omp_get_wtime()-tstart;} 
  	return (t);
}
double zeit(void){ //check total runtime for this step, since previous call
 	static double tima=0.;
	double tim, t=0.;
	tim=(double)omp_get_wtime(); t=tim-tima; tima=tim;
  	return (t);
}
#else
int runtime(int t0){ //check total runtime since start
	static int tstart;
	double t;
	if(t0==0){tstart=(int)clock()/CLOCKS_PER_SEC; t=0.;}
	else{t=(int)clock()/CLOCKS_PER_SEC-tstart;} 
	return (t);
}
double zeit(void){ //check total runtime for this step, since previous call
	static double tima=0.;
	double tim, t;
	tim=(double)clock()/CLOCKS_PER_SEC; t=tim-tima; tima=tim;
	return (t);
}
#endif
int sign(int x){ //sign function
	if(x>0) return 1;
	if(x==0) return 0;
	else return -1;
}
void intx1d(double *in, double *out, double fact){
	fftw_execute(dfxin); 
	for(ix=1;ix<nw+1;ix++){outk[ix]=fact*coef1*ink[ix]/ix;}
	for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
	fftw_execute(dfxout); 
}
void inty1d(double *in, double *out, double fact){
	fftw_execute(dfyin); 
	for(iy=1;iy<ny/2;iy++){outk[iy]=fact*coef1*ink[iy]/iy;}
	for(iy=1;iy<ny/2;iy++){outk[ny-iy]=0.;} outk[0]=0.; 
	fftw_execute(dfyout); 
}
void intz1d(double *in, double *out, double fact){
	fftw_execute(dfzin); 
	for(iz=1;iz<nz/2;iz++){outk[iz]=fact*coef1*ink[iz]/iz;}
	for(iz=1;iz<nz/2;iz++){outk[nz-iz]=0.;} outk[0]=0.; 
	fftw_execute(dfzout);  
}
void dfdx1d(double *in, double *out, double fact){
	fftw_execute(dfxin); 
	for(ix=1;ix<nw+1;ix++){outk[ix]=fact*coef2*ix*ink[ix];}
	for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
	fftw_execute(dfxout); 
}
void dfdx(double *fin, double *fout, double fact){
	for(iy=0;iy<ny;iy++)
	for(iz=0;iz<nz;iz++){ir=iy*nz+iz;
		for(ix=0;ix<nx;ix++){in[ix]=fin[ix*ny*nz+ir];} fftw_execute(dfxin); 
		for(ix=1;ix<nw+1;ix++){outk[ix]=fact*coef2*ix*ink[ix];}
		for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
		fftw_execute(dfxout); for(ix=0;ix<nx;ix++){fout[ix*ny*nz+ir]=out[ix];} 
	}
}
void d2fdx(double *fin, double *fout, double fact){
	for(ix=0;ix<nx;ix++){in[ix]=fin[ix];} fftw_execute(dfxin); 
	for(ix=1;ix<nw+1;ix++){outk[ix]=fact*coef4*ix*ix*ink[ix];}
	for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
	fftw_execute(dfxout); for(ix=0;ix<nx;ix++){fout[ix]=out[ix];} 
}
void dfdy(double *fin, double *fout, double fact){
	for(ix=0;ix<nx;ix++)
	for(iz=0;iz<nz;iz++){ir=ix*ny*nz+iz;
		for(iy=0;iy<ny;iy++){in[iy]=fin[iy*nz+ir];} fftw_execute(dfyin); 
		for(iy=1;iy<ny/2;iy++){outk[iy]=fact*coef2*iy*ink[iy];}
		for(iy=1;iy<ny/2;iy++){outk[ny-iy]=0.;} outk[0]=0.; 
		fftw_execute(dfyout); for(iy=0;iy<ny;iy++){fout[iy*nz+ir]=out[iy];} 
	}
}
void dfdz(double *fin, double *fout, double fact){
	for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++){ir=ix*ny*nz+iy*nz;
		for(iz=0;iz<nz;iz++){in[iz]=fin[iz+ir];} fftw_execute(dfzin); 
		for(iz=1;iz<nz/2;iz++){outk[iz]=fact*coef2*iz*ink[iz];}
		for(iz=1;iz<nz/2;iz++){outk[nz-iz]=0.;} outk[0]=0.; 
		fftw_execute(dfzout); for(iz=0;iz<nz;iz++){fout[iz+ir]=out[iz];} 
	}
}
void dFdx(double complex *fin, double complex *fout, int nv, double fact){
	int nvyz, nvxyz; //differentiation of complex array over x
	nvyz=(2*nv-1)*(2*nv-1); nvxyz=nv*nvyz;
	#if defined(_OPENMP)
	#pragma omp parallel for private(ix,ivx,iv,ir)
	#endif
	for(iv=0;iv<nvyz;iv++)
	for(ivx=0;ivx<nv;ivx++){ir=ivx*nvyz+iv; //differentiation of the real part
		for(ix=0;ix<nx;ix++){f1[iv*nx+ix]=creal(fin[ix*nvxyz+ir]);}
		fftw_execute(fdxf[iv]);	
		for(ix=0;ix<nw+1;ix++){f1[iv*nx+ix]=fact*coef2*ix*f2[iv*nx+ix];} 
		for(ix=nw+1;ix<nx;ix++){f1[iv*nx+ix]=0.;} //filter out higher harmonics
		//for(ix=nx-nw;ix<nx;ix++){f1[iv*nx+ix]=fact*coef2*ix*f2[iv*nx+ix];} 
		fftw_execute(fdxb[iv]); 
		for(ix=0;ix<nx;ix++){fout[ix*nvxyz+ir]=f2[iv*nx+ix];}
		//differentiation of the imaginary part
		for(ix=0;ix<nx;ix++){f1[iv*nx+ix]=cimag(fin[ix*nvxyz+ir]);}
		fftw_execute(fdxf[iv]);	
		for(ix=0;ix<nw+1;ix++){f1[iv*nx+ix]=fact*coef2*ix*f2[iv*nx+ix];} 
		for(ix=nw+1;ix<nx;ix++){f1[iv*nx+ix]=0.;} //filter out higher harmonics
		//for(ix=nx-nw;ix<nx;ix++){f1[iv*nx+ix]=fact*coef2*ix*f2[iv*nx+ix];} 
		fftw_execute(fdxb[iv]); 
		for(ix=0;ix<nx;ix++){fout[ix*nvxyz+ir]+=I*f2[iv*nx+ix];}
	}
}
#if defined(_3D)
void Laplace(double *fin, double *fout){  //3D Laplace spectral operator, gives fout(r) = Delta fin(r)
	for(iy=0;iy<ny;iy++)
	for(iz=0;iz<nz;iz++){ir=iy*nz+iz;
		for(ix=0;ix<nx;ix++){in[ix]=fin[ix*ny*nz+ir];} fftw_execute(dfxin); 
		for(ix=1;ix<nw+1;ix++){outk[ix]=coef4*ix*ix*ink[ix];}
		for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
		fftw_execute(dfxout); for(ix=0;ix<nx;ix++){fout[ix*ny*nz+ir]=out[ix];} 
	}
	for(ix=0;ix<nx;ix++)
	for(iz=0;iz<nz;iz++){ir=ix*ny*nz+iz;
		for(iy=0;iy<ny;iy++){in[iy]=fin[iy*nz+ir];} fftw_execute(dfyin); 
		for(iy=1;iy<ny/2;iy++){outk[iy]=coef4*iy*iy*ink[iy];}
		for(iy=1;iy<ny/2;iy++){outk[ny-iy]=0.;} outk[0]=0.; 
		fftw_execute(dfyout); for(iy=0;iy<ny;iy++){fout[iy*nz+ir]+=out[iy];} 
	}
	for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++){ir=ix*ny*nz+iy*nz;
		for(iz=0;iz<nz;iz++){in[iz]=fin[iz+ir];} fftw_execute(dfzin); 
		for(iz=1;iz<nz/2;iz++){outk[iz]=coef4*iz*iz*ink[iz];}
		for(iz=1;iz<nz/2;iz++){outk[nz-iz]=0.;} outk[0]=0.; 
		fftw_execute(dfzout); for(iz=0;iz<nz;iz++){fout[iz+ir]+=out[iz];} 
	}
}
#elif defined(_2D)
void Laplace(double *fin, double *fout){  //2D Laplace spectral operator, gives fout(r) = Delta fin(r)
	for(iy=0;iy<ny;iy++){
		for(ix=0;ix<nx;ix++){in[ix]=fin[ix*ny+iy];} fftw_execute(dfxin); 
		for(ix=1;ix<nw+1;ix++){outk[ix]=coef4*ix*ix*ink[ix];}
		for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
		fftw_execute(dfxout); for(ix=0;ix<nx;ix++){fout[ix*ny+iy]=out[ix];} 
	}
	for(ix=0;ix<nx;ix++){ir=ix*ny;
		for(iy=0;iy<ny;iy++){in[iy]=fin[iy+ir];} fftw_execute(dfyin); 
		for(iy=1;iy<ny/2;iy++){outk[iy]=coef4*iy*iy*ink[iy];}
		for(iy=1;iy<ny/2;iy++){outk[ny-iy]=0.;} outk[0]=0.; 
		fftw_execute(dfyout); for(iy=0;iy<ny;iy++){fout[iy+ir]+=out[iy];} 
	}
}
#elif defined(_1D)
void Laplace(double *fin, double *fout){ //1D Laplace spectral operator, gives fout(x) = d2fin(x)/dx2
	for(ix=0;ix<nx;ix++){in[ix]=fin[ix];} fftw_execute(dfxin); 
	for(ix=1;ix<nw+1;ix++){outk[ix]=coef4*ix*ix*ink[ix];}
	for(ix=nw+1;ix<nx;ix++){outk[ix]=0.;} outk[0]=0.; 
	fftw_execute(dfxout); for(ix=0;ix<nx;ix++){fout[ix]=out[ix];} 
}
#endif
#if defined(_3D)
void Poisson(double *fin, double *fout, double pfact){ //3D Poisson solver, solves Delta fout(r) = pfact*fin(r)

	trans1=fftw_plan_dft_r2c_3d(nx,ny,nz,fin,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_3d(nx,ny,nz,foutk,fout,FFTW_MEASURE);

	fftw_execute(trans1);
	for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++)
   	for(iz=0;iz<nz/2+1;iz++){i=(ix*ny+iy)*(nz/2+1)+iz;
         if((cx[ix]+cy[iy]+cz[iz]-3.)!=0){foutk[i]=pfact*fink[i]/(cx[ix]+cy[iy]+cz[iz]-3.);}
         else{foutk[i]=0.;}
		}
	fftw_execute(trans2);
}
#elif defined(_2D)
void Poisson(double *fin, double *fout, double pfact){ //2D Poisson solver

	trans1=fftw_plan_dft_r2c_2d(nx,ny,fin,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_2d(nx,ny,foutk,fout,FFTW_MEASURE);

	fftw_execute(trans1);
	for(ix=0;ix<nx;ix++)
   	for(iy=0;iy<ny/2+1;iy++){i=ix*(ny/2+1)+iy;
         if((cx[ix]+cy[iy]-2.)!=0){foutk[i]=pfact*fink[i]/(cx[ix]+cy[iy]-2.);}
         else{foutk[i]=0.;}
 		}
	fftw_execute(trans2);
}
#elif defined(_1D)
void Poisson(double *fin, double *fout, double pfact){ //1D Poisson solver
	
	trans1=fftw_plan_dft_r2c_1d(nx,fin,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,fout,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=pfact*fink[i]/(i*i);}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
}
#endif

void dumpfield(FILE *streamw, double *field){ //save the real field data in binary format
   fwrite(field, sizeof(double), nr, streamw); fclose(streamw);
}
void dumpfc(FILE *streamw, fftw_complex *f, int size){ //save the complex PDF in binary format
    int iout; 
	iout=nr/2*size; //3*size; //choose the middle point
	fwrite(&f[iout], sizeof(double complex), size, streamw); fclose(streamw);
}

void dumpfc2(FILE *streamw, fftw_complex *f, int size){ //save complex array of variable size (for testing)
    int i, iout;
	iout=nw*size;
    fprintf(streamw,"\n\n");
    for(i=0;i<size;i++){fprintf(streamw,"%le \n ", creal(f[i+iout]));}
    fprintf(streamw,"\n\n");
    for(i=0;i<size;i++){fprintf(streamw,"%le \n ", cimag(f[i+iout]));}
}
void dumpall(FILE *streamw){ //save the data for restart
	fwrite(&t, sizeof(double), 1, streamw);
	fwrite(&step, sizeof(int), 1, streamw);
	fwrite(&F1[0], sizeof(double complex), n11, streamw); 
	#if defined(_2POP)
	fwrite(&F2[0], sizeof(double complex), n21, streamw); 
	#elif defined(_3POP)
	fwrite(&F2[0], sizeof(double complex), n21, streamw); 
	fwrite(&F3[0], sizeof(double complex), n31, streamw); 
	#elif defined(_4POP)
	fwrite(&F2[0], sizeof(double complex), n21, streamw); 
	fwrite(&F3[0], sizeof(double complex), n31, streamw); 
	fwrite(&F4[0], sizeof(double complex), n41, streamw); 
	#elif defined(_5POP)
	fwrite(&F2[0], sizeof(double complex), n21, streamw); 
	fwrite(&F3[0], sizeof(double complex), n31, streamw); 
	fwrite(&F4[0], sizeof(double complex), n41, streamw); 
	fwrite(&F5[0], sizeof(double complex), n51, streamw); 
	#endif
	fwrite(Ax, sizeof(double), nr, streamw);
	fwrite(Ay, sizeof(double), nr, streamw);
	fwrite(Az, sizeof(double), nr, streamw);
	fwrite(Gx, sizeof(double), nr, streamw);
	fwrite(Gy, sizeof(double), nr, streamw);
	fwrite(Gz, sizeof(double), nr, streamw);
	fwrite(Bx, sizeof(double), nr, streamw);
	fwrite(By, sizeof(double), nr, streamw);
	fwrite(Bz, sizeof(double), nr, streamw);
	fwrite(Ex, sizeof(double), nr, streamw);
	fwrite(Ey, sizeof(double), nr, streamw);
	fwrite(Ez, sizeof(double), nr, streamw);
	#if !defined(_HYBRID)
	fwrite(phi, sizeof(double), nr, streamw);
	#endif
	#if defined(_NONINERT)
	fwrite(uyold, sizeof(double), nr, streamw);
	fwrite(uzold, sizeof(double), nr, streamw);
	fwrite(&udold, sizeof(double), 1, streamw);
	#endif
	fclose(streamw);
}	
void reload(FILE *streamw){ //restart the simulation from a pre-saved data check-point
	fread(&t,sizeof(double), 1, streamw);
	fread(&step, sizeof(int), 1, streamw);
	fread(&F1[0], sizeof(double complex), n11, streamw); 
	#if defined(_2POP)
	fread(&F2[0], sizeof(double complex), n21, streamw); 
	#elif defined(_3POP)
	fread(&F2[0], sizeof(double complex), n21, streamw); 
	fread(&F3[0], sizeof(double complex), n31, streamw); 
	#elif defined(_4POP)
	fread(&F2[0], sizeof(double complex), n21, streamw); 
	fread(&F3[0], sizeof(double complex), n31, streamw); 
	fread(&F4[0], sizeof(double complex), n41, streamw); 
	#elif defined(_5POP)
	fread(&F2[0], sizeof(double complex), n21, streamw); 
	fread(&F3[0], sizeof(double complex), n31, streamw); 
	fread(&F4[0], sizeof(double complex), n41, streamw); 
	fread(&F5[0], sizeof(double complex), n51, streamw); 
	#endif
	fread(Ax, sizeof(double), nr, streamw);
	fread(Ay, sizeof(double), nr, streamw);
	fread(Az, sizeof(double), nr, streamw);
	fread(Gx, sizeof(double), nr, streamw);
	fread(Gy, sizeof(double), nr, streamw);
	fread(Gz, sizeof(double), nr, streamw);
	fread(Bx, sizeof(double), nr, streamw);
	fread(By, sizeof(double), nr, streamw);
	fread(Bz, sizeof(double), nr, streamw);
	fread(Ex, sizeof(double), nr, streamw);
	fread(Ey, sizeof(double), nr, streamw);
	fread(Ez, sizeof(double), nr, streamw);
	#if !defined(_HYBRID)
	fread(phi, sizeof(double), nr, streamw);
	#endif
	#if defined(_NONINERT)
	fread(uyold, sizeof(double), nr, streamw);
	fread(uzold, sizeof(double), nr, streamw);
	fread(&udold,sizeof(double), 1, streamw);
	#endif
}
double max(double *field){ //max of real array
    double help=0.;
    int i;
    for(i=0;i<nr;i++){if(fabs(field[i])>help){help=fabs(field[i]);}}
    return(help);
}
double ave(double *field){ //average absolute value of a real array
    double help=0.;
    int i=0;
	#if defined(_OPENMP)
	#pragma omp parallel for private(i) reduction(+: help)
	#endif
    for(i=0;i<nr;i++){help+=fabs(field[i]);}
    return(help/nr);
}
double aver(double complex *field){ //average absolute value of real part of complex array
    double help=0.;
    int i=0;
	#if defined(_OPENMP)
	#pragma omp parallel for private(i) reduction(+: help)
	#endif
    for(i=0;i<nr;i++){help+=fabs(creal(field[i]));}
    return(help/nr);
}
double avei(double complex *field){ //average of imaginary part of complex array
    double help=0.;
    int i=0;
	#if defined(_OPENMP)
	#pragma omp parallel for private(i) reduction(+: help)
	#endif
    for(i=0;i<nr;i++){help+=fabs(cimag(field[i]));}
    return(help/nr);
}
double maxs(double *field){ //max signed of real array
    double h=0., h1=0.;
    int i;
    for(i=0;i<nr;i++){if(fabs(field[i])>h){h=fabs(field[i]);h1=field[i];}}
    return(h1);
}
double aves(double *field){ //signed average of real array
    double help=0.;
    int i=0;
	#if defined(_OPENMP)
	#pragma omp parallel for private(i) reduction(+: help)
	#endif
    for(i=0;i<nr;i++){help+=field[i];}
    return(help/nr);
}

void history(double *field){ //output absolute and signed averages of fields
	double average, averagesign;
   average=ave(field); averagesign=aves(field);
   fprintf(save,"\t %.12e", average); fprintf(saves,"\t %.12e", averagesign);
}

#if defined(_HYBRID)
void printparam( void ){ //save the main plasma parameters for hybrid codes
	FILE  *plasmastream;
	int iout;

  printf("\n\n Simulation parameters are written to plasma.dat... \n");
  for(iout=0;iout<2;iout++){
    if(iout==1){
    	plasmastream = fopen("plasma.dat","w");
    }
    if(iout==0){plasmastream=stdout;}
    fprintf(plasmastream,"\n wpp \n\t\t\t %1.3e",wpp);
    fprintf(plasmastream,"\n wpe \n\t\t\t %1.3e",wpe);
    fprintf(plasmastream,"\n wcp \n\t\t\t %1.3e",wcp);
    fprintf(plasmastream,"\n wce \n\t\t\t %1.3e",wce);
    fprintf(plasmastream,"\n vTp \n\t\t\t %1.3e",vTp);
#if defined(_2POP) || defined(_3POP)
    fprintf(plasmastream,"\n vT2 \n\t\t\t %1.3e",vT2);
#endif
#if defined(_3POP)
    fprintf(plasmastream,"\n vT3 \n\t\t\t %1.3e",vT3);
#endif
    fprintf(plasmastream,"\n vA  \n\t\t\t %1.3e",vA);
    fprintf(plasmastream,"\n vA/wcp \n\t\t\t %1.3e",L);
    fprintf(plasmastream,"\n lDp \n\t\t\t %1.3e",lDp);
    fprintf(plasmastream,"\n B0 \n\t\t\t %1.3e",B0);
    fprintf(plasmastream,"\n beta \n\t\t\t %1.3e",beta);
    fprintf(plasmastream,"\n dt \n\t\t\t %1.3e",dt);
    fprintf(plasmastream,"\n dx \n\t\t\t %1.3e",dx);
    fprintf(plasmastream,"\n savestep \n\t\t\t %i",savestep);
    fprintf(plasmastream,"\n endstep \n\t\t\t %i",endstep);
    fprintf(plasmastream,"\n histstep \n\t\t\t %i",histstep);
    fprintf(plasmastream,"\n pdfstep \n\t\t\t %i",pdfstep);
	if(iout==1){fclose(plasmastream);}
  }
printf("\n");
}
#else
void printparam( void ){ //save the main plasma parameters for electron and general codes
	FILE  *plasmastream;
	int iout;

  printf("\n\n Simulation parameters are written to plasma.dat... \n");
  for(iout=0;iout<2;iout++){
    if(iout==1){
    	plasmastream = fopen("plasma.dat","w");
    }
    if(iout==0){plasmastream=stdout;}
    fprintf(plasmastream,"\n wpp \n\t\t\t %1.3e",wpp);
    fprintf(plasmastream,"\n wpe \n\t\t\t %1.3e",wpe);
    fprintf(plasmastream,"\n wcp \n\t\t\t %1.3e",wcp);
    fprintf(plasmastream,"\n wce \n\t\t\t %1.3e",wce);
    fprintf(plasmastream,"\n vTp \n\t\t\t %1.3e",vTp);
    fprintf(plasmastream,"\n vTe \n\t\t\t %1.3e",vTe);
    fprintf(plasmastream,"\n vA  \n\t\t\t %1.3e",vA);
    fprintf(plasmastream,"\n vs  \n\t\t\t %1.3e",vs);
    fprintf(plasmastream,"\n lDp \n\t\t\t %1.3e",lDp);
    fprintf(plasmastream,"\n lDe \n\t\t\t %1.3e",lDe);
    fprintf(plasmastream,"\n vA/wcp \n\t\t\t %1.3e",L);
    fprintf(plasmastream,"\n B0 \n\t\t\t %1.3e",B0);
    fprintf(plasmastream,"\n beta \n\t\t\t %1.3e",beta);
    fprintf(plasmastream,"\n dt \n\t\t\t %1.3e",dt);
    fprintf(plasmastream,"\n dx \n\t\t\t %1.3e",dx);
    fprintf(plasmastream,"\n dve (cm/sec)\n\t\t\t %1.3e",dVe);
    fprintf(plasmastream,"\n savestep \n\t\t\t %i",savestep);
    fprintf(plasmastream,"\n endstep \n\t\t\t %i",endstep);
    fprintf(plasmastream,"\n histstep \n\t\t\t %i",histstep);
    fprintf(plasmastream,"\n pdfstep \n\t\t\t %i",pdfstep);
	if(iout==1){fclose(plasmastream);}
  }
printf("\n");
}
#endif

#if defined(_3POP)
void dumpdist(void){ //PDF diagnostics in codes with 3 populations
	#pragma omp parallel sections
   {
	#pragma omp section
   {
    sf1=fopen("f1.dat","a+b"); dumpfc(sf1,F1,nv13); 
    }
	#pragma omp section
   {
    sf2=fopen("f2.dat","a+b"); dumpfc(sf2,F2,nv23); 
   }
	#pragma omp section
   {
    sf3=fopen("f3.dat","a+b"); dumpfc(sf3,F3,nv33); 
   }
   }
}
#elif defined(_2POP)
void dumpdist(void){ //or 2 populations
	#pragma omp parallel sections
   {
	#pragma omp section
   {
    sf1=fopen("f1.dat","a+b"); dumpfc(sf1,F1,nv13); 
    }
	#pragma omp section
   {
    sf2=fopen("f2.dat","a+b"); dumpfc(sf2,F2,nv23); 
   }
   }
}
#elif defined(_1POP)
void dumpdist(void){ //or 1 population
    sf1=fopen("f1.dat","a+b"); dumpfc(sf1,F1,nv13); 
}
#endif

void dumpdata(void){ //save the field data
	#pragma omp parallel sections
   {
#if defined(_2POP) || defined(_3POP)
	#pragma omp section
   {
    srho2=fopen("rho2.dat","a+b"); dumpfield(srho2,rho2);
   }
	#pragma omp section
   {
    sjx2=fopen("jx2.dat","a+b"); dumpfield(sjx2,jx2);
   }
	#pragma omp section
   {
    sjy2=fopen("jy2.dat","a+b"); dumpfield(sjy2,jy2);
   }
	#pragma omp section
   {
    sjz2=fopen("jz2.dat","a+b"); dumpfield(sjz2,jz2);
   }
#endif
#if defined(_3POP)
	#pragma omp section
   {
    srho3=fopen("rho3.dat","a+b"); dumpfield(srho3,rho3);
   }
	#pragma omp section
   {
    sjx3=fopen("jx3.dat","a+b"); dumpfield(sjx3,jx3);
   }
	#pragma omp section
   {
    sjy3=fopen("jy3.dat","a+b"); dumpfield(sjy3,jy3);
   }
	#pragma omp section
   {
    sjz3=fopen("jz3.dat","a+b"); dumpfield(sjz3,jz3);
   }
#endif
	/*#pragma omp section
   {
    srho=fopen("rho.dat","a+b"); dumpfield(srho,rho);
   }
   	#pragma omp section
   {
    sjx=fopen("jx.dat","a+b"); dumpfield(sjx,jx);
   }
	#pragma omp section
   {
    sjy=fopen("jy.dat","a+b"); dumpfield(sjy,jy);
   }
	#pragma omp section
   {
    sjz=fopen("jz.dat","a+b"); dumpfield(sjz,jz);
   }
*/
	#pragma omp section
   {
    srho1=fopen("rho1.dat","a+b"); dumpfield(srho1,rho1);
   }
	#pragma omp section
   {
    sjx1=fopen("jx1.dat","a+b"); dumpfield(sjx1,jx1);
   }
	#pragma omp section
   {
    sjy1=fopen("jy1.dat","a+b"); dumpfield(sjy1,jy1);
   }
	#pragma omp section
   {
    sjz1=fopen("jz1.dat","a+b"); dumpfield(sjz1,jz1);
   }
	#pragma omp section
   {
    sEx=fopen("Ex.dat","a+b"); dumpfield(sEx,Ex);
   }
	#pragma omp section
   {
    sEy=fopen("Ey.dat","a+b"); dumpfield(sEy,Ey);
   }
	#pragma omp section
   {
    sEz=fopen("Ez.dat","a+b"); dumpfield(sEz,Ez);
   }
	#pragma omp section
   {
    sBx=fopen("Bx.dat","a+b"); dumpfield(sBx,Bx);
   }
	#pragma omp section
   {
    sBy=fopen("By.dat","a+b"); dumpfield(sBy,By);
   }
	#pragma omp section
   {
    sBz=fopen("Bz.dat","a+b"); dumpfield(sBz,Bz);
   }
	#if !defined(_HYBRID)
	#pragma omp section
   {
    sphi=fopen("phi.dat","a+b"); dumpfield(sphi,phi);
   }
   #endif
	#if defined(_NONINERT)
   #pragma omp section
   {
    suy=fopen("uy.dat","a+b"); dumpfield(suy,uyave);
   }
	#pragma omp section
   {
    suz=fopen("uz.dat","a+b"); dumpfield(suz,uzave);
   }
   #endif
   }
}
#if defined(_3POP)
void histdata(void){ //save the average simulation data
	save=fopen("timeave.dat","a+b"); fprintf(save,"\n %g ",t);
	saves = fopen("timeavesign.dat","a+b"); fprintf(saves,"\n %g ",t);

	history(rho1); history(rho2); history(rho3); history(rho); history(jx1); history(jx2); history(jx3);
	history(jy1); history(jy2); history(jy3);  history(jz1); history(jz2); history(jz3); history(Ex); 
	history(Ey); history(Ez); history(Bx); history(By); history(Bz); 
	history(t1x); history(t2x); history(t3x); history(t1y); history(t2y);  history(t3y); 
	history(t1z); history(t2z); history(t3z);

	fprintf(save,"\t %.12Le", Et1); fprintf(save,"\t %.12Le", Et2); fprintf(save,"\t %.12Le", Et3);
	fprintf(save,"\t %.12Le", Ek1); fprintf(save,"\t %.12Le", Ek2); fprintf(save,"\t %.12Le", Ek3);
	fprintf(save,"\t %.12Le", EE); fprintf(save,"\t %.12Le", EB); fprintf(save,"\t %.12Le", Etot);
	fprintf(save,"\t %.12Le", Px1); fprintf(save,"\t %.12Le", Px2); fprintf(save,"\t %.12Le", Px3);
	fprintf(save,"\t %.12Le", Pw); fprintf(save,"\t %.12Le", Pxtot); 
	#if defined(_NONINERT)
	fprintf(save,"\t %.12le", udold); history(uyold); history(uzold);
	#endif
	fclose(save); fclose(saves); 
}
#elif defined(_2POP)
void histdata(void){
	save=fopen("timeave.dat","a+b"); fprintf(save,"\n %g ",t);
	saves = fopen("timeavesign.dat","a+b"); fprintf(saves,"\n %g ",t);

	history(rho1); history(rho2); history(rho); history(jx1); history(jx2);
	history(jy1); history(jy2);  history(jz1); history(jz2); history(Ex); 
	history(Ey); history(Ez); history(Bx); history(By); history(Bz); 
	history(t1x); history(t2x); history(t1y); history(t2y); history(t1z); history(t2z);

	fprintf(save,"\t %.12Le", Et1); fprintf(save,"\t %.12Le", Et2);
	fprintf(save,"\t %.12Le", Ek1); fprintf(save,"\t %.12Le", Ek2);
	fprintf(save,"\t %.12Le", EE); fprintf(save,"\t %.12Le", EB); fprintf(save,"\t %.12Le", Etot);
	fprintf(save,"\t %.12Le", Px1); fprintf(save,"\t %.12Le", Px2); fprintf(save,"\t %.12Le", Pw);
	fprintf(save,"\t %.12le", Pxtot); 
	#if defined(_NONINERT)
	fprintf(save,"\t %.12le", udold); history(uyold); history(uzold);
	#endif
	fclose(save); fclose(saves); 
}
#elif defined(_1POP)
void histdata(void){
	save=fopen("timeave.dat","a+b"); fprintf(save,"\n %g ",t);
	saves = fopen("timeavesign.dat","a+b"); fprintf(saves,"\n %g ",t);

	history(rho1); history(rho); history(jx1); 
	history(jy1); history(jz1); history(Ex); 
	history(Ey); history(Ez); history(Bx); history(By); history(Bz); 
	history(t1x); history(t1y); history(t1z); 

	fprintf(save,"\t %.12Le", Et1); fprintf(save,"\t %.12Le", Ek1); 
	fprintf(save,"\t %.12Le", EE); fprintf(save,"\t %.12Le", EB); fprintf(save,"\t %.12Le", Etot);
	fprintf(save,"\t %.12Le", Px1); fprintf(save,"\t %.12Le", Pw); fprintf(save,"\t %.12Le", Pxtot);
	fclose(save); fclose(saves); 
}
#endif
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// PDF MESH REFINEMENT //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void refine(double complex *f, int nv){ //tricubic interpolation of PDF for simulations of strongly heated populations
	int nvy, nvyz, nvxyz;
	
	//nvy=2*nv-1; nv2=nvy*nvz; nv3=nv*nv2;
	nvy=2*nv-1; nvyz=nvy*nvy; nvxyz=nv*nvyz;
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i,i0)
	for(ix=0;ix<nx;ix++){
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++)
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0]=f[i];
		//if(i0==2262712){printf("\n this is place #1, G[%i]=%1.3e+I*%1.3e",i0,creal(G[i0]),cimag(G[i0]));}
	}
	// interpolate in x-direction at vx>0
	for(ivx=1;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++)
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
		//if(i0+nvyz==2262712){printf("\n this is place #2, G[%i]=%1.3e+I*%1.3e",i0+nvyz,creal(G[i0]),cimag(G[i0]));}
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivy=0;ivy<nv;ivy++)
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	//obtain G0, G1, G2 for the y-interpolation boundaries
	
	ivy=0;
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(-1+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G0[i0]=f[i];
	}
	// interpolate in x-direction at vx>0
	for(ivx=1;ivx<(nv-1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(-1+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G0[i0+nvy]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+(-1+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+ivz;
		G0[i0+nvy]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-0)*nvy+nv+(nv-3)/2-ivz])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-0)*nvy+nv+(nv-3)/2-ivz]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	
	ivy=nv;
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G1[i0]=f[i];
	}
	// interpolate in x-direction at vx>0
	for(ivx=1;ivx<(nv-1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G1[i0+nvy]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+ivz;
		G1[i0+nvy]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	
	ivy=nv+1;
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G2[i0]=f[i];
	}
	// interpolate in x-direction at vx>0
	for(ivx=1;ivx<(nv+1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+ivz;
		G2[i0+nvy]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+ivz;
		G2[i0+nvy]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-ivz]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}

	// interpolate in y-direction 
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++)
	for(ivz=0;ivz<nv;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvy]=-.0625*f[i-nvy]+.5625*f[i]+.5625*f[i+nvy]-.0625*f[i+2*nvy];
		//if(i0+nvy==2262712){printf("\n this is place #3, G[%i]=%1.3e+I*%1.3e",i0+nvy,creal(G[i0]),cimag(G[i0]));}
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=1;ivy<nv-2;ivy++)
	for(ivz=0;ivz<nv;ivz++){i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+nvy]=-.0625*G[i0+nvyz-2*nvy]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2*nvy]-.0625*G[i0+nvyz+4*nvy];
		//if(i0+nvyz+nvy==2262712){printf("\n this is place #4, G[%i]=%1.3e+I*%1.3e",i0+nvyz+nvy,creal(G[i0]),cimag(G[i0]));}
	}
	//need values of G at vy=0, nv+1 and nv+2 for all x 
	//boundaries at vy=0 and vy>=nv-2	
	//at vy=0	
	ivy=0; 
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nv*nvy+2*ivx*nvy+ivz;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivz;
		G[i0+nvyz+nvy]=-.0625*G0[i+nvy]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2*nvy]-.0625*G[i0+nvyz+4*nvy];
		//if(i0+nvyz+nvy==2262712){printf("\n this is place #5, G[%i]=%1.3e+I*%1.3e",i0+nvyz+nvy,creal(G[i0]),cimag(G[i0]));}
	}
	ivy=nv-2;
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nv*nvy+2*ivx*nvy+ivz;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+nvy]=-.0625*G[i0+nvyz-2*nvy]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2*nvy]-.0625*G1[i+nvy];
	}
	ivy=nv-1;
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivz=0;ivz<nv;ivz++){i=ix*nv*nvy+2*ivx*nvy+ivz;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+nvy]=-.0625*G[i0+nvyz-2*nvy]+.5625*G[i0+nvyz]+.5625*G1[i+nvy]-.0625*G2[i+nvy];
	}
	

	// interpolate in z-direction 
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++)
	for(ivz=0;ivz<nv-1;ivz++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0]=f[i];
		G[i0+1]=-.0625*f[i-1]+.5625*f[i]+.5625*f[i+1]-.0625*f[i+2];
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++)
	for(ivz=1;ivz<nv-2;ivz++){
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+1]=-.0625*G[i0+nvyz-2]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2]-.0625*G[i0+nvyz+4];
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++)
	for(ivz=1;ivz<nv-2;ivz++){//i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+nvy+1]=-.0625*G[i0+nvyz+nvy-2]+.5625*G[i0+nvyz+nvy]+.5625*G[i0+nvyz+nvy+2]-.0625*G[i0+nvyz+nvy+4];
	}
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++)
	for(ivz=1;ivz<nv-2;ivz++){//i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvy+1]=-.0625*G[i0+nvy-2]+.5625*G[i0+nvy]+.5625*G[i0+nvy+2]-.0625*G[i0+nvy+4];
	}
	//need G at vz=0, nv+1 and vz+2 for all x and y
	//boundaries at vz=0 and vz>=nv-2
	//at vz=0
	
		//prepare G0, G1, G2 for the vz-interpolation at vz=0, nv+1, nv+2
	//ivz=0; shifted by dvx/2 in x-direction
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy-1+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0]=f[i]; //this is not even necessary
	}
	// interpolate in x-direction at vx>0, shifted by 1/2 grid in y
	for(ivx=1;ivx<(nv+1)/2;ivx++) //this is necessary
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy-1+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0+nvy]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+(ivy+(nv-1)/2)*nvy-1+(nv-1)/2;
		i0=ix*nv*nvy+2*ivy;
		G0[i0+nvy]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-1])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+nv+(nv-3)/2-1]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	ivz=0; 
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy; i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy;
		G[i0+nvyz+1]=-.0625*G0[i+nvy]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2]-.0625*G[i0+nvyz+4];
	}
	
	
	
	//ivz=0
	//ivz=0; shifted by dvy/2 in y-direction and by dvx/2 in x-direction
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0+nvy]=G0[i0+nvy]; 
	}
	// interpolate in x-direction at vx>0, shifted by 1/2 grid in y
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0+nvy+1]=-.0625*G0[i0+nvy-2]+.5625*G0[i0+nvy]+.5625*G0[i0+nvy+2]-.0625*G0[i0+4+nvy];
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy;
		G[i0+nvy+nvyz+1]=-.0625*G1[i+nvy+1]+.5625*G[i0+nvyz+nvy]+.5625*G[i0+nvyz+nvy+2]-.0625*G[i0+nvyz+nvy+4];
	}
	
	//prepare G0, G1, G2 for the vz-interpolation at vz=0, nv+1, nv+2
	//ivz=0; shifted by dvy/2 in y-direction
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy-1+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0]=f[i]; //this is not even necessary
	}
	// interpolate in y-direction 
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy-1+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0+nvy]=-.0625*f[i-nvy]+.5625*f[i]+.5625*f[i+nvy]-.0625*f[i+2*nvy];
	}
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy; i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy;
		G[i0+nvy+1]=-.0625*G0[i+nvy]+.5625*G[i0+nvy]+.5625*G[i0+nvy+2]-.0625*G[i0+nvy+4];
	}
		
	//prepare G0, G1, G2 for the vz-interpolation at vz=nv-2
	ivz=nv-2; //shifted by dvx/2 in x-direction ivz=nv-2
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0]=f[i]; //this is not even necessary
	}
	// interpolate in x-direction at vx>0, shifted by 1/2 grid in y
	for(ivx=1;ivx<(nv+1)/2;ivx++) //this is necessary
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G0[i0+nvy]=-.0625*f[i-nvyz]+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	// interpolate in x-direction at vx=0, extrapolate f_i-1 into the vx<0 region
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivy;
		G0[i0+nvy]=-.0625*(creal(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+(nv-1)/2])
			-I*cimag(f[ix*nvxyz+nvyz+(nv+(nv-3)/2-ivy)*nvy+(nv-1)/2]))
			+.5625*f[i]+.5625*f[i+nvyz]-.0625*f[i+2*nvyz];
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy;
		i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+1]=-.0625*G[i0+nvyz-2]+.5625*G[i0+nvyz]+.5625*G[i0+nvyz+2]-.0625*G0[i+nvy];
	}
	
	//prepare G0, G1, G2 for the vz-interpolation at vz=0, nv+1, nv+2
	//ivz=0; shifted by dvx/2 in x-direction ivz=nv-2
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0]=G0[i0+nvy]; 
	}
	// interpolate in y-direction 
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0+nvy+1]=-.0625*G0[i0+nvy-2]+.5625*G0[i0+nvy]+.5625*G0[i0+nvy+2]-.0625*G0[i0+nvy+4];
	}
	for(ivx=0;ivx<(nv-1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy; i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvyz+nvy+1]=-.0625*G[i0+nvyz+nvy-2]+.5625*G[i0+nvyz+nvy]+.5625*G[i0+nvyz+nvy+2]-.0625*G1[i+nvy+1];
	}
	
	//ivz=nv-2; //shifted by dvx/2 in x-direction ivz=nv-2
	for(ivx=0;ivx<(nv+1)/2;ivx++) // record all coarse grid values
	for(ivy=0;ivy<nv;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0]=f[i]; //this is not even necessary
	}
	// interpolate in y-direction 
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nvxyz+ivx*nvyz+(ivy+(nv-1)/2)*nvy+ivz+(nv-1)/2;
		i0=ix*nv*nvy+2*ivx*nvy+2*ivy;
		G1[i0+1]=-.0625*f[i-nvy]+.5625*f[i]+.5625*f[i+nvy]-.0625*f[i+2*nvy];
	}
	for(ivx=0;ivx<(nv+1)/2;ivx++) 
	for(ivy=0;ivy<nv;ivy++){i=ix*nv*nvy+2*ivx*nvy+2*ivy; i0=ix*nvxyz+2*ivx*nvyz+2*ivy*nvy+2*ivz;
		G[i0+nvy+1]=-.0625*G[i0+nvy-2]+.5625*G[i0+nvy]+.5625*G[i0+nvy+2]-.0625*G1[i+1];
	}
	
	// THE BOUNDARY CONDITIONS AT VZ MIN AND VZ MAX ARE NOT FINISHED!!!
	

	}
	#pragma omp parallel for private(i)
	for(i=0;i<nx*nvxyz;i++){f[i]=G[i];}
}
/*
void refinelin(double complex *f, int nv){ //linear interpolation, for testing purposes only!
	int nv2, nv3, nvy, nvy;
	
	nvy=nvz=2*nv-1; nv2=nvy*nvz; nv3=nv*nv2;
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i,i0)
	for(ix=0;ix<nx;ix++){
	for(ivx=0;ivx<(nv-1)/2;ivx++)
	for(ivy=0;ivy<nv-1;ivy++)
	for(ivz=0;ivz<nv-1;ivz++){i=ix*nv3+ivx*nv2+(ivy+(nv-1)/2)*nvz+ivz+(nv-1)/2;
		i0=ix*nv3+2*ivx*nv2+2*ivy*nvz+2*ivz;
		G[i0]=f[i];
		G[i0+1]=.5*(f[i]+f[i+1]);
		G[i0+nvz]=.5*(f[i]+f[i+nvz]);
		G[i0+nv2]=.5*(f[i]+f[i+nv2]);
		G[i0+nvz+1]=.25*(f[i]+f[i+1]+f[i+nvz]+f[i+nvz+1]);
		G[i0+nv2+1]=.25*(f[i]+f[i+1]+f[i+nv2]+f[i+nv2+1]);
		G[i0+nv2+nvz]=.25*(f[i]+f[i+nvz]+f[i+nv2]+f[i+nv2+nvz]);
		G[i0+nv2+nvz+1]=.125*(f[i]+f[i+1]+f[i+nvz]+f[i+nv2]+f[i+nvz+1]+f[i+nvz+nv2]+f[i+nv2+1]+f[i+nv2+nvz]);
	}
	//boundary at vx=vxmax
	for(ivy=0;ivy<nv-1;ivy++)
	for(ivz=0;ivz<nv-1;ivz++){i=ix*nv3+(nv-1)*nv2/2+(ivy+(nv-1)/2)*nvz+ivz+(nv-1)/2;
		i0=ix*nv3+(nv-1)*nv2+2*ivy*nvz+2*ivz;
		G[i0]=f[i];
		G[i0+1]=.5*(f[i]+f[i+1]);
		G[i0+nvz]=.5*(f[i]+f[i+nvz]);
		G[i0+nvz+1]=.25*(f[i]+f[i+1]+f[i+nvz]+f[i+nvz+1]);
	}
	//boundary at vy=vymax
	for(ivx=0;ivx<(nv-1)/2;ivx++)
	for(ivz=0;ivz<nv-1;ivz++){i=ix*nv3+ivx*nv2+(nv-1+(nv-1)/2)*nvz+ivz+(nv-1)/2;
		i0=ix*nv3+2*ivx*nv2+2*(nv-1)*nvz+2*ivz;
		G[i0]=f[i];
		G[i0+1]=.5*(f[i]+f[i+1]);
		G[i0+nv2]=.5*(f[i]+f[i+nv2]);
		G[i0+nv2+1]=.25*(f[i]+f[i+1]+f[i+nv2]+f[i+nv2+1]);
	}
	//boundary at vz=vzmax
	for(ivx=0;ivx<(nv-1)/2;ivx++)
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nv3+ivx*nv2+(ivy+(nv-1)/2)*nvz+nv-1+(nv-1)/2;
		i0=ix*nv3+2*ivx*nv2+2*ivy*nvz+2*(nv-1);
		G[i0]=f[i];
		G[i0+nvz]=.5*(f[i]+f[i+nvz]);
		G[i0+nv2]=.5*(f[i]+f[i+nv2]);
		G[i0+nv2+nvz]=.25*(f[i]+f[i+nvz]+f[i+nv2]+f[i+nv2+nvz]);
	}
	//boundary at vx=vxmax & vy=vymax
	for(ivz=0;ivz<nv-1;ivz++){i=ix*nv3+(nv-1)*nv2/2+(nv-1+(nv-1)/2)*nvz+ivz+(nv-1)/2;
		i0=ix*nv3+(nv-1)*nv2+2*(nv-1)*nvz+2*ivz;
		G[i0]=f[i];
		G[i0+1]=.5*(f[i]+f[i+1]);
	}	
	//boundary at vx=vxmax & vz=vzmax
	for(ivy=0;ivy<nv-1;ivy++){i=ix*nv3+(nv-1)*nv2/2+(ivy+(nv-1)/2)*nvz+nv-1+(nv-1)/2;
		i0=ix*nv3+(nv-1)*nv2+2*ivy*nvz+2*(nv-1);
		G[i0]=f[i];
		G[i0+nvz]=.5*(f[i]+f[i+nvz]);
	}
	//boundary at vy=vymax & vz=vzmax
	for(ivx=0;ivx<(nv-1)/2;ivx++){i=ix*nv3+ivx*nv2+(nv-1+(nv-1)/2)*nvz+nv-1+(nv-1)/2;
		i0=ix*nv3+2*ivx*nv2+2*(nv-1)*nvz+2*(nv-1);
		G[i0]=f[i];
		G[i0+nv2]=.5*(f[i]+f[i+nv2]);
	}
	//boundary at vx=vxmax & vy=vymax & vz=vzmax
	i=ix*nv3+(nv-1)*nv2/2+(nv-1+(nv-1)/2)*nvz+nv-1+(nv-1)/2;
	i0=ix*nv3+(nv-1)*nv2+2*(nv-1)*nvz+2*(nv-1);
	G[i0]=f[i];
	}
	#pragma omp parallel for private(i)
	for(i=0;i<nx*nv3;i++){f[i]=G[i];}
}*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////// INTEGRATION OF MOMENTS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if defined(_3POP)
void moments(double complex *f1, double complex *f2, double complex *f3){
	dv1=Q1*TPI3/(60.*dV1); dv2=Q2*TPI3/(60.*dV2);  dv3=Q3*TPI3/(60.*dV3);  	
	//#pragma omp parallel for private(ir,ie0)
	#if defined(_NONINERT)
	ud=0.;
	#endif
	for(ir=0;ir<nr;ir++){ 
		i10=ir*nv13+(nv1y+1)*(nv1-1); i20=ir*nv23+(nv2y+1)*(nv2-1); i30=ir*nv33+(nv3y+1)*(nv3-1); 
		rho1[ir]=Q1*TPI3*creal(f1[i10]); rho2[ir]=Q2*TPI3*creal(f2[i20]); rho3[ir]=Q3*TPI3*creal(f3[i30]);

		jx1[ir]=2.*dv1*(cimag(45.*f1[i10+nv12]-9.*f1[i10+2*nv12]+f1[i10+3*nv12]));
		jy1[ir]=dv1*(cimag(f1[i10+3*nv1z]-9.*f1[i10+2*nv1z]+45.*f1[i10+nv1z]
		-45.*f1[i10-nv1z]+9.*f1[i10-2*nv1z]-f1[i10-3*nv1z]));
		jz1[ir]=dv1*(cimag(f1[i10+3]-9.*f1[i10+2]+45.*f1[i10+1]
		-45.*f1[i10-1]+9.*f1[i10-2]-f1[i10-3]));
		
		jx2[ir]=2.*dv2*(cimag(45.*f2[i20+nv22]-9.*f2[i20+2*nv22]+f2[i20+3*nv22]));
		jy2[ir]=dv2*(cimag(f2[i20+3*nv2z]-9.*f2[i20+2*nv2z]+45.*f2[i20+nv2z]
		-45.*f2[i20-nv2z]+9.*f2[i20-2*nv2z]-f2[i20-3*nv2z]));
		jz2[ir]=dv2*(cimag(f2[i20+3]-9.*f2[i20+2]+45.*f2[i20+1]
		-45.*f2[i20-1]+9.*f2[i20-2]-f2[i20-3]));
		
		#if defined(_NONINERT)
		jx3[ir]=rho3[ir]*udold+2.*dv3*(cimag(45.*f3[i30+nv32]-9.*f3[i30+2*nv32]+f3[i30+3*nv32]));
		jy3[ir]=rho3[ir]*uyold[ir]+dv3*(cimag(f3[i30+3*nv3z]-9.*f3[i30+2*nv3z]+45.*f3[i30+nv3z]
		-45.*f3[i30-nv3z]+9.*f3[i30-2*nv3z]-f3[i30-3*nv3z]));
		jz3[ir]=rho3[ir]*uzold[ir]+dv3*(cimag(f3[i30+3]-9.*f3[i30+2]+45.*f3[i30+1]
		-45.*f3[i30-1]+9.*f3[i30-2]-f3[i30-3]));
		ud+=jx3[ir]/rho3[ir]; uy[ir]=jy3[ir]/rho3[ir]; uz[ir]=jz3[ir]/rho3[ir];
		#else
		jx3[ir]=2.*dv3*(cimag(45.*f3[i30+nv32]-9.*f3[i30+2*nv32]+f3[i30+3*nv32]));
		jy3[ir]=dv3*(cimag(f3[i30+3*nv3z]-9.*f3[i30+2*nv3z]+45.*f3[i30+nv3z]
		-45.*f3[i30-nv3z]+9.*f3[i30-2*nv3z]-f3[i30-3*nv3z]));
		jz3[ir]=dv3*(cimag(f3[i30+3]-9.*f3[i30+2]+45.*f3[i30+1]
		-45.*f3[i30-1]+9.*f3[i30-2]-f3[i30-3]));
		#endif
		
		#if defined(_HYBRID)
		rho[ir]=0.; rhoe[ir]=rho1[ir]+rho2[ir]+rho3[ir];
		uex[ir]=(jx1[ir]+jx2[ir]+jx3[ir])/rhoe[ir]; 
		ubx[ir]=(jx2[ir]+jx3[ir])/rhoe[ir]; 
		jx[ir]=jx1[ir]+jx2[ir]+jx3[ir]-rhoe[ir]*uex[ir]; 
		#elif defined(_ELEC)
		rho[ir]=n1+n2+n3+rho1[ir]+rho2[ir]+rho3[ir]; jx[ir]=jx1[ir]+jx2[ir]+jx3[ir]; 		
		#else
		rho[ir]=rho1[ir]+rho2[ir]+rho3[ir]; jx[ir]=jx1[ir]+jx2[ir]+jx3[ir]; 
		#endif
		jy[ir]=jy1[ir]+jy2[ir]+jy3[ir]; jz[ir]=jz1[ir]+jz2[ir]+jz3[ir];
	}
	#if defined(_NONINERT)
	ud=ud/nr;
	#endif
#ifdef _1D	
	//filter out high harmonics from all particle moments
	trans1=fftw_plan_dft_r2c_1d(nx,rho,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rho,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jy,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jy,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jz,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);

	#if defined(_HYBRID)
	trans1=fftw_plan_dft_r2c_1d(nx,rhoe,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rhoe,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,uex,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uex,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,ubx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,ubx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
	#if defined(_NONINERT)
	trans1=fftw_plan_dft_r2c_1d(nx,uy,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uy,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,uz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uz,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
#endif
}
#elif defined(_2POP)
void moments(double complex *f1, double complex *f2){
	dv1=Q1*TPI3/(60.*dV1); dv2=Q2*TPI3/(60.*dV2);    	
	#if defined(_NONINERT)
	ud=0.;
	#endif
	for(ir=0;ir<nr;ir++){ 
		i10=ir*nv13+(nv1y+1)*(nv1-1); i20=ir*nv23+(nv2y+1)*(nv2-1);  
		rho1[ir]=Q1*TPI3*creal(f1[i10]); rho2[ir]=Q2*TPI3*creal(f2[i20]); 
		
		jx1[ir]=2.*dv1*(cimag(45.*f1[i10+nv12]-9.*f1[i10+2*nv12]+f1[i10+3*nv12]));
		jy1[ir]=dv1*(cimag(f1[i10+3*nv1z]-9.*f1[i10+2*nv1z]+45.*f1[i10+nv1z]
		-45.*f1[i10-nv1z]+9.*f1[i10-2*nv1z]-f1[i10-3*nv1z]));
		jz1[ir]=dv1*(cimag(f1[i10+3]-9.*f1[i10+2]+45.*f1[i10+1]
		-45.*f1[i10-1]+9.*f1[i10-2]-f1[i10-3]));
		
		
		#if defined(_NONINERT)
		jx2[ir]=rho2[ir]*udold+2.*dv2*(cimag(45.*f2[i20+nv22]-9.*f2[i20+2*nv22]+f2[i20+3*nv22]));
		jy2[ir]=rho2[ir]*uyold[ir]+dv2*(cimag(f2[i20+3*nv2z]-9.*f2[i20+2*nv2z]+45.*f2[i20+nv2z]
		-45.*f2[i20-nv2z]+9.*f2[i20-2*nv2z]-f2[i20-3*nv2z]));
		jz2[ir]=rho2[ir]*uzold[ir]+dv2*(cimag(f2[i20+3]-9.*f2[i20+2]+45.*f2[i20+1]
		-45.*f2[i20-1]+9.*f2[i20-2]-f2[i20-3]));
		ud+=jx2[ir]/rho2[ir]; uy[ir]=jy2[ir]/rho2[ir]; uz[ir]=jz2[ir]/rho2[ir];
		#else
		jx2[ir]=2.*dv2*(cimag(45.*f2[i20+nv22]-9.*f2[i20+2*nv22]+f2[i20+3*nv22]));
		jy2[ir]=dv2*(cimag(f2[i20+3*nv2z]-9.*f2[i20+2*nv2z]+45.*f2[i20+nv2z]
		-45.*f2[i20-nv2z]+9.*f2[i20-2*nv2z]-f2[i20-3*nv2z]));
		jz2[ir]=dv2*(cimag(f2[i20+3]-9.*f2[i20+2]+45.*f2[i20+1]
		-45.*f2[i20-1]+9.*f2[i20-2]-f2[i20-3]));
		#endif
		
		#if defined(_HYBRID)
		rho[ir]=0.; rhoe[ir]=rho1[ir]+rho2[ir];
		uex[ir]=(jx1[ir]+jx2[ir])/rhoe[ir]; ubx[ir]=(jx2[ir])/rhoe[ir]; 
		jx[ir]=jx1[ir]+jx2[ir]-rhoe[ir]*uex[ir]; //total parallel current should be 0
		#elif defined(_ELEC)
		rho[ir]=n1+n2+rho1[ir]+rho2[ir]; jx[ir]=jx1[ir]+jx2[ir]; //for electron codes rho1<0
		#else
		rho[ir]=rho1[ir]+rho2[ir]; jx[ir]=jx1[ir]+jx2[ir]; 
		#endif
		jy[ir]=jy1[ir]+jy2[ir]; jz[ir]=jz1[ir]+jz2[ir];
	}
	#if defined(_NONINERT)
	ud=ud/nr;
	#endif
#ifdef _1D	
	trans1=fftw_plan_dft_r2c_1d(nx,rho,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rho,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jy,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jy,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jz,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	#if defined(_HYBRID)
	trans1=fftw_plan_dft_r2c_1d(nx,rhoe,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rhoe,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,uex,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uex,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,ubx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,ubx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
	#if defined(_NONINERT)
	trans1=fftw_plan_dft_r2c_1d(nx,uy,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uy,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw/2+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw/2+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,uz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uz,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw/2+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw/2+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
#endif
}
#elif defined(_1POP)
void moments(double complex *f1){
	dv1=Q1*TPI3/(60.*dV1);   	
	for(ir=0;ir<nr;ir++){ i10=ir*nv13+(nv1y+1)*(nv1-1);  
		rho1[ir]=Q1*TPI3*creal(f1[i10]); 

		jx1[ir]=2.*dv1*(cimag(45.*f1[i10+nv12]-9.*f1[i10+2*nv12]+f1[i10+3*nv12]));
		jy1[ir]=dv1*(cimag(f1[i10+3*nv1z]-9.*f1[i10+2*nv1z]+45.*f1[i10+nv1z]
		-45.*f1[i10-nv1z]+9.*f1[i10-2*nv1z]-f1[i10-3*nv1z]));
		jz1[ir]=dv1*(cimag(f1[i10+3]-9.*f1[i10+2]+45.*f1[i10+1]
		-45.*f1[i10-1]+9.*f1[i10-2]-f1[i10-3]));
		//1POP codes can only be either hybrid or electron!
		#if defined(_HYBRID)
		rho[ir]=0.; rhoe[ir]=rho1[ir];
		uex[ir]=jx1[ir]/rhoe[ir]; ubx[ir]=0.;
		jx[ir]=jx1[ir]-rhoe[ir]*uex[ir]; //return electron current
		jy[ir]=jy1[ir]; jz[ir]=jz1[ir];
		#else
		rho[ir]=n1+rho1[ir]; jx[ir]=jx1[ir]; //for electron codes rho1<0
		jy[ir]=jy1[ir]; jz[ir]=jz1[ir];		
		#endif
	}
#ifdef _1D	
	for(i=nw+1;i<nx;i++){foutk[i]=0.;} //filter out high harmonics
	trans1=fftw_plan_dft_r2c_1d(nx,rho,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rho,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jy,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jy,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,jz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,jz,FFTW_MEASURE);

	fftw_execute(trans1);
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#if defined(_HYBRID)
	trans1=fftw_plan_dft_r2c_1d(nx,rhoe,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,rhoe,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,uex,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,uex,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,ubx,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,ubx,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////////
///////////////////// INTEGRATION OF 2-nd ORDER MOMENTS ////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if defined(_3POP)
void moments2(void){
	double vx, vy, vz;
	E1=E2=E3=Ek1=Ek2=Ek3=Et1=Et2=Et3=Px1=Px2=Px3=0.;
	dv1=Q1*TPI3*M1/(12.*dV1*dV1); dv2=Q2*TPI3*M2/(12.*dV2*dV2); dv3=Q3*TPI3*M3/(12.*dV3*dV3); 
	//#pragma omp parallel for private(ir,vx,vy,vz,ie0,Ee,Ethe,Ethi)
	for(ir=0;ir<nr;ir++){
	i10=ir*nv13+(nv1y+1)*(nv1-1); i20=ir*nv23+(nv2y+1)*(nv2-1); i30=ir*nv33+(nv3y+1)*(nv3-1);
	
	vx=jx1[ir]/rho1[ir]; vy=jy1[ir]/rho1[ir]; vz=jz1[ir]/rho1[ir];  
	
    t1x[ir]=-2.*dv1*creal(-15.*F1[i10]+16.*F1[i10+nv12]/cos(Vx1[1]*vx)
		-F1[i10+2*nv12]/cos(Vx1[2]*vx))/rho1[ir];
	
    t1y[ir]=-dv1*creal(-F1[i10-2*nv1z]/cos(Vy1[nv1-3]*vy)
	+16.*F1[i10-nv1z]/cos(Vy1[nv1-2]*vy)-30.*F1[i10]
	+16.*F1[i10+nv1z]/cos(Vy1[nv1]*vy)-F1[i10+2*nv1z]/cos(Vy1[nv1+1]*vy))/rho1[ir];
	
    t1z[ir]=-dv1*creal(-F1[i10-2]/cos(Vz1[nv1-3]*vz)
	+16.*F1[i10-1]/cos(Vz1[nv1-2]*vz)-30.*F1[i10]
	+16.*F1[i10+1]/cos(Vz1[nv1]*vz)-F1[i10+2]/cos(Vz1[nv1+1]*vz))/rho1[ir];

    E1+=-.5*dv1*(2.*creal(-15.*F1[i10]+16.*F1[i10+nv12]-F1[i10+2*nv12])+
	creal(-F1[i10-2*nv1z]+16.*F1[i10-nv1z]-30.*F1[i10]+16.*F1[i10+nv1z]-F1[i10+2*nv1z])+
	creal(-F1[i10-2]+16.*F1[i10-1]-30.*F1[i10]+16.*F1[i10+1]-F1[i10+2]))/Q1;

	vx=jx2[ir]/rho2[ir]; vy=jy2[ir]/rho2[ir]; vz=jz2[ir]/rho2[ir];  
	
    t2x[ir]=-2.*dv2*creal(-15.*F2[i20]+16.*F2[i20+nv22]/cos(Vx2[1]*vx)
		-F2[i20+2*nv22]/cos(Vx2[2]*vx))/rho2[ir];
	
    t2y[ir]=-dv2*creal(-F2[i20-2*nv2z]/cos(Vy2[nv2-3]*vy)
	+16.*F2[i20-nv2z]/cos(Vy2[nv2-2]*vy)-30.*F2[i20]
	+16.*F2[i20+nv2z]/cos(Vy2[nv2]*vy)-F2[i20+2*nv2z]/cos(Vy2[nv2+1]*vy))/rho2[ir];
	
    t2z[ir]=-dv2*creal(-F2[i20-2]/cos(Vz2[nv2-3]*vz)
	+16.*F2[i20-1]/cos(Vz2[nv2-2]*vz)-30.*F2[i20]
	+16.*F2[i20+1]/cos(Vz2[nv2]*vz)-F2[i20+2]/cos(Vz2[nv2+1]*vz))/rho2[ir];

    E2+=-.5*dv2*(2.*creal(-15.*F2[i20]+16.*F2[i20+nv22]-F2[i20+2*nv22])+
	creal(-F2[i20-2*nv2z]+16.*F2[i20-nv2z]-30.*F2[i20]+16.*F2[i20+nv2z]-F2[i20+2*nv2z])+
	creal(-F2[i20-2]+16.*F2[i20-1]-30.*F2[i20]+16.*F2[i20+1]-F2[i20+2]))/Q2;

#if defined(_NONINERT)
	vx=jx3[ir]/rho3[ir]-udold; vy=jy3[ir]/rho3[ir]-uyold[ir]; vz=jz3[ir]/rho3[ir]-uzold[ir];  

    t3x[ir]=-2.*dv3*creal(-15.*F3[i30]+16.*F3[i30+nv32]/cos(Vx3[1]*vx)
		-F3[i30+2*nv32]/cos(Vx3[2]*vx))/rho3[ir];
	
    t3y[ir]=-dv3*creal(-F3[i30-2*nv3z]/cos(Vy3[nv3-3]*vy)
	+16.*F3[i30-nv3z]/cos(Vy3[nv3-2]*vy)-30.*F3[i30]
	+16.*F3[i30+nv3z]/cos(Vy3[nv3]*vy)-F3[i30+2*nv3z]/cos(Vy3[nv3+1]*vy))/rho3[ir];
	
    t3z[ir]=-dv3*creal(-F3[i30-2]/cos(Vz3[nv3-3]*vz)
	+16.*F3[i30-1]/cos(Vz3[nv3-2]*vz)-30.*F3[i30]
	+16.*F3[i30+1]/cos(Vz3[nv3]*vz)-F3[i30+2]/cos(Vz3[nv3+1]*vz))/rho3[ir];

    E3+=-.5*dv3*(2.*creal(-15.*F3[i30]+16.*F3[i30+nv32]-F3[i30+2*nv32])+
	creal(-F3[i30-2*nv3z]+16.*F3[i30-nv3z]-30.*F3[i30]+16.*F3[i30+nv3z]-F3[i30+2*nv3z])+
	creal(-F3[i30-2]+16.*F3[i30-1]-30.*F3[i30]+16.*F3[i30+1]-F3[i30+2]))/Q3
	+(udold*jx3[ir]-.5*udold*udold*rho3[ir]
	+uyold[ir]*jy3[ir]-.5*uyold[ir]*uyold[ir]*rho3[ir]
	+uzold[ir]*jz3[ir]-.5*uzold[ir]*uzold[ir]*rho3[ir]
	)/QM3;
#else
	vx=jx3[ir]/rho3[ir]; vy=jy3[ir]/rho3[ir]; vz=jz3[ir]/rho3[ir];  

    t3x[ir]=-2.*dv3*creal(-15.*F3[i30]+16.*F3[i30+nv32]/cos(Vx3[1]*vx)
		-F3[i30+2*nv32]/cos(Vx3[2]*vx))/rho3[ir];
	
    t3y[ir]=-dv3*creal(-F3[i30-2*nv3z]/cos(Vy3[nv3-3]*vy)
	+16.*F3[i30-nv3z]/cos(Vy3[nv3-2]*vy)-30.*F3[i30]
	+16.*F3[i30+nv3z]/cos(Vy3[nv3]*vy)-F3[i30+2*nv3z]/cos(Vy3[nv3+1]*vy))/rho3[ir];
	
    t3z[ir]=-dv3*creal(-F3[i30-2]/cos(Vz3[nv3-3]*vz)
	+16.*F3[i30-1]/cos(Vz3[nv3-2]*vz)-30.*F3[i30]
	+16.*F3[i30+1]/cos(Vz3[nv3]*vz)-F3[i30+2]/cos(Vz3[nv3+1]*vz))/rho3[ir];

    E3+=-.5*dv3*(2.*creal(-15.*F3[i30]+16.*F3[i30+nv32]-F3[i30+2*nv32])+
	creal(-F3[i30-2*nv3z]+16.*F3[i30-nv3z]-30.*F3[i30]+16.*F3[i30+nv3z]-F3[i30+2*nv3z])+
	creal(-F3[i30-2]+16.*F3[i30-1]-30.*F3[i30]+16.*F3[i30+1]-F3[i30+2]))/Q3;
#endif

	Et1+=(t1x[ir]+t1y[ir]+t1z[ir])*rho1[ir]; Px1+=jx1[i];
	Et2+=(t2x[ir]+t2y[ir]+t2z[ir])*rho2[ir]; Px2+=jx2[i];
	Et3+=(t3x[ir]+t3y[ir]+t3z[ir])*rho3[ir]; Px3+=jx3[i];
	}
	Et1=.5*Et1/Q1; Et2=.5*Et2/Q2; Et3=.5*Et3/Q3; Ek1=E1-Et1; Ek2=E2-Et2; Ek3=E3-Et3; 
	Etot=E1+E2+E3+EE+EB; Px1=Px1/QM1; Px2=Px2/QM2; Px3=Px3/QM3; Pxtot=Px1+Px2+Px3;
}
#elif defined(_2POP)
void moments2(void){
	double vx, vy, vz;
	E1=E2=Ek1=Ek2=Et1=Et2=Px1=Px2=0.;
	dv1=Q1*TPI3*M1/(12.*dV1*dV1); dv2=Q2*TPI3*M2/(12.*dV2*dV2);  
	//#pragma omp parallel for private(ir,vx,vy,vz,ie0,Ee,Ethe,Ethi)
	for(ir=0;ir<nr;ir++){
	i10=ir*nv13+(nv1y+1)*(nv1-1); i20=ir*nv23+(nv2y+1)*(nv2-1); 
	
	vx=jx1[ir]/rho1[ir]; vy=jy1[ir]/rho1[ir]; vz=jz1[ir]/rho1[ir];  
	
    t1x[ir]=-2.*dv1*creal(-15.*F1[i10]+16.*F1[i10+nv12]/cos(Vx1[1]*vx)
		-F1[i10+2*nv12]/cos(Vx1[2]*vx))/rho1[ir];
	
    t1y[ir]=-dv1*creal(-F1[i10-2*nv1z]/cos(Vy1[nv1-3]*vy)
	+16.*F1[i10-nv1z]/cos(Vy1[nv1-2]*vy)-30.*F1[i10]
	+16.*F1[i10+nv1z]/cos(Vy1[nv1]*vy)-F1[i10+2*nv1z]/cos(Vy1[nv1+1]*vy))/rho1[ir];
	
    t1z[ir]=-dv1*creal(-F1[i10-2]/cos(Vz1[nv1-3]*vz)
	+16.*F1[i10-1]/cos(Vz1[nv1-2]*vz)-30.*F1[i10]
	+16.*F1[i10+1]/cos(Vz1[nv1]*vz)-F1[i10+2]/cos(Vz1[nv1+1]*vz))/rho1[ir];

    E1+=-.5*dv1*(2.*creal(-15.*F1[i10]+16.*F1[i10+nv12]-F1[i10+2*nv12])+
	creal(-F1[i10-2*nv1z]+16.*F1[i10-nv1z]-30.*F1[i10]+16.*F1[i10+nv1z]-F1[i10+2*nv1z])+
	creal(-F1[i10-2]+16.*F1[i10-1]-30.*F1[i10]+16.*F1[i10+1]-F1[i10+2]))/Q1;

#if defined(_NONINERT)
	vx=jx2[ir]/rho2[ir]-udold; vy=jy2[ir]/rho2[ir]-uyold[ir]; vz=jz2[ir]/rho2[ir]-uzold[ir];  

    t2x[ir]=-2.*dv2*creal(-15.*F2[i20]+16.*F2[i20+nv22]/cos(Vx2[1]*vx)
		-F2[i20+2*nv22]/cos(Vx2[2]*vx))/rho2[ir];
	
    t2y[ir]=-dv2*creal(-F2[i20-2*nv2z]/cos(Vy2[nv2-3]*vy)
	+16.*F2[i20-nv2z]/cos(Vy2[nv2-2]*vy)-30.*F2[i20]
	+16.*F2[i20+nv2z]/cos(Vy2[nv2]*vy)-F2[i20+2*nv2z]/cos(Vy2[nv2+1]*vy))/rho2[ir];
	
    t2z[ir]=-dv2*creal(-F2[i20-2]/cos(Vz2[nv2-3]*vz)
	+16.*F2[i20-1]/cos(Vz2[nv2-2]*vz)-30.*F2[i20]
	+16.*F2[i20+1]/cos(Vz2[nv2]*vz)-F2[i20+2]/cos(Vz2[nv2+1]*vz))/rho2[ir];

    E2+=-.5*dv2*(2.*creal(-15.*F2[i20]+16.*F2[i20+nv22]-F2[i20+2*nv22])+
	creal(-F2[i20-2*nv2z]+16.*F2[i20-nv2z]-30.*F2[i20]+16.*F2[i20+nv2z]-F2[i20+2*nv2z])+
	creal(-F2[i20-2]+16.*F2[i20-1]-30.*F2[i20]+16.*F2[i20+1]-F2[i20+2]))/Q2
	+(udold*jx2[ir]-.5*udold*udold*rho2[ir]
	+uyold[ir]*jy2[ir]-.5*uyold[ir]*uyold[ir]*rho2[ir]
	+uzold[ir]*jz2[ir]-.5*uzold[ir]*uzold[ir]*rho2[ir]
	)/QM2;
#else
	vx=jx2[ir]/rho2[ir]; vy=jy2[ir]/rho2[ir]; vz=jz2[ir]/rho2[ir];  
	
    t2x[ir]=-2.*dv2*creal(-15.*F2[i20]+16.*F2[i20+nv22]/cos(Vx2[1]*vx)
		-F2[i20+2*nv22]/cos(Vx2[2]*vx))/rho2[ir];
	
    t2y[ir]=-dv2*creal(-F2[i20-2*nv2z]/cos(Vy2[nv2-3]*vy)
	+16.*F2[i20-nv2z]/cos(Vy2[nv2-2]*vy)-30.*F2[i20]
	+16.*F2[i20+nv2z]/cos(Vy2[nv2]*vy)-F2[i20+2*nv2z]/cos(Vy2[nv2+1]*vy))/rho2[ir];
	
    t2z[ir]=-dv2*creal(-F2[i20-2]/cos(Vz2[nv2-3]*vz)
	+16.*F2[i20-1]/cos(Vz2[nv2-2]*vz)-30.*F2[i20]
	+16.*F2[i20+1]/cos(Vz2[nv2]*vz)-F2[i20+2]/cos(Vz2[nv2+1]*vz))/rho2[ir];

    E2+=-.5*dv2*(2.*creal(-15.*F2[i20]+16.*F2[i20+nv22]-F2[i20+2*nv22])+
	creal(-F2[i20-2*nv2z]+16.*F2[i20-nv2z]-30.*F2[i20]+16.*F2[i20+nv2z]-F2[i20+2*nv2z])+
	creal(-F2[i20-2]+16.*F2[i20-1]-30.*F2[i20]+16.*F2[i20+1]-F2[i20+2]))/Q2;
#endif
	Et1+=(t1x[ir]+t1y[ir]+t1z[ir])*rho1[ir]; Px1+=jx1[i];
	Et2+=(t2x[ir]+t2y[ir]+t2z[ir])*rho2[ir]; Px2+=jx2[i];
	}
	Et1=.5*Et1/Q1; Et2=.5*Et2/Q2; Ek1=E1-Et1; Ek2=E2-Et2; Etot=E1+E2+EE+EB;
	Px1=Px1/QM1; Px2=Px2/QM2; Pxtot=Px1+Px2;
}
#elif defined(_1POP)
void moments2(void){
	double vx, vy, vz;
	E1=Ek1=Et1=Px1=0.; dv1=Q1*TPI3*M1/(12.*dV1*dV1); 
	//#pragma omp parallel for private(ir,vx,vy,vz,ie0,Ee,Ethe,Ethi)
	for(ir=0;ir<nr;ir++){i10=ir*nv13+(nv1y+1)*(nv1-1);  
	
	vx=jx1[ir]/rho1[ir]; vy=jy1[ir]/rho1[ir]; vz=jz1[ir]/rho1[ir];  
	
    t1x[ir]=-2.*dv1*creal(-15.*F1[i10]+16.*F1[i10+nv12]/cos(Vx1[1]*vx)
		-F1[i10+2*nv12]/cos(Vx1[2]*vx))/rho1[ir];
	
    t1y[ir]=-dv1*creal(-F1[i10-2*nv1z]/cos(Vy1[nv1-3]*vy)
	+16.*F1[i10-nv1z]/cos(Vy1[nv1-2]*vy)-30.*F1[i10]
	+16.*F1[i10+nv1z]/cos(Vy1[nv1]*vy)-F1[i10+2*nv1z]/cos(Vy1[nv1+1]*vy))/rho1[ir];
	
    t1z[ir]=-dv1*creal(-F1[i10-2]/cos(Vz1[nv1-3]*vz)
	+16.*F1[i10-1]/cos(Vz1[nv1-2]*vz)-30.*F1[i10]
	+16.*F1[i10+1]/cos(Vz1[nv1]*vz)-F1[i10+2]/cos(Vz1[nv1+1]*vz))/rho1[ir];

    E1+=-.5*dv1*(2.*creal(-15.*F1[i10]+16.*F1[i10+nv12]-F1[i10+2*nv12])+
	creal(-F1[i10-2*nv1z]+16.*F1[i10-nv1z]-30.*F1[i10]+16.*F1[i10+nv1z]-F1[i10+2*nv1z])+
	creal(-F1[i10-2]+16.*F1[i10-1]-30.*F1[i10]+16.*F1[i10+1]-F1[i10+2]))/Q1;

	Et1+=(t1x[ir]+t1y[ir]+t1z[ir])*rho1[ir]; Px1+=jx1[i];
	}
	Et1=.5*Et1/Q1; Ek1=E1-Et1; Etot=E1+EE+EB; Px1=Px1/QM1; Pxtot=Px1;
}
#endif
////////////////////////////////////////////////////////////////////////////////
///////////////////// Calculation of the electric fields ///////////////////////
////////////////////////////////////////////////////////////////////////////////
#if defined(_HYBRID)
void fields0(void){ //fields of the initial current perturbation, hybrid case

	fftw_execute(transjyk); fftw_execute(transjzk);
	for(i=1;i<nw+1;i++){
		ayk[i]=jyk[i]/(i*i*TPI*TPI/nx/nx/dx/dx-i*TPI/dx/nx*rhoe[i]*(uex[i]-1.2*vA)/Bx0)/nx; 
		azk[i]=jzk[i]/(i*i*TPI*TPI/nx/nx/dx/dx-i*TPI/dx/nx*rhoe[i]*(uex[i]-1.2*vA)/Bx0)/nx;}
	for(i=nw+1;i<nx;i++){ayk[i]=0.; azk[i]=0.;}
	ayk[0]=0.; azk[0]=0.; fftw_execute(transay0); fftw_execute(transaz0);	
	
	for(i=0;i<nx;i++){pe[i]=.01*rhoe[i]*T1;}

	dfdx(Ay,Bz,1.); dfdx(Az,By,-1.); dfdx(pe,dpe,1.);
		
	EE=EB=Pw=0.;
	for(i=0;i<nx;i++){Ey[i]=-vA*Bz[i]; Ez[i]=vA*By[i];
		//Ex[i]=-(Ey[i]*By[i]+Ez[i]*Bz[i])/Bx0+dpe[i]/rhoe[i];
		Ex[i]=-(Ey[i]*By[i]+Ez[i]*Bz[i])/Bx[i];
		By[i]=By[i]+By0; Bz[i]=Bz[i]+Bz0;
		Pw+=(Ey[i]*Bz[i]-Ez[i]*By[i]); 
		EE+=Ex[i]*Ex[i]+Ey[i]*Ey[i]+Ez[i]*Ez[i]; EB+=By[i]*By[i]+Bz[i]*Bz[i];}
	EE=.5*EE/c2; EB=.5*EB; Pw=Pw/c2; Pxtot+=Pw;
	#if defined(_1D)
	trans1=fftw_plan_dft_r2c_1d(nx,pe,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,pe,FFTW_MEASURE);

	fftw_execute(trans1); //the DC pe field does not have to be 0!
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ex,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ex,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ex field does not have to be 0!
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ey,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ey,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ey = 0
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ez,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ez,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ez = 0
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,By,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,By,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Bz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Bz,FFTW_MEASURE);

	fftw_execute(trans1);
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif	
}
#else
fields0(void){ //fields of the initial current perturbation, generic/electron case
	
	fftw_execute(transjxk); 
	fftw_execute(transjyk); fftw_execute(transjzk);
	for(i=1;i<nw+1;i++){
		//ayk[i]=jyk[i]/(i*i*TPI*TPI/nx/nx/dx/dx*(c*c-.25*vA*vA))/nx; 
		//azk[i]=jzk[i]/(i*i*TPI*TPI/nx/nx/dx/dx*(c*c-.25*vA*vA))/nx;}
		//outk[i]=-ink[i]/(2.*i*i*TPI*TPI/nx/nx/dx/dx*c*c)/nx; 
		//outk[i]=-ink[i]/(i*i*TPI*TPI/nx/nx/dx/dx*c*c-1.)/nx; 
		outk[i]=-ink[i]/(i*i*TPI*TPI/nx/nx/dx/dx*c2)/nx; 
		ayk[i]=jyk[i]/(i*i*TPI*TPI/nx/nx/dx/dx*c2)/nx; 
		azk[i]=jzk[i]/(i*i*TPI*TPI/nx/nx/dx/dx*c2)/nx;}
	for(i=nw+1;i<nx;i++){ayk[i]=0.; azk[i]=0.;}
	ayk[0]=0.; azk[0]=0.; fftw_execute(transay0); fftw_execute(transaz0);	
	outk[0]=0.; fftw_execute(transax0); //this initializes Ax out of jx
	
	dfdx(Ax,Gx,1.); //this prepares Gx with the assumpthion of omega = w_pe
	
	dfdx(Gx,dGx,1.);
		
	for(i=0;i<nr;i++){dphi[i]=-rho[i]-dGx[i];}
	
#if defined(_3D)
	Poisson(dphi,phi,-TPI/nr); 
	dfdy(phi,Ey,-1.); dfdz(phi,Ez,-1.); dfdy(Az,Bx,1.); 
	dfdz(Ay,dAydz,1.); dfdy(Ax,dAxdy,1.); dfdz(Ax,dAxdz,1.);
#elif defined(_1D)
	Poisson(dphi,phi,coef3); 
#endif	
	
	dfdx(phi,Ex,-1.); dfdx(Ay,Bz,1.); dfdx(Az,By,-1.);
	
 	for(i=0;i<nr;i++){
#if defined(_3D)
	Bx[i]+=Bx0-dAydz[i]; By[i]+=By0+dAxdz[i]; Bz[i]+=Bz0-dAxdy[i];
#elif defined(_1D)
	Ex[i]-=Gx[i]; Ey[i]=-vA*Bz[i]; Ez[i]=vA*By[i]; 
	Bx[i]=Bx0; By[i]+=By0; Bz[i]+=Bz0;
#endif
	Gy[i]=-Ey[i]; Gz[i]=-Ez[i]; }	
}
#endif

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// MAXWELL SOLVER ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if defined(_HYBRID)
void fields(double *ay, double *az, double *pe, double *pet){
	dfdx(az,By,-1); dfdx(ay,Bz,1); d2fdx(az,d2Az,1); d2fdx(ay,d2Ay,1);
	dfdx(pe,dpe,1); dfdx(ubx,dubx,1.);	
	
	EE=EB=Pw=0.;
	for(i=0;i<nr;i++){ // calcualte E 		
		Ey[i]=uex[i]*Bz[i]-Bx[i]/rhoe[i]*(jz[i]+d2Az[i]); 
		Ez[i]=-uex[i]*By[i]+Bx[i]/rhoe[i]*(jy[i]+d2Ay[i]);
		//Ex[i]=-(Ey[i]*By[i]+Ez[i]*Bz[i])/Bx[i]+dpe[i]/rhoe[i];
		Ex[i]=-(Ey[i]*By[i]+Ez[i]*Bz[i])/Bx[i];
		//Ex[i]=0.;
		By[i]=By[i]+By0; Bz[i]=Bz[i]+Bz0;
		pet[i]=-(dubx[i]*pe[i]+ge*ubx[i]*dpe[i]);
		Pw+=(Ey[i]*Bz[i]-Ez[i]*By[i]); 
		EE+=Ex[i]*Ex[i]+Ey[i]*Ey[i]+Ez[i]*Ez[i]; EB+=By[i]*By[i]+Bz[i]*Bz[i];}
	EE=.5*EE/c2; EB=.5*EB; Pw=Pw/c2; Pxtot+=Pw;
		
	#if defined(_1D)	
	for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	trans1=fftw_plan_dft_r2c_1d(nx,pet,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,pet,FFTW_MEASURE);

	fftw_execute(trans1); //the DC pet field does not have to be 0!
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ex,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ex,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ex field does not have to be 0!
	//foutk[0]=0.; 
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ey,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ey,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ey = 0
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Ez,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Ez,FFTW_MEASURE);

	fftw_execute(trans1); //the DC Ez = 0
	foutk[0]=0.; for(i=1;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,By,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,By,FFTW_MEASURE);

	fftw_execute(trans1); //average By != 0
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	
	trans1=fftw_plan_dft_r2c_1d(nx,Bz,fink,FFTW_MEASURE);
	trans2=fftw_plan_dft_c2r_1d(nx,foutk,Bz,FFTW_MEASURE);

	fftw_execute(trans1); //average Bz != 0
	for(i=0;i<nw+1;i++){foutk[i]=fink[i]/nx;}
	//for(i=nw+1;i<nx;i++){foutk[i]=0.;}
	fftw_execute(trans2);
	#endif
}

#else
//the new Maxwell field solver for the electron/general case goes here

#endif
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// VLASOV SOLVER ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void vlasov(double complex *fnew, double complex *f, double *Vx, double *Vy, double *Vz, double QM, double dV, int nv, int nvy, int nvz, double *ay, double *az){
	int nwav, nvyz, nvxyz, nzvxyz, nyzvxyz, nvyyz, nzvyyz, nyzvyyz;
	//old 	  nv2    nv3     n2       n1     n3       n4 		n5
	double diss;
	bool flag1, flag2;
	
	fact2=QM*I; dv=.5/dV; dv2=1./60./dV; dv3=3./dV; diss=0.0001/dx/dx;
	nvyz=nvy*nvz; nvxyz=nv*nvyz; nzvxyz=nz*nvxyz; nyzvxyz=ny*nzvxyz; // dimensions of f
	nvyyz=nvy*nvyz; nzvyyz=nz*nvyyz; nyzvyyz=ny*nzvyyz; //dimensions of G
	//zeit();
	#pragma omp parallel for private(ix,ivy,ivz,i,b)
	for(ix=0;ix<nx;ix++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ix*nvyz+ivy*nvz+ivz; 
		b=-QM*(az[ix]*Vz[ivz]+Vy[ivy]*ay[ix]);
		bp[i]=cos(b)+I*sin(b); bm[i]=cos(b)-I*sin(b);
	}
	//printf("\n b[i] initialized in %1.3e sec ", zeit());

	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){kx1[ivy*nvz+ivz]=QM*(By0*Vz[ivz]-Vy[ivy]*Bz0)/nx;}

	#pragma omp parallel for private(ix,ivx,iv,i,il,ir)
	for(iv=0;iv<nvyz;iv++){
	for(ivx=0;ivx<nv;ivx++){ir=ivx*nvyz+iv;
		il=(ivx+nv-1)*nvyz+iv; 
		for(ix=0;ix<nx;ix++){i=ix*nvyz+iv; f1[iv*nx+ix]=f[ix*nyzvxyz+ir]*bm[i];}
		fftw_execute(fdxf[iv]);
		for(ix=0;ix<nw+1;ix++){G[ix*nyzvyyz+il]=I*f2[iv*nx+ix]*(kx[ix]-kx1[iv]);}
		for(ix=nx-nw;ix<nx;ix++){G[ix*nyzvyyz+il]=I*f2[iv*nx+ix]*(kx[ix]-kx1[iv]);} 
	}
	}
	
	//printf("\n fourier transform of f taken in %1.3e sec ", zeit());

	#pragma omp parallel for private(ix,ivx,ivy,ivz,ir,i,i0)
	for(ivx=0;ivx<nv;ivx++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){ir=ivx*nvyz+ivy*nvz+ivz;
		//for kx=0, ix=0
		i0=(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[ir]=creal(G[i0])-I*cimag(G[i0]);
		for(ix=1;ix<nw+1;ix++){i=ix*nyzvyyz+ir; 
		i0=(nx-ix)*nyzvyyz+(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[i]=creal(G[i0])-I*cimag(G[i0]);}
		for(ix=nx-nw;ix<nx;ix++){i=ix*nyzvyyz+ir; 
		i0=(nx-ix)*nyzvyyz+(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[i]=creal(G[i0])-I*cimag(G[i0]);}
	}	

	//printf("\n f extrapolated to -vx in %1.3e sec ", zeit());

	#pragma omp parallel for private(ix,iy,iz,ivx,ivy,ivz,i,i0,A,B,C,D,bet,nwav,fact1,flag1,flag2)
 	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){fact1=QM*(By0*Vz[ivz]-Vy[ivy]*Bz0); 
		
  	for(ix=0;ix<nw+1;ix++){
		if(ix==0){
			if(fact1<0){nwav=ceil(fact1*nx*dx/TPI);}
			if(fact1>0){nwav=-ceil(-fact1*nx*dx/TPI);}
			if(fact1==0){nwav=0;}
		}
		else{nwav=(int)(fact1*nx*dx/TPI);}

		if(abs(nwav)>nx/2-1){nwav=sign(nwav)*(nx/2-1);}
		
		if(ix==0){ //boundary conditions for modes kx=0
			if(nwav>0){flag1=true; G[ivy*nvz+ivz]=0.; flag2=false;}
			else if(nwav<0){flag2=true; G[(nvy-1)*nvyz+ivy*nvz+ivz]=0.; flag1=false;}
			else if(nwav==0){flag1=flag2=false;}
		}
		else{ //boundary conditions for the mode kx!=0
			if(nwav>=0){ 
				if(ix>=nwav&&ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
				if(ix<nwav||ix>nx/2+1){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}
			}
			else{
				if(ix>nx/2+1&&ix<=nx-1+nwav){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}	
				if(ix>=nx-1+nwav||ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
			}
		}

		i0=ix*nvyyz+ivy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 &&flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]=dft[i];}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]+=I*dft[i];}
	}
			
   	for(ix=nx-nw;ix<nx;ix++){
		nwav=(int)(fact1*nx*dx/TPI);

		if(abs(nwav)>nx/2-1){nwav=sign(nwav)*(nx/2-1);}
		
		if(ix==0){ //boundary conditions for modes kx=0
			if(nwav>0){flag1=true; G[ivy*nvz+ivz]=0.; flag2=false;}
			else if(nwav<0){flag2=true; G[(nvy-1)*nvyz+ivy*nvz+ivz]=0.; flag1=false;}
			else if(nwav==0){flag1=flag2=false;}
		}
		else{ //boundary conditions for the mode kx!=0
			if(nwav>=0){ 
				if(ix>=nwav&&ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
				if(ix<nwav||ix>nx/2+1){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}
			}
			else{
				if(ix>nx/2+1&&ix<=nx-1+nwav){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}	
				if(ix>=nx-1+nwav||ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
			}
		}

		i0=ix*nvyyz+ivy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 &&flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]=dft[i];}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]+=I*dft[i];}
	}
	}
	
	//printf("\n f differentiated over vx in %1.3e sec ", zeit());

	#pragma omp parallel for private(ix,ivx,iv,il,ir)
	for(iv=0;iv<nvyz;iv++)
 	for(ivx=0;ivx<nv;ivx++){ir=ivx*nvyz+iv; il=(ivx+nv-1)*nvyz+iv;
		//for(ix=0;ix<nx;ix++){f1[iv*nx+ix]=dfdv[ix*n5+il];} fftw_execute(fdxb[iv]); 
		for(ix=0;ix<nw+1;ix++){f1[iv*nx+ix]=dfdv[ix*nyzvyyz+il];} 
		for(ix=nw+1;ix<nx-nw;ix++){f1[iv*nx+ix]=0.;} 
		for(ix=nx-nw;ix<nx;ix++){f1[iv*nx+ix]=dfdv[ix*nyzvyyz+il];} fftw_execute(fdxb[iv]); 
		for(ix=0;ix<nx;ix++){dfdvy[ix*nyzvxyz+ir]=I*f2[iv*nx+ix]*bp[ix*nvyz+iv];}
	}	
	
	//printf("\n df/dvx fourier-transormed in %1.3e sec ", zeit());
	
	//add dissipative term d6f/dx6*delta*dx4
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i)
   	for(ix=3;ix<nx-3;ix++) //differentiation of f over x   
	for(ivx=0;ivx<nv;ivx++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ix*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=0
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i+(nx-2)*nvxyz]+15.*f[i+(nx-1)*nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=1
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i+(nx-2)*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=2
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=2*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-1
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-1)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i-(nx-1)*nvxyz]-6.*f[i-(nx-2)*nvxyz]+f[i-(nx-3)*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-2
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-2)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i-(nx-2)*nvxyz]+f[i-(nx-3)*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-3
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-3)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i-(nx-3)*nvxyz]);}
	
	
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i,i0,A,B,C,D,bet,flag1,flag2,fact1)
   	for(ix=0;ix<nx;ix++) //differentiation of G over vz   
	for(ivx=0;ivx<nv;ivx++){
	for(ivy=0;ivy<nvy;ivy++){ fact1=QM*(Bx[ix]*Vy[ivy]-Vx[ivx]*By[ix]); 
		
		if(fact1>0){ //ivy<nvy/2
			flag1=false; flag2=true; f[ix*nvxyz+ivx*nvyz+ivy*nvz+nvz-1]=0.;}
		else if(fact1<0){ //ivy>nvy/2
			flag2=false; flag1=true; f[ix*nvxyz+ivx*nvyz+ivy*nvz]=0.;}
		else if(fact1==0){flag1=flag2=false;}
		
		i0=ix*nvxyz+ivx*nvyz+ivy*nvz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*f[i0]-4.*f[i0+1]-f[i0+2]);}
		bet=B; dft[i0]=D/bet;
		for(ivz=1;ivz<nvz;ivz++){i=i0+ivz;  //foward substitution
			if(ivz>1 && ivz<nvz-1){A=1.;B=4.;C=1.;D=dv3*(creal(f[i+1]-f[i-1]));}
			else if(ivz==nvz-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivz==nvz-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*f[i]-4.*f[i-1]-f[i-2]);}
			else if(ivz==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(f[i+1]-f[i-1]));}
			else if(ivz==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(f[i+1]-f[i-1]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-1])/bet;
		}	//backward substitutioin
      	for(ivz=nvz-2;ivz>=0;ivz--){i=i0+ivz; dft[i]-=gam[i+1]*dft[i+1];}
		for(ivz=0;ivz<nvz;ivz++){i=i0+ivz; dfdvy[i]+=dft[i]*fact1;}		
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*f[i0]-4.*f[i0+1]-f[i0+2]);}
		bet=B; dft[i0]=D/bet;
		for(ivz=1;ivz<nvz;ivz++){i=i0+ivz;  //foward substitution
			if(ivz>1 && ivz<nvz-1){A=1.;B=4.;C=1.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			else if(ivz==nvz-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivz==nvz-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*f[i]-4.*f[i-1]-f[i-2]);}
			else if(ivz==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			else if(ivz==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-1])/bet;
		}	//backward substitutioin
      	for(ivz=nvz-2;ivz>=0;ivz--){i=i0+ivz; dft[i]-=gam[i+1]*dft[i+1];}
		for(ivz=0;ivz<nvz;ivz++){i=i0+ivz; dfdvy[i]+=I*dft[i]*fact1;}
	}

	for(ivz=0;ivz<nvz;ivz++){ fact1=QM*(Bz[ix]*Vx[ivx]-Vz[ivz]*Bx[ix]); 
		
		if(fact1<0){
			flag1=true; f[ix*nvxyz+ivx*nvyz+ivz]=0.; flag2=false;}
		else if(fact1>0){
			flag2=true; f[ix*nvxyz+ivx*nvyz+(nvy-1)*nvz+ivz]=0.; flag1=false;}
		else if(fact1==0){flag1=flag2=false;}
		
		i0=ix*nvxyz+ivx*nvy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*f[i0]-4.*f[i0+nvz]-f[i0+2*nvz]);}
		bet=B; dft[i0]=D/bet;
		for(ivy=1;ivy<nvy;ivy++){i=i0+ivy*nvz;  //foward substitution
			if(ivy>1 && ivy<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			else if(ivy==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivy==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*f[i]-4.*f[i-nvz]-f[i-2*nvz]);}
			else if(ivy==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			else if(ivy==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvz])/bet;
		}	//backward substitutioin
      	for(ivy=nvy-2;ivy>=0;ivy--){i=i0+ivy*nvz; dft[i]-=gam[i+nvz]*dft[i+nvz];}
		for(ivy=0;ivy<nvy;ivy++){i=i0+ivy*nvz; dfdvy[i]+=dft[i]*fact1;}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*f[i0]-4.*f[i0+nvz]-f[i0+2*nvz]);}
		bet=B; dft[i0]=D/bet;
		for(ivy=1;ivy<nvy;ivy++){i=i0+ivy*nvz;  //foward substitution
			if(ivy>1 && ivy<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			else if(ivy==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivy==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*f[i]-4.*f[i-nvz]-f[i-2*nvz]);}
			else if(ivy==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			else if(ivy==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvz])/bet;
		}	//backward substitutioin
      	for(ivy=nvy-2;ivy>=0;ivy--){i=i0+ivy*nvz; dft[i]-=gam[i+nvz]*dft[i+nvz];}
		for(ivy=0;ivy<nvy;ivy++){i=i0+ivy*nvz; 
		fnew[i]=dfdvy[i]+I*dft[i]*fact1+fact2*(Ex[ix]*Vx[ivx]+Ey[ix]*Vy[ivy]+Ez[ix]*Vz[ivz])*f[i];}
	}
	}
	
	//printf("\n f differentiated over vy and vz and added to fnew in %1.3e sec ", zeit());
}

#ifdef _NONINERT
////////////////////////////////////////////////////////////////////////////////
///////////////// VLASOV SOLVER FOR THE DRIFTING POPULATION ////////////////////
////////////////////////////////////////////////////////////////////////////////
void vlasov2(double complex *fnew, double complex *f, double *Vx, double *Vy, double *Vz, double QM, double dV, int nv, int nvy, int nvz, double *ay, double *az){
	int nwav, nvyz, nvxyz, nzvxyz, nyzvxyz, nvyyz, nzvyyz, nyzvyyz;
	//old 	  nv2    nv3     n2       n1     n3       n4      n5
	double diff, diss;
	bool flag1, flag2;
	
	fact2=QM*I; dv=.5/dV; dv2=1./60./dV; dv3=3./dV; diff=0.5*udave/dx; diss=0.0001/dx/dx;
	nvyz=nvy*nvz; nvxyz=nv*nvyz; nzvxyz=nz*nvxyz; nyzvxyz=ny*nzvxyz; // dimensions of f
	nvyyz=nvy*nvyz; nzvyyz=nz*nvyyz; nyzvyyz=ny*nzvyyz; //dimensions of G

	//zeit();
	
	#pragma omp parallel for private(ix,ivy,ivz,i,b)
	for(ix=0;ix<nx;ix++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ix*nvyz+ivy*nvz+ivz; 
		b=-QM*(az[ix]*Vz[ivz]+Vy[ivy]*ay[ix])-(Vy[ivy]*uyave[ix]+Vz[ivz]*uzave[ix]);		
		bp[i]=cos(b)+I*sin(b); bm[i]=cos(b)-I*sin(b);
	}
	//printf("\n b[i] initialized in %1.3e sec ", zeit());

	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){kx1[ivy*nvz+ivz]=QM*(By0*Vz[ivz]-Vy[ivy]*Bz0)/nx;}

	#pragma omp parallel for private(ix,ivx,iv,i,il,ir)
	for(iv=0;iv<nvyz;iv++){
	for(ivx=0;ivx<nv;ivx++){ir=ivx*nvyz+iv;
		il=(ivx+nv-1)*nvyz+iv; 
		for(ix=0;ix<nx;ix++){i=ix*nvyz+iv; f1[iv*nx+ix]=f[ix*nyzvxyz+ir]*bm[i];}
		fftw_execute(fdxf[iv]);
		for(ix=0;ix<nw+1;ix++){G[ix*nyzvyyz+il]=I*f2[iv*nx+ix]*(kx[ix]-kx1[iv]);}
		for(ix=nx-nw;ix<nx;ix++){G[ix*nyzvyyz+il]=I*f2[iv*nx+ix]*(kx[ix]-kx1[iv]);} 
	}
	}

	//printf("\n fourier transform of f taken in %1.3e sec ", zeit());

	//dFdx(f,dfdvy,nve,diff);
	
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i)
   	for(ix=1;ix<nx-1;ix++) //differentiation of f over x   
	for(ivx=0;ivx<nv;ivx++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ix*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]=-diff*(f[i+nvxyz]-f[i-nvxyz]);}

	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++)   	//differentiation of f over x   
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-1)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]=-diff*(f[i-(nx-1)*nvxyz]-f[i-nvxyz]);}

	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++)   	//differentiation of f over x   
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]=-diff*(f[i+nvxyz]-f[i+(nx-1)*nvxyz]);}
	
	//add dissipative term d6f/dx6*delta*dx4
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i)
   	for(ix=3;ix<nx-3;ix++) //differentiation of f over x   
	for(ivx=0;ivx<nv;ivx++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ix*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=0
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i+(nx-2)*nvxyz]+15.*f[i+(nx-1)*nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=1
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i+(nx-2)*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=2
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=2*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i+(nx-3)*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i+3*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-1
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-1)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i-(nx-1)*nvxyz]-6.*f[i-(nx-2)*nvxyz]+f[i-(nx-3)*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-2
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-2)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i-(nx-2)*nvxyz]+f[i-(nx-3)*nvxyz]);}
	
	#pragma omp parallel for private(ivx,ivy,ivz,i)
	for(ivx=0;ivx<nv;ivx++) //ix=nx-3
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){i=(nx-3)*nvxyz+ivx*nvyz+ivy*nvz+ivz;
		dfdvy[i]+=diss*(f[i-3*nvxyz]-6.*f[i-2*nvxyz]+15.*f[i-nvxyz]-20.*f[i]+15.*f[i+nvxyz]-6.*f[i+2*nvxyz]+f[i-(nx-3)*nvxyz]);}

	#pragma omp parallel for private(ix,ivx,ivy,ivz,ir,i,i0)
	for(ivx=0;ivx<nv;ivx++)
	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){ir=ivx*nvyz+ivy*nvz+ivz;
		//for kx=0, ix=0
		i0=(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[ir]=creal(G[i0])-I*cimag(G[i0]);
		for(ix=1;ix<nw+1;ix++){i=ix*nyzvyyz+ir; 
		i0=(nx-ix)*nyzvyyz+(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[i]=creal(G[i0])-I*cimag(G[i0]);}
		for(ix=nx-nw;ix<nx;ix++){i=ix*nyzvyyz+ir; 
		i0=(nx-ix)*nyzvyyz+(nvy-1-ivx)*nvyz+(nvy-1-ivy)*nvz+nvz-1-ivz;
		G[i]=creal(G[i0])-I*cimag(G[i0]);}
	}
	
	//printf("\n f extrapolated to -vx in %1.3e sec ", zeit());

	#pragma omp parallel for private(ix,iy,iz,ivx,ivy,ivz,i,i0,A,B,C,D,bet,nwav,fact1,flag1,flag2)
 	for(ivy=0;ivy<nvy;ivy++)
	for(ivz=0;ivz<nvz;ivz++){fact1=QM*(By0*Vz[ivz]-Vy[ivy]*Bz0); 
		
  	for(ix=0;ix<nw+1;ix++){
		if(ix==0){
			if(fact1<0){nwav=ceil(fact1*nx*dx/TPI);}
			if(fact1>0){nwav=-ceil(-fact1*nx*dx/TPI);}
			if(fact1==0){nwav=0;}
		}
		else{nwav=(int)(fact1*nx*dx/TPI);}

		if(abs(nwav)>nx/2-1){nwav=sign(nwav)*(nx/2-1);}
		
		if(ix==0){ //boundary conditions for modes kx=0
			if(nwav>0){flag1=true; G[ivy*nvz+ivz]=0.; flag2=false;}
			else if(nwav<0){flag2=true; G[(nvy-1)*nvyz+ivy*nvz+ivz]=0.; flag1=false;}
			else if(nwav==0){flag1=flag2=false;}
		}
		else{ //boundary conditions for the mode kx!=0
			if(nwav>=0){ 
				if(ix>=nwav&&ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
				if(ix<nwav||ix>nx/2+1){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}
			}
			else{
				if(ix>nx/2+1&&ix<=nx-1+nwav){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}	
				if(ix>=nx-1+nwav||ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
			}
		}

		i0=ix*nvyyz+ivy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 &&flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]=dft[i];}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]+=I*dft[i];}
	}
			
   	for(ix=nx-nw;ix<nx;ix++){
		nwav=(int)(fact1*nx*dx/TPI);

		if(abs(nwav)>nx/2-1){nwav=sign(nwav)*(nx/2-1);}
		
		if(ix==0){ //boundary conditions for modes kx=0
			if(nwav>0){flag1=true; G[ivy*nvz+ivz]=0.; flag2=false;}
			else if(nwav<0){flag2=true; G[(nvy-1)*nvyz+ivy*nvz+ivz]=0.; flag1=false;}
			else if(nwav==0){flag1=flag2=false;}
		}
		else{ //boundary conditions for the mode kx!=0
			if(nwav>=0){ 
				if(ix>=nwav&&ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
				if(ix<nwav||ix>nx/2+1){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}
			}
			else{
				if(ix>nx/2+1&&ix<=nx-1+nwav){flag2=true; G[ix*nyzvyyz+(nvy-1)*nvyz+ivy*nvz+ivz]=0.;}
				else{flag2=false;}	
				if(ix>=nx-1+nwav||ix<nx/2-1){flag1=true; G[ix*nyzvyyz+ivy*nvz+ivz]=0.;}
				else{flag1=false;}
			}
		}

		i0=ix*nvyyz+ivy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 &&flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]=dft[i];}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*G[i0]-4.*G[i0+nvyz]-G[i0+2*nvyz]);}
		bet=B; dft[i0]=D/bet;
		for(ivx=1;ivx<nvy;ivx++){i=i0+ivx*nvyz;  //foward substitution
			if(ivx>1 && ivx<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivx==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*G[i]-4.*G[i-nvyz]-G[i-2*nvyz]);}
			else if(ivx==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			else if(ivx==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(G[i+nvyz]-G[i-nvyz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvyz])/bet;
		}	//backward substitutioin
      	for(ivx=nvy-2;ivx>=0;ivx--){i=i0+ivx*nvyz; dft[i]-=gam[i+nvyz]*dft[i+nvyz];}
		for(ivx=0;ivx<nvy;ivx++){i=i0+ivx*nvyz; dfdv[i]+=I*dft[i];}
	}
	}
	

	//printf("\n f differentiated over vx in %1.3e sec ", zeit());

	#pragma omp parallel for private(ix,ivx,iv,il,ir)
	for(iv=0;iv<nvyz;iv++)
 	for(ivx=0;ivx<nv;ivx++){ir=ivx*nvyz+iv; il=(ivx+nv-1)*nvyz+iv;
		for(ix=0;ix<nw+1;ix++){f1[iv*nx+ix]=dfdv[ix*nyzvyyz+il];} 
		for(ix=nw+1;ix<nx-nw;ix++){f1[iv*nx+ix]=0.;} 
		for(ix=nx-nw;ix<nx;ix++){f1[iv*nx+ix]=dfdv[ix*nyzvyyz+il];} fftw_execute(fdxb[iv]); 
		for(ix=0;ix<nx;ix++){dfdvy[ix*nyzvxyz+ir]+=I*f2[iv*nx+ix]*bp[ix*nvyz+iv];}
	}	
	
	//printf("\n df/dvx fourier-transormed in %1.3e sec ", zeit());
		
	#pragma omp parallel for private(ix,ivx,ivy,ivz,i,i0,A,B,C,D,bet,flag1,flag2,fact1)
   	for(ix=0;ix<nx;ix++) //differentiation of G over vz   
	for(ivx=0;ivx<nv;ivx++){
	for(ivy=0;ivy<nvy;ivy++){ fact1=QM*(Bx[ix]*Vy[ivy]-Vx[ivx]*By[ix]); 
		
		if(fact1>0){ //ivy<nvy/2
			flag1=false; flag2=true; f[ix*nvxyz+ivx*nvyz+ivy*nvz+nvz-1]=0.;}
		else if(fact1<0){ //ivy>nvy/2
			flag2=false; flag1=true; f[ix*nvxyz+ivx*nvyz+ivy*nvz]=0.;}
		else if(fact1==0){flag1=flag2=false;}
		
		i0=ix*nvxyz+ivx*nvyz+ivy*nvz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*f[i0]-4.*f[i0+1]-f[i0+2]);}
		bet=B; dft[i0]=D/bet;
		for(ivz=1;ivz<nvz;ivz++){i=i0+ivz;  //foward substitution
			if(ivz>1 && ivz<nvz-1){A=1.;B=4.;C=1.;D=dv3*(creal(f[i+1]-f[i-1]));}
			else if(ivz==nvz-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivz==nvz-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*f[i]-4.*f[i-1]-f[i-2]);}
			else if(ivz==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(f[i+1]-f[i-1]));}
			else if(ivz==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(f[i+1]-f[i-1]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-1])/bet;
		}	//backward substitutioin
      	for(ivz=nvz-2;ivz>=0;ivz--){i=i0+ivz; dft[i]-=gam[i+1]*dft[i+1];}
		for(ivz=0;ivz<nvz;ivz++){i=i0+ivz; dfdvy[i]+=dft[i]*fact1;}		
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*f[i0]-4.*f[i0+1]-f[i0+2]);}
		bet=B; dft[i0]=D/bet;
		for(ivz=1;ivz<nvz;ivz++){i=i0+ivz;  //foward substitution
			if(ivz>1 && ivz<nvz-1){A=1.;B=4.;C=1.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			else if(ivz==nvz-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivz==nvz-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*f[i]-4.*f[i-1]-f[i-2]);}
			else if(ivz==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			else if(ivz==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(f[i+1]-f[i-1]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-1])/bet;
		}	//backward substitutioin
      	for(ivz=nvz-2;ivz>=0;ivz--){i=i0+ivz; dft[i]-=gam[i+1]*dft[i+1];}
		for(ivz=0;ivz<nvz;ivz++){i=i0+ivz; dfdvy[i]+=I*dft[i]*fact1;}
	}

	for(ivz=0;ivz<nvz;ivz++){ fact1=QM*(Bz[ix]*Vx[ivx]-Vz[ivz]*Bx[ix]); 
		//differentiation of G over vy
		if(fact1<0){
			flag1=true; f[ix*nvxyz+ivx*nvyz+ivz]=0.; flag2=false;}
		else if(fact1>0){
			flag2=true; f[ix*nvxyz+ivx*nvyz+(nvy-1)*nvz+ivz]=0.; flag1=false;}
		else if(fact1==0){flag1=flag2=false;}
		
		i0=ix*nvxyz+ivx*nvy*nvz+ivz; //real part
		if(flag1){A=0.;B=1.; C=-1.; D=0.;}
		else{A=0.;B=1.; C=2.; D=-dv*creal(5.*f[i0]-4.*f[i0+nvz]-f[i0+2*nvz]);}
		bet=B; dft[i0]=D/bet;
		for(ivy=1;ivy<nvy;ivy++){i=i0+ivy*nvz;  //foward substitution
			if(ivy>1 && ivy<nvy-1){A=1.;B=4.;C=1.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			else if(ivy==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivy==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*creal(5.*f[i]-4.*f[i-nvz]-f[i-2*nvz]);}
			else if(ivy==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			else if(ivy==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(creal(f[i+nvz]-f[i-nvz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvz])/bet;
		}	//backward substitutioin
      	for(ivy=nvy-2;ivy>=0;ivy--){i=i0+ivy*nvz; dft[i]-=gam[i+nvz]*dft[i+nvz];}
		for(ivy=0;ivy<nvy;ivy++){i=i0+ivy*nvz; dfdvy[i]+=dft[i]*fact1;}
				
		if(flag1){A=0.;B=1.; C=-1.; D=0.;} //imaginary part
		else{A=0.;B=1.; C=2.; D=-dv*cimag(5.*f[i0]-4.*f[i0+nvz]-f[i0+2*nvz]);}
		bet=B; dft[i0]=D/bet;
		for(ivy=1;ivy<nvy;ivy++){i=i0+ivy*nvz;  //foward substitution
			if(ivy>1 && ivy<nvy-1){A=1.;B=4.;C=1.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			else if(ivy==nvy-1 && flag2){A=-1.;B=1.;C=1.; D=0.;}
			else if(ivy==nvy-1 && !flag2){A=2.;B=1.;C=1.; D=dv*cimag(5.*f[i]-4.*f[i-nvz]-f[i-2*nvz]);}
			else if(ivy==1 && !flag1){C=2.;A=1.;B=4.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			else if(ivy==1 && flag1){C=-1.; A=1.;B=4.;D=dv3*(cimag(f[i+nvz]-f[i-nvz]));}
			gam[i]=C/bet; bet=B-A*gam[i]; dft[i]=(D-A*dft[i-nvz])/bet;
		}	//backward substitutioin
      	for(ivy=nvy-2;ivy>=0;ivy--){i=i0+ivy*nvz; dft[i]-=gam[i+nvz]*dft[i+nvz];}
		for(ivy=0;ivy<nvy;ivy++){i=i0+ivy*nvz; 
		fnew[i]=dfdvy[i]+I*dft[i]*fact1
			+(fact2*(Ex[ix]*Vx[ivx]+Ey[ix]*Vy[ivy]+Ez[ix]*Vz[ivz]
			+udave*(By[ix]*Vz[ivz]-Bz[ix]*Vy[ivy])
			+uyave[ix]*(Bz[ix]*Vx[ivx]-Bx[ix]*Vz[ivz])+uzave[ix]*(Bx[ix]*Vy[ivy]-By[ix]*Vx[ivx])
			)
			-I*(Vx[ivx]*uxt
			+Vy[ivy]*(duy[ix]*udave+uyt[ix])+Vz[ivz]*(duz[ix]*udave+uzt[ix])
			)
			)*f[i];}
	}
	}
	
	//printf("\n f differentiated over vy and vz and added to fnew in %1.3e sec ", zeit());
}
#endif


#if defined(_3POP)
void integrateRK4(){	//Runge-Kutta time advancement 

	moments(F1, F2, F3); 
	#if defined(_NONINERT)
	udnew=udave=udold;
	for(i=0;i<nx;i++){
		uynew[i]=uyold[i]; uznew[i]=uzold[i]; uyave[i]=uyold[i]; uzave[i]=uzold[i];
	}
	dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	#endif
	for(i=0;i<nr;i++){pe0[i]=pe[i]; Ay0[i]=Ay[i]; Az0[i]=Az[i]; Ayt[i]=Ay[i]; Azt[i]=Az[i];}
		
	fields(Ay, Az, pe, pet);	
	
	vlasov(Ft1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	vlasov(Ft2, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft3, F3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#else
	vlasov(Ft3, F3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]=F1[i]+dt2*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];} 
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]=F2[i]+dt2*Ft2[i]; Ft2[i]=F2[i]+dt1*Ft2[i];} 
	#pragma omp parallel for private(i)
	for(i=0;i<n31;i++){Fn3[i]=F3[i]+dt2*Ft3[i]; Ft3[i]=F3[i]+dt1*Ft3[i];} 

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]=pe0[i]+dt2*pet[i]; 	
		Ayn[i]=Ay0[i]-dt2*Ey[i]; Azn[i]=Az0[i]-dt2*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}

	moments(Ft1, Ft2, Ft3); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt1*uxt; udnew+=dt2*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt2*uyt[i]; uznew[i]+=dt2*uzt[i];
		uyave[i]=uyold[i]+dt1*uyt[i]; uzave[i]=uzold[i]+dt1*uzt[i];
	}
	dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#else
	vlasov(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]+=dt3*Ft2[i]; Ft2[i]=F2[i]+dt1*Ft2[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n31;i++){Fn3[i]+=dt3*Ft3[i]; Ft3[i]=F3[i]+dt1*Ft3[i];}
	
	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}
	
	moments(Ft1, Ft2, Ft3); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt1*uxt; udnew+=dt3*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt3*uyt[i]; uznew[i]+=dt3*uzt[i];
		uyave[i]=uyold[i]+dt1*uyt[i]; uzave[i]=uzold[i]+dt1*uzt[i];
	}
	dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#else
	vlasov(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt*Ft1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]+=dt3*Ft2[i]; Ft2[i]=F2[i]+dt*Ft2[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n31;i++){Fn3[i]+=dt3*Ft3[i]; Ft3[i]=F3[i]+dt*Ft3[i];}

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt*Ey[i]; Azt[i]=Az0[i]-dt*Ez[i];}
	
	moments(Ft1, Ft2, Ft3); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt*uxt; udnew+=dt3*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt3*uyt[i]; uznew[i]+=dt3*uzt[i];
		uyave[i]=uyold[i]+dt*uyt[i]; uzave[i]=uzold[i]+dt*uzt[i];
	}
	dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#else
	vlasov(Ft3, Ft3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ayt, Azt);  
	#endif
 	#pragma omp parallel for private(i)
  	for(i=0;i<n11;i++){F1[i]=Fn1[i]+dt2*Ft1[i];} 
 	#pragma omp parallel for private(i)
  	for(i=0;i<n21;i++){F2[i]=Fn2[i]+dt2*Ft2[i];} 
 	#pragma omp parallel for private(i)
  	for(i=0;i<n31;i++){F3[i]=Fn3[i]+dt2*Ft3[i];} 

	moments(F1, F2, F3); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udold=udnew+dt2*uxt;
	for(i=0;i<nr;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uyold[i]=uynew[i]+dt2*uyt[i]; uzold[i]=uznew[i]+dt2*uzt[i];
		//uyold[i]=0.; uzold[i]=0.;
	}
	#endif
	for(i=0;i<nr;i++){
		pe[i]=pen[i]+dt2*pet[i];	
		Ay[i]=Ayn[i]-dt2*Ey[i]; Az[i]=Azn[i]-dt2*Ez[i];}
}
#elif defined(_2POP)
void integrateRK4(){	//Runge-Kutta time advancement 
	moments(F1, F2); 
	#if defined(_NONINERT)
	udave=udold; udnew=udold; //uxt=(ud-udold)/dt; 
	for(i=0;i<nx;i++){
		//uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;	
		uynew[i]=uyold[i]; uznew[i]=uzold[i];
		uyave[i]=uyold[i]; uzave[i]=uzold[i];
	}
	//dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	dfdx(uyave,duy,1); dfdx(uzave,duz,1); 
	#endif
	for(i=0;i<nr;i++){pe0[i]=pe[i]; Ay0[i]=Ay[i]; Az0[i]=Az[i]; Ayt[i]=Ay[i]; Azt[i]=Az[i];}
		
	fields(Ay, Az, pe, pet);	
	
	vlasov(Ft1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft2, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#else
	vlasov(Ft2, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]=F1[i]+dt2*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]=F2[i]+dt2*Ft2[i]; Ft2[i]=F2[i]+dt1*Ft2[i];}

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]=pe0[i]+dt2*pet[i]; 	
		Ayn[i]=Ay0[i]-dt2*Ey[i]; Azn[i]=Az0[i]-dt2*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}

	moments(Ft1, Ft2); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt1*uxt; udnew+=dt2*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt2*uyt[i]; uznew[i]+=dt2*uzt[i];
		uyave[i]=uyold[i]+dt1*uyt[i]; uzave[i]=uzold[i]+dt1*uzt[i];
	}
	//dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	dfdx(uyave,duy,1); dfdx(uzave,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#else
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]+=dt3*Ft2[i]; Ft2[i]=F2[i]+dt1*Ft2[i];}
	
	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}
	
	moments(Ft1, Ft2); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt1*uxt; udnew+=dt3*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt3*uyt[i]; uznew[i]+=dt3*uzt[i];
		uyave[i]=uyold[i]+dt1*uyt[i]; uzave[i]=uzold[i]+dt1*uzt[i];
	}
	//dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	dfdx(uyave,duy,1); dfdx(uzave,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#else
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt*Ft1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){Fn2[i]+=dt3*Ft2[i]; Ft2[i]=F2[i]+dt*Ft2[i];}

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt*Ey[i]; Azt[i]=Az0[i]-dt*Ez[i];}
	
	moments(Ft1, Ft2); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udave=udold+dt*uxt; udnew+=dt3*uxt;
	for(i=0;i<nx;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uynew[i]+=dt3*uyt[i]; uznew[i]+=dt3*uzt[i];
		uyave[i]=uyold[i]+dt*uyt[i]; uzave[i]=uzold[i]+dt*uzt[i];
	}
	//dfdx(uyold,duy,1); dfdx(uzold,duz,1); 
	dfdx(uyave,duy,1); dfdx(uzave,duz,1); 
	#endif
	fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ayt, Azt);   
	#if defined(_NONINERT)
	vlasov2(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#else
	vlasov(Ft2, Ft2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ayt, Azt);   
	#endif
 	#pragma omp parallel for private(i)
  	for(i=0;i<n11;i++){F1[i]=Fn1[i]+dt2*Ft1[i];} 
 	#pragma omp parallel for private(i)
  	for(i=0;i<n21;i++){F2[i]=Fn2[i]+dt2*Ft2[i];} 
			
	moments(F1, F2); 
	#if defined(_NONINERT)
	uxt=(ud-udold)/dt; udold=udnew+dt2*uxt;
	for(i=0;i<nr;i++){
		uyt[i]=(uy[i]-uyold[i])/dt; uzt[i]=(uz[i]-uzold[i])/dt;
		uyold[i]=uynew[i]+dt2*uyt[i]; uzold[i]=uznew[i]+dt2*uzt[i];
	}
	#endif
	for(i=0;i<nr;i++){
		pe[i]=pen[i]+dt2*pet[i];	
		Ay[i]=Ayn[i]-dt2*Ey[i]; Az[i]=Azn[i]-dt2*Ez[i];}
}
#elif defined(_1POP)
void integrateRK4(){	//Runge-Kutta time advancement 
	
	moments(F1);  
	
	for(i=0;i<nr;i++){pe0[i]=pe[i]; Ay0[i]=Ay[i]; Az0[i]=Az[i]; Ayt[i]=Ay[i]; Azt[i]=Az[i];}
		
	fields(Ay, Az, pe, pet);	
	
	vlasov(Ft1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]=F1[i]+dt2*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];}

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]=pe0[i]+dt2*pet[i]; 	
		Ayn[i]=Ay0[i]-dt2*Ey[i]; Azn[i]=Az0[i]-dt2*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}

	moments(Ft1); fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt1*Ft1[i];}
	
	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt1*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt1*Ey[i]; Azt[i]=Az0[i]-dt1*Ez[i];}
	
	moments(Ft1); fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){Fn1[i]+=dt3*Ft1[i]; Ft1[i]=F1[i]+dt*Ft1[i];}

	for(i=0;i<nr;i++){
		pe[i]=pe0[i]+dt*pet[i]; pen[i]+=dt3*pet[i]; 	
		Ayn[i]-=dt3*Ey[i]; Azn[i]-=dt3*Ez[i];
		Ayt[i]=Ay0[i]-dt*Ey[i]; Azt[i]=Az0[i]-dt*Ez[i];}
	
	moments(Ft1); fields(Ayt, Azt, pe, pet);	

	vlasov(Ft1, Ft1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
 	#pragma omp parallel for private(i)
  	for(i=0;i<n11;i++){F1[i]=Fn1[i]+dt2*Ft1[i];} 

	for(i=0;i<nr;i++){
		pe[i]=pen[i]+dt2*pet[i];	
		Ay[i]=Ayn[i]-dt2*Ey[i]; Az[i]=Azn[i]-dt2*Ez[i];}

}
#endif
#if defined(_AB4)
#if defined(_3POP)
void integrateAB4(){	//Adams-Bashforth time advancement 
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF4[i]=dF3[i]; dF3[i]=dF2[i]; dF2[i]=dF1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){dF24[i]=dF23[i]; dF23[i]=dF22[i]; dF22[i]=dF21[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n31;i++){dF34[i]=dF33[i]; dF33[i]=dF32[i]; dF32[i]=dF31[i];}

	#if defined(_NONINERT)
	ud0=udold; dud4=dud3; dud3=dud2; dud2=dud1;
	for(i=0;i<nr;i++){uy0[i]=uyold[i]; uz0[i]=uzold[i];}
	#endif
	moments(F1,F2,F3); fields(Ay, Az, pe, pet);	
	vlasov(dF1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	vlasov(dF21, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ay, Az);   
	#if defined(_NONINERT)
	//ud0=udold; dud4=dud3; dud3=dud2; dud2=dud1;
	//for(i=0;i<nr;i++){uy0[i]=uyold[i]; uz0[i]=uzold[i];}
	vlasov2(dF31, F3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ay, Az);   
	#else
	vlasov(dF31, F3, Vx3, Vy3, Vz3, QM3, dV3, nv3, nv3y, nv3z, Ay, Az);   
	#endif
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF1[i]=dF1[i]*dt; 
		F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){dF21[i]=dF21[i]*dt; 
		F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n31;i++){dF31[i]=dF31[i]*dt; 
		F3[i]+=B1*dF31[i]+B2*dF32[i]+B3*dF33[i]+B4*dF34[i];}
	for(i=0;i<nr;i++){
		dAy4[i]=dAy3[i]; dAy3[i]=dAy2[i]; dAy2[i]=dAy1[i]; dAy1[i]=-dt*Ey[i]; 
		dAz4[i]=dAz3[i]; dAz3[i]=dAz2[i]; dAz2[i]=dAz1[i]; dAz1[i]=-dt*Ez[i];
		Ay[i]+=B1*dAy1[i]+B2*dAy2[i]+B3*dAy3[i]+B4*dAy4[i]; 
		Az[i]+=B1*dAz1[i]+B2*dAz2[i]+B3*dAz3[i]+B4*dAz4[i];}
	#if defined(_NONINERT)
	dud1=ud-ud0; udold=ud0+B1*dud1+B2*dud2+B3*dud3+B4*dud4;
	for(i=0;i<nr;i++){
		duy4[i]=duy3[i]; duy3[i]=duy2[i]; duy2[i]=duy1[i]; duy1[i]=uy[i]-uy0[i]; 
		duz4[i]=duz3[i]; duz3[i]=duz2[i]; duz2[i]=duz1[i]; duz1[i]=uz[i]-uz0[i];
		uyold[i]=ud0[i]+B1*duy1[i]+B2*duy2[i]+B3*duy3[i]+B4*duy4[i];
		uzold[i]=ud0[i]+B1*duz1[i]+B2*duz2[i]+B3*duz3[i]+B4*duz4[i];}
	#endif
}
#elif defined(_2POP)
void integrateAB4(){	//Adams-Bashforth time advancement 
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF4[i]=dF3[i]; dF3[i]=dF2[i]; dF2[i]=dF1[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){dF24[i]=dF23[i]; dF23[i]=dF22[i]; dF22[i]=dF21[i];}
	#if defined(_NONINERT)
	ud0=udold; dud4=dud3; dud3=dud2; dud2=dud1;
	for(i=0;i<nr;i++){uy0[i]=uyold[i]; uz0[i]=uzold[i];}
	#endif
	moments(F1,F2); fields(Ay, Az, pe, pet);	
	vlasov(dF1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	#if defined(_NONINERT)
	//ud0=udold; dud4=dud3; dud3=dud2; dud2=dud1;
	//for(i=0;i<nr;i++){uy0[i]=uyold[i]; uz0[i]=uzold[i];}
	vlasov2(dF21, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ay, Az);   
	#else
	vlasov(dF21, F2, Vx2, Vy2, Vz2, QM2, dV2, nv2, nv2y, nv2z, Ay, Az);   
	#endif	
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF1[i]=dF1[i]*dt; 
		F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
	#pragma omp parallel for private(i)
	for(i=0;i<n21;i++){dF21[i]=dF21[i]*dt; 
		F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
	for(i=0;i<nr;i++){
		dAy4[i]=dAy3[i]; dAy3[i]=dAy2[i]; dAy2[i]=dAy1[i]; dAy1[i]=-dt*Ey[i]; 
		dAz4[i]=dAz3[i]; dAz3[i]=dAz2[i]; dAz2[i]=dAz1[i]; dAz1[i]=-dt*Ez[i];
		Ay[i]+=B1*dAy1[i]+B2*dAy2[i]+B3*dAy3[i]+B4*dAy4[i]; 
		Az[i]+=B1*dAz1[i]+B2*dAz2[i]+B3*dAz3[i]+B4*dAz4[i];}
	#if defined(_NONINERT)
	dud1=ud-ud0; udold=ud0+B1*dud1+B2*dud2+B3*dud3+B4*dud4;
	for(i=0;i<nr;i++){
		duy4[i]=duy3[i]; duy3[i]=duy2[i]; duy2[i]=duy1[i]; duy1[i]=uy[i]-uy0[i]; 
		duz4[i]=duz3[i]; duz3[i]=duz2[i]; duz2[i]=duz1[i]; duz1[i]=uz[i]-uz0[i];
		uyold[i]=uy0[i]+B1*duy1[i]+B2*duy2[i]+B3*duy3[i]+B4*duy4[i];
		uzold[i]=uz0[i]+B1*duz1[i]+B2*duz2[i]+B3*duz3[i]+B4*duz4[i];}
	#endif
}
#elif defined(_1POP)
void integrateAB4(){	//Adams-Bashforth time advancement 
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF4[i]=dF3[i]; dF3[i]=dF2[i]; dF2[i]=dF1[i];}
	moments(F1); fields(Ay, Az, pe, pet);	
	vlasov(dF1, F1, Vx1, Vy1, Vz1, QM1, dV1, nv1, nv1y, nv1z, Ay, Az);   
	#pragma omp parallel for private(i)
	for(i=0;i<n11;i++){dF1[i]=dF1[i]*dt; 
		F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
	for(i=0;i<nr;i++){
		dAy4[i]=dAy3[i]; dAy3[i]=dAy2[i]; dAy2[i]=dAy1[i]; dAy1[i]=-dt*Ey[i]; 
		dAz4[i]=dAz3[i]; dAz3[i]=dAz2[i]; dAz2[i]=dAz1[i]; dAz1[i]=-dt*Ez[i];
		Ay[i]+=B1*dAy1[i]+B2*dAy2[i]+B3*dAy3[i]+B4*dAy4[i]; 
		Az[i]+=B1*dAz1[i]+B2*dAz2[i]+B3*dAz3[i]+B4*dAz4[i];}
}
#endif
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////// MAIN FUNCTION ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
main(void){
////////////////////////////////////////////////////////////////////////////////
////////////////////////// INITIAL PARAMETERS READING //////////////////////////
////////////////////////////////////////////////////////////////////////////////
if((sinit = fopen("param.init","r")) == NULL){ //check if init file is there
	printf("\n ERROR: can't find param.init \n"); exit(1);}
//scan the generic simulation parameters first
fscanf(sinit,"%s %i",&tag[0],&nx);			//number of grids in x
fscanf(sinit,"%s %i",&tag[0],&ny);			//number of grids in y
fscanf(sinit,"%s %i",&tag[0],&nz);			//number of grids in z
fscanf(sinit,"%s %lf",&tag[0],&xmax);		//size of domain in x
fscanf(sinit,"%s %lf",&tag[0],&end_t);		//total runtime in plasma/gyro-periods
fscanf(sinit,"%s %lf",&tag[0],&tres);		//time resolution in plasma/gyro-periods
fscanf(sinit,"%s %lf",&tag[0],&hist_t);	//time step for averaged diagnostics
fscanf(sinit,"%s %lf",&tag[0],&save_t);	//time step for field diagnostics
fscanf(sinit,"%s %lf",&tag[0],&pdf_t);		//time step for PDF diagnostics
fscanf(sinit,"%s %lf",&tag[0],&c);			//speed of light
fscanf(sinit,"%s %lf",&tag[0],&mratio);	//mp/me ratio
fscanf(sinit,"%s %lf",&tag[0],&beta);		//plasma beta
fscanf(sinit,"%s %lf",&tag[0],&alpha);		//angle between B and x axis in x-y plane
fscanf(sinit,"%s %lf",&tag[0],&delta);		//angle between B and y axis in y-z plane
fscanf(sinit,"%s %lf",&tag[0],&psi0);		//noise amplitude
fscanf(sinit,"%s %i",&tag[0],&trigger);	//noise type
fscanf(sinit,"%s %lf",&tag[0],&dummy1);		//reserve variable 1
fscanf(sinit,"%s %lf",&tag[0],&dummy2);		//reserve variable 2
fscanf(sinit,"%s %lf",&tag[0],&dummy3);		//reserve variable 3
//scan the parameters of the 1st particle population
fscanf(sinit,"%s %lf",&tag[0],&n1);	//density
fscanf(sinit,"%s %lf",&tag[0],&ud1);	//drift speed
fscanf(sinit,"%s %lf",&tag[0],&T1);	//(parallel) temeperature
fscanf(sinit,"%s %lf",&tag[0],&A1);	// temperature anisotropy per/par
fscanf(sinit,"%s %lf",&tag[0],&M1);	//mass
fscanf(sinit,"%s %lf",&tag[0],&Q1);	//charge
fscanf(sinit,"%s %lf",&tag[0],&k1);	//kappa index, kappa >1.5
fscanf(sinit,"%s %lf",&tag[0],&vmax1);	//v_max in velocity space
fscanf(sinit,"%s %i",&tag[0],&nv1);	//number of grids per [0, v_max]
#if defined(_2POP) || defined(_3POP) || defined(_4POP)
//scan the parameters of the 2nd particle population
fscanf(sinit,"%s %lf",&tag[0],&n2);	//density
fscanf(sinit,"%s %lf",&tag[0],&ud2);	//drift speed
fscanf(sinit,"%s %lf",&tag[0],&T2);	//(parallel) temeperature
fscanf(sinit,"%s %lf",&tag[0],&A2);	// temperature anisotropy per/par
fscanf(sinit,"%s %lf",&tag[0],&M2);	//mass
fscanf(sinit,"%s %lf",&tag[0],&Q2);	//charge
fscanf(sinit,"%s %lf",&tag[0],&k2);	//kappa index, kappa >1.5
fscanf(sinit,"%s %lf",&tag[0],&vmax2);	//v_max in velocity space
fscanf(sinit,"%s %i",&tag[0],&nv2);	//number of grids per [0, v_max]
#endif
#if defined(_3POP) || defined(_4POP)
//scan the parameters of the 3rd particle population
fscanf(sinit,"%s %lf",&tag[0],&n3);	//density
fscanf(sinit,"%s %lf",&tag[0],&ud3);	//drift speed
fscanf(sinit,"%s %lf",&tag[0],&T3);	//(parallel) temeperature
fscanf(sinit,"%s %lf",&tag[0],&A3);	// temperature anisotropy per/par
fscanf(sinit,"%s %lf",&tag[0],&M3);	//mass
fscanf(sinit,"%s %lf",&tag[0],&Q3);	//charge
fscanf(sinit,"%s %lf",&tag[0],&k3);	//kappa index, kappa >1.5
fscanf(sinit,"%s %lf",&tag[0],&vmax3);	//v_max in velocity space
fscanf(sinit,"%s %i",&tag[0],&nv3);	//number of grids per [0, v_max]
#endif
#if defined(_4POP)
//scan the parameters of the 4th particle population
fscanf(sinit,"%s %lf",&tag[0],&n4);	//density
fscanf(sinit,"%s %lf",&tag[0],&ud4);	//drift speed
fscanf(sinit,"%s %lf",&tag[0],&T4);	//(parallel) temeperature
fscanf(sinit,"%s %lf",&tag[0],&A4);	// temperature anisotropy per/par
fscanf(sinit,"%s %lf",&tag[0],&M4);	//mass
fscanf(sinit,"%s %lf",&tag[0],&Q4);	//charge
fscanf(sinit,"%s %lf",&tag[0],&k4);	//kappa index, kappa >1.5
fscanf(sinit,"%s %lf",&tag[0],&vmax4);	//v_max in velocity space
fscanf(sinit,"%s %i",&tag[0],&nv4);	//number of grids per [0, v_max]
#endif
fclose(sinit);
////////////////////////////////////////////////////////////////////////////////
////////////////////////// PHYSICAL VALUES CALCULATION /////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifdef _1D
ny=1; nz=1;
#endif
#ifdef _2D
nz=1;
#endif
ge=5./3.; c2=c*c; nr=nx*ny*nz; nw=nx/2-4; nwy=ny/4; nwz=nz/4;
alpha=alpha*TPI/180.; delta=delta*TPI/180.;

if(nx<16){printf("\n ERROR: the number of grids in the x direction must be no less than 16! \n ");
exit(1);}
if(n1==0){printf("\n ERROR: the density of the first particle population cannot be 0! \n ");
exit(1);}
if(Q1==0){printf("\n ERROR: the charge of the first particle population cannot be 0! \n ");
exit(1);}
if(M1==0){printf("\n ERROR: the mass of the first particle population cannot be 0! \n ");
exit(1);}
QM1=Q1/M1; nv1y=nv1z=2*nv1-1; nv12=nv1y*nv1z; nv13=nv1*nv12; n11=nr*nv13; n12=nr*nv12*nv1y;

#if defined(_NONINERT) && defined(_1POP)
printf("\n ERROR: the code was compiled with flag -D_NONINERT to enable moving non-inertial reference frame, \n this option needs at least 2 particle populations! \n ");
exit(1);
#endif

#if defined(_NONINERT) && defined(_2POP)
if(ud2==0){printf("\n ERROR: the code was compiled with flag -D_NONINERT to enable moving non-inertial reference frame, \n however drift velocity ud2 is 0! \n ");
	exit(1);}
#endif

#if defined(_NONINERT) && defined(_3POP)
if(ud3==0){printf("\n ERROR: the code was compiled with flag -D_NONINERT to enable moving non-inertial reference frame, \n however drift velocity ud3 is 0! \n ");
	exit(1);}
#endif

#if defined(_HYBRID)
if(Q1<0. || fabs(Q1)>1.){printf("\n ERROR: the code was compiled with flag -D_HYBRID for hybrid simulation, \n so the first particle population must be protons, Q1=%f! \n ", Q1);
	exit(1);}
#endif

#if !defined(_HYBRID) && !defined(_ELEC)
#if defined (_2POP)
if(Q1*n1+Q2*n2!=0){printf("\n ERROR: the charge neutrality is violated sum Qi*ni must be 0! \n "); exit(1);}
#elif defined (_3POP)
if(Q1*n1+Q2*n2+Q3*n3!=0){printf("\n ERROR: the charge neutrality is violated sum Qi*ni must be 0! \n "); exit(1);}
#endif
#endif

#if defined(_3POP)
if(n2==0||n3==0){printf("\n ERROR: the code was compiled with flag -D_3POP to enable three populations, \n however one of the three population densities is 0! \n ");
	exit(1);}
QM3=Q3/M3; QM2=Q2/M2;
nv2y=nv2z=2*nv2-1; nv22=nv2y*nv2z; nv23=nv2*nv22; n21=nr*nv23; n22=nr*nv22*nv2y;
nv3y=nv3z=2*nv3-1; nv32=nv3y*nv3z; nv33=nv3*nv32; n31=nr*nv33; n32=nr*nv32*nv3y;
#elif defined(_2POP)
if(n2==0){printf("\n ERROR: the code was compiled with flag -D_2POP to enable the second population, \n however it's density is 0! \n ");
	exit(1);}
QM2=Q2/M2; nv2y=nv2z=2*nv2-1; nv22=nv2y*nv2z; nv23=nv2*nv22; n21=nr*nv23; n22=nr*nv22*nv2y;
#endif



 
#if defined(_HYBRID)
//normalization for hybrid codes, mu0=1, epsilon0=1/c2
	wp1=wpp=c; wpe=wpp*sqrt(mratio); 
	T1=.5*beta; B0=sqrt(n1); wc1=wcp=B0; vA=B0; wce=wcp*mratio;
	vT1=vTp=sqrt(T1/M1); dV1=vmax1/nv1/vT1;L=vA/wcp; lDp=vTp/wpp;

	#if defined(_2POP) || defined(_3POP)
	T2=T2*T1; vT2=sqrt(T2/M2); dV2=vmax2/nv2/vT2;
	#endif
	#if defined(_3POP)
	T3=T3*T1; vT3=sqrt(T3/M3); dV3=vmax3/nv3/vT3;
	#endif
	
	
	
	
#else
//normalization for the electron and general codes, epsilon0=1, mu0=1/c2
//B0=sqrt(2.*n1*T1/beta)/c; //in SI units, if epsilon_0=1, mu_0 = 1/c^2
//vA=v_thi/sqrt(.5*betai); wci=wpi*vA/c; wce=wci*mratio;
#endif



Bx0=B0*cos(alpha); By0=B0*sin(alpha)*cos(delta); Bz0=B0*sin(alpha)*sin(delta);
xmax=xmax*L; dx=xmax/nx;  

coef=-2.*I*dx; coef1=-dx*I/TPI; //integration
coef2=TPI*I/(nx*nx*dx); //differentiation
coef3=-nx*dx*dx/(TPI*TPI); //second order integration
//coef3=-nx*nx*nx*dx*dx/(4.*TPI*TPI*TPI); //second order integration
coef4=-TPI*TPI/(nx*nx*nx*dx*dx); //second order differentiation

#if defined(_HYBRID)
dt=tres/wcp;  //time is in terms of proton gyrofrequency!!!
//printf("\n dt=%f  tres=%f   /  wcp = %f \n", dt, tres, wcp);
if(dx<c*dt){dt=.5*dx/c;} //CFL condition in config space
#else
dt=tres/wpe;  //time is in terms of electron plasma frequency!!!
if(dx<c*dt){dt=.5*dx/c;} //CFL condition in config space
#endif

//if(dVe<fabs(QMe)*dt*vmaxeinvt*B0){dt=.25*fabs(QMe)*B0/nve;} //CFL condition in phase space

//if(dt>.01/wci){dt=.01/wci;} //CFL condition for cyclotron rotation

dt1=.5*dt; dt2=dt/6.; dt3=dt/3.;

//printf("\n nx=%i \t ny=%i \t nz=%i", nx, ny, nz);
//printf("\n nve=%i \t nvi=%i \t nr=%i ", nve, nvi, nr);
//printf("\n nvey=%i \t nviy=%i \t nvez=%i \t nviz=%i", nvey, nviy, nvez, nviz);
//printf("\n ne=%i \t ni=%i \t ne2=%i \t ni2=%i ", ne, ni, ne2, ni2);

//printf("\n dt=%1.3e \t dVi=%1.3e \t dVe=%1.3e", dt, dVi, dVe);
//printf("\n Bx0=%1.3e \t By0=%1.3e \t Bz0=%1.3e", Bx0, By0, Bz0);

savestep=(int)(save_t/dt/wcp); if(savestep < 1){savestep=1;}
endstep=(int)(end_t/dt/wcp); if(endstep < 2){endstep=2;}
histstep=(int)(hist_t/dt/wcp); if(histstep < 1){histstep=1;}
pdfstep=(int)(pdf_t/dt/wcp); if(pdfstep < 1){pdfstep=1;}

////////////////////////////////////////////////////////////////////////////////
////////////////////////// DYNAMIC MEMORY ALLOCATION ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
//generic and auxilary variables allocated first
#if defined(_2D) || defined(_3D)
cx= (double *) calloc(nx, sizeof(double));
cy= (double *) calloc(ny, sizeof(double));
#endif
#if defined(_3D)
cz= (double *) calloc(nz, sizeof(double));
#endif
Ex= (double *) calloc(nr, sizeof(double));
Ey= (double *) calloc(nr, sizeof(double));
Ez= (double *) calloc(nr, sizeof(double));
Bx= (double *) calloc(nr, sizeof(double));
By= (double *) calloc(nr, sizeof(double));
Bz= (double *) calloc(nr, sizeof(double));
rho= (double *) calloc(nr, sizeof(double));
jx= (double *) calloc(nr, sizeof(double));
jy= (double *) calloc(nr, sizeof(double));
jz= (double *) calloc(nr, sizeof(double));
fluc= (double *) calloc(nr, sizeof(double));
flucx= (double *) calloc(nr, sizeof(double));
flucy= (double *) calloc(nr, sizeof(double));
flucz= (double *) calloc(nr, sizeof(double));
Ax= (double *) calloc(nr, sizeof(double));
Ay= (double *) calloc(nr, sizeof(double));
Az= (double *) calloc(nr, sizeof(double));
Axt= (double *) calloc(nr, sizeof(double));
Ayt= (double *) calloc(nr, sizeof(double));
Azt= (double *) calloc(nr, sizeof(double));
Axn= (double *) calloc(nr, sizeof(double));
Ayn= (double *) calloc(nr, sizeof(double));
Azn= (double *) calloc(nr, sizeof(double));
Ax0= (double *) calloc(nr, sizeof(double));
Ay0= (double *) calloc(nr, sizeof(double));
Az0= (double *) calloc(nr, sizeof(double));
dAxdy= (double *) calloc(nr, sizeof(double));
dAxdz= (double *) calloc(nr, sizeof(double));
dAydz= (double *) calloc(nr, sizeof(double));
d2Ax= (double *) calloc(nr, sizeof(double));
d2Ay= (double *) calloc(nr, sizeof(double));
d2Az= (double *) calloc(nr, sizeof(double));
dGx= (double *) calloc(nr, sizeof(double));
dGy= (double *) calloc(nr, sizeof(double));
dGz= (double *) calloc(nr, sizeof(double));
Gx= (double *) calloc(nr, sizeof(double));
Gy= (double *) calloc(nr, sizeof(double));
Gz= (double *) calloc(nr, sizeof(double));
Gxt= (double *) calloc(nr, sizeof(double));
Gyt= (double *) calloc(nr, sizeof(double));
Gzt= (double *) calloc(nr, sizeof(double));
Gxn= (double *) calloc(nr, sizeof(double));
Gyn= (double *) calloc(nr, sizeof(double));
Gzn= (double *) calloc(nr, sizeof(double));
Gx0= (double *) calloc(nr, sizeof(double));
Gy0= (double *) calloc(nr, sizeof(double));
Gz0= (double *) calloc(nr, sizeof(double));
#if !defined(_HYBRID)
phi= (double *) calloc(nr, sizeof(double));
dphi= (double *) calloc(nr, sizeof(double));
#else
rhoe= (double *) calloc(nr, sizeof(double));
uex= (double *) calloc(nr, sizeof(double));
pe0= (double *) calloc(nr, sizeof(double));
pe= (double *) calloc(nr, sizeof(double));
dpe= (double *) calloc(nr, sizeof(double));
pet= (double *) calloc(nr, sizeof(double));
pen= (double *) calloc(nr, sizeof(double));
ubx= (double *) calloc(nr, sizeof(double));
dubx= (double *) calloc(nr, sizeof(double));
#endif

in= (double *) calloc(nx, sizeof(double));
out= (double *) calloc(nx, sizeof(double));


//all particle specific variables allocated now
Vx1= (double *) calloc(nv1, sizeof(double));
Vy1= (double *) calloc(nv1y, sizeof(double));
Vz1= (double *) calloc(nv1z, sizeof(double));
rho1= (double *) calloc(nr, sizeof(double));
jx1= (double *) calloc(nr, sizeof(double));
jy1= (double *) calloc(nr, sizeof(double));
jz1= (double *) calloc(nr, sizeof(double));
t1x= (double *) calloc(nr, sizeof(double));
t1y= (double *) calloc(nr, sizeof(double));
t1z= (double *) calloc(nr, sizeof(double));


F1=(double complex*)fftw_malloc(sizeof(double complex)*n11);
Fn1=(double complex*)fftw_malloc(sizeof(double complex)*n11);
Ft1=(double complex*)fftw_malloc(sizeof(double complex)*n11);
rhok=(double complex*)fftw_malloc(sizeof(double complex)*nr);
jxk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
jyk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
jzk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
df2=(double complex*)fftw_malloc(sizeof(double complex)*nx);
kx=(double complex*)fftw_malloc(sizeof(double complex)*nx);
ky=(double complex*)fftw_malloc(sizeof(double complex)*ny);
kz=(double complex*)fftw_malloc(sizeof(double complex)*nz);
fbk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
axk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
ayk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
azk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
daxk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
dayk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
dazk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
ink=(double complex*)fftw_malloc(sizeof(double complex)*nx);
outk=(double complex*)fftw_malloc(sizeof(double complex)*nx);
fink=(double complex*)fftw_malloc(sizeof(double complex)*nr);
foutk=(double complex*)fftw_malloc(sizeof(double complex)*nr);

#if defined(_AB4)
//allocation of arrays for Adams-Bashforth scheme
F0=(double complex*)fftw_malloc(sizeof(double complex)*n11);
dF1=(double complex*)fftw_malloc(sizeof(double complex)*n11);
dF2=(double complex*)fftw_malloc(sizeof(double complex)*n11);
dF3=(double complex*)fftw_malloc(sizeof(double complex)*n11);
dF4=(double complex*)fftw_malloc(sizeof(double complex)*n11);
dAy1= (double *) calloc(nr, sizeof(double));
dAy2= (double *) calloc(nr, sizeof(double));
dAy3= (double *) calloc(nr, sizeof(double));
dAy4= (double *) calloc(nr, sizeof(double));
dAz1= (double *) calloc(nr, sizeof(double));
dAz2= (double *) calloc(nr, sizeof(double));
dAz3= (double *) calloc(nr, sizeof(double));
dAz4= (double *) calloc(nr, sizeof(double));
	#if defined(_2POP) || defined(_3POP)
	F20=(double complex*)fftw_malloc(sizeof(double complex)*n21);
	dF21=(double complex*)fftw_malloc(sizeof(double complex)*n21);
	dF22=(double complex*)fftw_malloc(sizeof(double complex)*n21);
	dF23=(double complex*)fftw_malloc(sizeof(double complex)*n21);
	dF24=(double complex*)fftw_malloc(sizeof(double complex)*n21);
	#endif
	#if defined(_3POP)
	F30=(double complex*)fftw_malloc(sizeof(double complex)*n31);
	dF31=(double complex*)fftw_malloc(sizeof(double complex)*n31);
	dF32=(double complex*)fftw_malloc(sizeof(double complex)*n31);
	dF33=(double complex*)fftw_malloc(sizeof(double complex)*n31);
	dF34=(double complex*)fftw_malloc(sizeof(double complex)*n31);
	#endif
	#if defined(_NONINERT)
	uy0= (double *) calloc(nr, sizeof(double));
	uz0= (double *) calloc(nr, sizeof(double));
	duy1= (double *) calloc(nr, sizeof(double));
	duy2= (double *) calloc(nr, sizeof(double));
	duy3= (double *) calloc(nr, sizeof(double));
	duy4= (double *) calloc(nr, sizeof(double));
	duz1= (double *) calloc(nr, sizeof(double));
	duz2= (double *) calloc(nr, sizeof(double));
	duz3= (double *) calloc(nr, sizeof(double));
	duz4= (double *) calloc(nr, sizeof(double));
	#endif
#endif

nvmax=nv1; //determine the size of auxilary arrays
#if defined(_2POP) || defined(_3POP)
if(nv2>nvmax){nvmax=nv2;} //if nv2 is larger, we'll use it
F2=(double complex*)fftw_malloc(sizeof(double complex)*n21);
Fn2=(double complex*)fftw_malloc(sizeof(double complex)*n21);
Ft2=(double complex*)fftw_malloc(sizeof(double complex)*n21);
G0=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv2*nv2y);
G1=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv2*nv2y);
G2=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv2*nv2y);
rho2= (double *) calloc(nx, sizeof(double));
jx2= (double *) calloc(nx, sizeof(double));
jy2= (double *) calloc(nx, sizeof(double));
jz2= (double *) calloc(nx, sizeof(double));
t2x= (double *) calloc(nr, sizeof(double));
t2y= (double *) calloc(nr, sizeof(double));
t2z= (double *) calloc(nr, sizeof(double));
Vx2= (double *) calloc(nv2, sizeof(double));
Vy2= (double *) calloc(nv2y, sizeof(double));
Vz2= (double *) calloc(nv2z, sizeof(double));
flucy2= (double *) calloc(nr, sizeof(double));
flucz2= (double *) calloc(nr, sizeof(double));
#if defined(_NONINERT)
uy= (double *) calloc(nx, sizeof(double));
uz= (double *) calloc(nx, sizeof(double));
uyave= (double *) calloc(nx, sizeof(double));
uzave= (double *) calloc(nx, sizeof(double));
uyold= (double *) calloc(nx, sizeof(double));
uzold= (double *) calloc(nx, sizeof(double));
uynew= (double *) calloc(nx, sizeof(double));
uznew= (double *) calloc(nx, sizeof(double));
uyt= (double *) calloc(nx, sizeof(double));
uzt= (double *) calloc(nx, sizeof(double));
duy= (double *) calloc(nx, sizeof(double));
duz= (double *) calloc(nx, sizeof(double));
#endif
#endif
#if defined(_3POP)
if(nv3>nvmax){nvmax=nv3;} //if nv3 is larger, we'll use it
F3=(double complex*)fftw_malloc(sizeof(double complex)*n31);
Fn3=(double complex*)fftw_malloc(sizeof(double complex)*n31);
Ft3=(double complex*)fftw_malloc(sizeof(double complex)*n31);
G0=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv3*nv3y);
G1=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv3*nv3y);
G2=(double complex*)fftw_malloc(sizeof(double complex)*nx*nv3*nv3y);
rho3= (double *) calloc(nx, sizeof(double));
jx3= (double *) calloc(nx, sizeof(double));
jy3= (double *) calloc(nx, sizeof(double));
jz3= (double *) calloc(nx, sizeof(double));
t3x= (double *) calloc(nr, sizeof(double));
t3y= (double *) calloc(nr, sizeof(double));
t3z= (double *) calloc(nr, sizeof(double));
Vx3= (double *) calloc(nv3, sizeof(double));
Vy3= (double *) calloc(nv3y, sizeof(double));
Vz3= (double *) calloc(nv3z, sizeof(double));
flucy3= (double *) calloc(nr, sizeof(double));
flucz3= (double *) calloc(nr, sizeof(double));
#endif
nvymax=2*nvmax-1; nvmax2=nvymax*nvymax;
bp=(double complex*)fftw_malloc(sizeof(double complex)*nr*nvmax2);
bm=(double complex*)fftw_malloc(sizeof(double complex)*nr*nvmax2);
f1=(double complex*)fftw_malloc(sizeof(double complex)*nx*nvmax2);
f2=(double complex*)fftw_malloc(sizeof(double complex)*nx*nvmax2);
kx1= (double *) calloc(nvmax2, sizeof(double));
fdxf= (fftw_plan *) calloc(nvmax2, sizeof(fftw_plan));
fdxb= (fftw_plan *) calloc(nvmax2, sizeof(fftw_plan));
G=(double complex*)fftw_malloc(sizeof(double complex)*nr*nvymax*nvmax2);
dfdv=(double complex*)fftw_malloc(sizeof(double complex)*nr*nvymax*nvmax2);
dft=(double *)fftw_malloc(sizeof(double)*nr*nvymax*nvmax2);
gam=(double *)fftw_malloc(sizeof(double)*nr*nvymax*nvmax2);
dfdvy=(double complex*)fftw_malloc(sizeof(double complex)*nr*nvmax*nvmax2);


transjxk=fftw_plan_dft_r2c_1d(nx,jx,jxk,FFTW_MEASURE);
transjyk=fftw_plan_dft_r2c_1d(nx,jy,jyk,FFTW_MEASURE);
transjzk=fftw_plan_dft_r2c_1d(nx,jz,jzk,FFTW_MEASURE);

transay0=fftw_plan_dft_c2r_1d(nx,ayk,Ay,FFTW_MEASURE);
transaz0=fftw_plan_dft_c2r_1d(nx,azk,Az,FFTW_MEASURE);

/* this entire block seems obsolete
transaxk=fftw_plan_dft_r2c_1d(nx,Ax,axk,FFTW_MEASURE);
transax=fftw_plan_dft_c2r_1d(nx,daxk,d2Ax,FFTW_MEASURE);

transayk=fftw_plan_dft_r2c_1d(nx,Ay,ayk,FFTW_MEASURE);
transay=fftw_plan_dft_c2r_1d(nx,dayk,d2Ay,FFTW_MEASURE);

transazk=fftw_plan_dft_r2c_1d(nx,Az,azk,FFTW_MEASURE);
transaz=fftw_plan_dft_c2r_1d(nx,dazk,d2Az,FFTW_MEASURE);

transrho=fftw_plan_dft_r2c_3d(nx,ny,nz,rho,rhok,FFTW_MEASURE);

*/
dxf=fftw_plan_dft_1d(nx,df2,fbk,FFTW_FORWARD,FFTW_MEASURE);
dxb=fftw_plan_dft_1d(nx,df2,fbk,FFTW_BACKWARD,FFTW_MEASURE);

dyf=fftw_plan_dft_1d(ny,df2,fbk,FFTW_FORWARD,FFTW_MEASURE);
dyb=fftw_plan_dft_1d(ny,df2,fbk,FFTW_BACKWARD,FFTW_MEASURE);

dzf=fftw_plan_dft_1d(nz,df2,fbk,FFTW_FORWARD,FFTW_MEASURE);
dzb=fftw_plan_dft_1d(nz,df2,fbk,FFTW_BACKWARD,FFTW_MEASURE);

db1=fftw_plan_dft_1d(nx,df2,fbk,FFTW_FORWARD,FFTW_MEASURE);
db2=fftw_plan_dft_1d(nx,fbk,df2,FFTW_BACKWARD,FFTW_MEASURE);


dfxin=fftw_plan_dft_r2c_1d(nx,in,ink,FFTW_MEASURE);
dfxout=fftw_plan_dft_c2r_1d(nx,outk,out,FFTW_MEASURE);

dfyin=fftw_plan_dft_r2c_1d(ny,in,ink,FFTW_MEASURE);
dfyout=fftw_plan_dft_c2r_1d(ny,outk,out,FFTW_MEASURE);

dfzin=fftw_plan_dft_r2c_1d(nz,in,ink,FFTW_MEASURE);
dfzout=fftw_plan_dft_c2r_1d(nz,outk,out,FFTW_MEASURE);


printf("\n memory allocated");

for(i=0;i<nvmax2;i++){
	fdxf[i]=fftw_plan_dft_1d(nx,&f1[nx*i],&f2[nx*i],FFTW_FORWARD,FFTW_MEASURE);
	fdxb[i]=fftw_plan_dft_1d(nx,&f1[nx*i],&f2[nx*i],FFTW_BACKWARD,FFTW_MEASURE);
}


for(ir=0;ir<nr;ir++){Bx[ir]=Bx0; By[ir]=By0; Bz[ir]=Bz0;} //initialize the external magnetic field

for(ix=0;ix<nx/2;ix++){kx[ix]=TPI/(nx*nx*dx)*ix; kx[nx-1-ix]=-TPI/(nx*nx*dx)*(ix+1);} 
for(ix=0;ix<3;ix++){kx[nx/2+ix]=0.; kx[nx/2-ix]=0.;} //filter out the harmonics close to Nyquist

#if defined(_2D) || defined(_3D)
for(iy=0;iy<ny/2;iy++){ky[iy]=TPI*I/(ny*ny*dx)*iy; ky[ny-1-iy]=-TPI*I/(ny*ny*dx)*(iy+1);} 
for(iy=0;iy<3;iy++){ky[ny/2-1+iy]=0.; ky[ny/2-1-iy]=0.;} 
for(i=0;i<nx;i++){cx[i]=cos(TPI*i/nx);}
for(i=0;i<ny;i++){cy[i]=cos(TPI*i/ny);}
#endif
#if defined(_3D)
for(iz=0;iz<nz/2;iz++){kz[iz]=TPI*I/(nz*nz*dx)*iz; kz[nz-1-iz]=-TPI*I/(nz*nz*dx)*(iz+1);} 
for(iz=0;iz<3;iz++){kz[nz/2-1+iz]=0.; kz[nz/2-1-iz]=0.;} 
for(i=0;i<nz;i++){cz[i]=cos(TPI*i/nz);}
#endif

/*printf("\n");
for(i=0;i<nx;i++){printf("kx[%i]=%1.3e \t",i,creal(kx[i]));}
printf("\n");
for(i=0;i<ny;i++){printf("ky[%i]=%1.3e \t",i,creal(ky[i]));}
printf("\n");
for(i=0;i<nz;i++){printf("kz[%i]=%1.3e \t",i,creal(kz[i]));}
*/

for(i=0;i<nv1;i++){Vx1[i]=dV1*i;} //  printf("\n Vx1[%i]=%f", i, Vxe[i]);}
for(i=0;i<nv1y;i++){Vy1[i]=Vz1[i]=dV1*(i-nv1+1);} // printf("\n Vy1[%i]=%f", i, Vy1[i]);}
#if defined(_2POP) || defined(_3POP)
for(i=0;i<nv2;i++){Vx2[i]=dV2*i;}
for(i=0;i<nv2y;i++){Vy2[i]=Vz2[i]=dV2*(i-nv2+1);}
/*if(setup==2){
	//printf("\n Bx0=%1.3f", Bx0);
	ptot=.5*Bx0*Bx0+(nb-n0)*(Ti/3.+2.*Ti*tper/3.)+n0*Tb;
	for(i=0;i<nx;i++){tix[i]=Ti; 
	if(i<nx/4 || i>3*nx/4){tiy[i]=tiz[i]=Ti;}
	if(i>=nx/4 && i<=3*nx/4){tiy[i]=tiz[i]=Ti*(2.-cos(4.*PI*(i-nx/4)/nx));}
	
	Bx[i]=sqrt(2.*(ptot-n0*Tb-(nb-n0)*(tix[i]+tiy[i]+tiz[i])/3.));
	//printf("\n Tx=%1.3f \t Ty=%1.3f \t Tz=%1.3f \t Bx=%1.3f", tix[i], tiy[i], tiz[i], Bx[i]);
	}
}*/
#endif
#if defined(_3POP)
for(i=0;i<nv3;i++){Vx3[i]=dV3*i;}
for(i=0;i<nv3y;i++){Vy3[i]=Vz3[i]=dV3*(i-nv3+1);}
/*if(setup==2){
	//printf("\n Bx0=%1.3f", Bx0);
	ptot=.5*Bx0*Bx0+no*To+(nb-n0-no)*(Ti/3.+2.*Ti*tper/3.)+n0*Tb;
	for(i=0;i<nx;i++){tix[i]=Ti; 
	if(i<nx/4 || i>3*nx/4){tiy[i]=tiz[i]=Ti;}
	if(i>=nx/4 && i<=3*nx/4){tiy[i]=tiz[i]=Ti*(2.-cos(4.*PI*(i-nx/4)/nx));}
	
	Bx[i]=sqrt(2.*(ptot-no*To-n0*Tb-(nb-n0-no)*(tix[i]+tiy[i]+tiz[i])/3.));
	//printf("\n Tx=%1.3f \t Ty=%1.3f \t Tz=%1.3f \t Bx=%1.3f", tix[i], tiy[i], tiz[i], Bx[i]);
	}
}*/
#endif

if((fopen("restart0.dat","r")) == NULL){ //if starting from scratch, not from the pre-saved data set

//dump=fopen("restart0.dat","wb"); //create the file for the data at the end of simulation

sprintf(tag,"restart0.dat"); //create the name for the data file

for(iout=0;iout<2;iout++){	//write the header in the stdout and in the param.dat file
    if(iout==1){streamw = fopen("param.dat","w");}
    if(iout==0){streamw=stdout;}

	#if defined(_HYBRID)
	sprintf(header,"Hybrid ");
	#endif
	#if defined(_ELEC)
	sprintf(header,"Electron ");
	#endif
	#if defined(_3D)
	strcat(header,"3D3V");
	#elif defined(_2D)
	strcat(header,"2D3V");
	#elif defined(_1D)
	strcat(header,"1D3V");
	#endif
	strcat(header," Fourier-Vlasov code vers. 27.06.2011"); //print out the date of the code version!!!
	
	
	#if defined(_OPENMP)
	fprintf(streamw,"\n %s with %i OMP threads\n", header,omp_get_max_threads());
	#else
	fprintf(streamw,"\n %s \n", header);
	#endif

#if defined(_AB4)
	fprintf(streamw,"\n 4th-order Adams-Bashforth time integration");
#else
	fprintf(streamw,"\n 4th-order Runge-Kutta time integration");
#endif
#if defined(_NONINERT)
	fprintf(streamw,"\n The beam reference frame is moving at variable speeds in all velocity directions");
#else
	fprintf(streamw,"\n The beam reference frame is stationary");
#endif
	fprintf(streamw,"\n Adaptive mesh refinement with tricubic interpolation is used");
	fprintf(streamw,"\n The number of resolved wave harmonics nw=%i",nw);
    fprintf(streamw,"\n xmax  = %.1f  \t  nx  = %i \t  ny = %i \t  nz = %i", xmax,nx,ny,nz);
    fprintf(streamw,"\n n1=%1.3f \t  ud1=%1.3f \t  T1=%.1f \t  M1=%.1f \t Q1=%.1f \t A1=%.1f \t nv1=%i",n1,ud1,T1,M1,Q1,A1,nv1);
	#if defined(_2POP) || defined(_3POP)
    fprintf(streamw,"\n n2=%1.3f \t  ud2=%1.3f \t  T2=%.1f \t  M2=%.1f \t Q2=%.1f \t A2=%.1f \t nv2=%i",n2,ud2,T2,M2,Q2,A2,nv2);
	#endif
	#if defined(_3POP)
    fprintf(streamw,"\n n3=%1.3f \t  ud3=%1.3f \t  T3=%.1f \t  M3=%.1f \t Q3=%.1f \t A3=%.1f \t nv3=%i",n3,ud3,T3,M3,Q3,A3,nv3);
	#endif
    fprintf(streamw,"\n\t mp/me = %1.2e \t\t\t  c = %1.2e",mratio, c);
    fprintf(streamw,"\n hist step = %.3f \t  fields saved   = %.3f  \t  end   = %.3f", hist_t,save_t,end_t);
    fprintf(streamw,"\n dt = %1.4f \t  PDF saved   = %.3f", dt, pdf_t);
	if(iout==1){fclose(streamw);}
}
////////////////////////////////////////////////////////////////////////////////
////////////////////// various types of initial perturbations //////////////////
////////////////////////////////////////////////////////////////////////////////
switch(trigger){ 
	case 0: //sinusoidal perturbation
	{for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++)
	for(iz=0;iz<nz;iz++){ir=ix*ny*nz+iy*nz+iz;
		//fluc[ir]=1.+psi0*sin(4.*TPI*ix/nx)*sin(2.*TPI*iy/ny)*sin(TPI*iz/nz);
		#if defined(_3D)
		fluc[ir]=1.+psi0*sin(TPI*ix/nx)*sin(TPI*iy/ny)*sin(TPI*iz/nz);
		flucx[ir]=TPI/dx/nx*psi0*cos(TPI*ix/nx)*sin(TPI*iy/ny)*sin(TPI*iz/nz);
		flucy[ir]=TPI/dx/ny*psi0*sin(TPI*ix/nx)*cos(TPI*iy/ny)*sin(TPI*iz/nz);
		flucz[ir]=TPI/dx/nz*psi0*sin(TPI*ix/nx)*sin(TPI*iy/ny)*cos(TPI*iz/nz);
		#elif defined(_1D)
		fluc[ir]=1.+psi0*sin(TPI*ix/nx);
		flucx[ir]=TPI/dx/nx*psi0*cos(TPI*ix/nx);
		#endif
		}} break;
	case 1:{            //noise
		th= (double *) calloc(nw, sizeof(double));
		for(k=0;k<2;k++){th[k]=Noise(TPI);}// printf("\n %1.3e", th[k]);}
		for(ix=0;ix<nx;ix++)
		for(iy=0;iy<ny;iy++)
		for(iz=0;iz<nz;iz++){ir=ix*ny*nz+iy*nz+iz;
			for(k=1;k<2;k++){
			fluc[ir]+=psi0*sin(k*TPI*ix/nx+th[k]);}}
		for(i=0;i<nr;i++){fluc[i]+=1.; flucy[i]=flucz[i]=0.;} //center around 1
	 } break;
	 case 2:{            //noise in perpendicular current
		for(i=0;i<nx;i++){flucy[i]=psi0*cos(TPI*i/nx); 
		flucz[i]=psi0*sin(TPI*i/nx);}
		for(i=0;i<nr;i++){fluc[i]=1.;} 
	 } break;
	 case 3:{            //noise in perpendicular current, waves going in both directions
		th= (double *) calloc(2*nw+1, sizeof(double));
		for(k=0;k<2*nw+1;k++){th[k]=Noise(TPI);}// printf("\n %1.3e", th[k]);}
#if defined(_3D)
		for(iy=0;iy<ny;iy++)
		for(iz=0;iz<nz;iz++){
		for(ix=0;ix<nx;ix++){i=ix*ny*nz+iy*nz+iz;
		for(k=1;k<nw+1;k++){
		flucy[i]+=psi0*cos(k*TPI*ix/nx+th[k]);flucz[i]+=psi0*sin(k*TPI*ix/nx+th[k]);}}	
		}
#elif defined(_1D)
		for(i=0;i<nx;i++){for(k=1;k<nw+1;k++){
		flucy[i]+=psi0*cos(k*TPI*i/nx+th[k]); flucz[i]+=psi0*sin(k*TPI*i/nx+th[k]);}}	
		for(i=0;i<nx;i++){for(k=nw+1;k<2*nw+1;k++){
		flucy[i]+=psi0*cos((k-nw-1)*TPI*i/nx+th[k]); flucz[i]-=psi0*sin((k-nw-1)*TPI*i/nx+th[k]);}}
	#if defined(_3POP)	
		for(k=0;k<2*nw+1;k++){th[k]=Noise(TPI);}
		for(i=0;i<nx;i++){for(k=1;k<nw+1;k++){
		flucy2[i]+=psi0*cos(k*TPI*i/nx+th[k]); flucz2[i]+=psi0*sin(k*TPI*i/nx+th[k]);}}	
		for(i=0;i<nx;i++){for(k=nw+1;k<2*nw+1;k++){
		flucy2[i]+=psi0*cos((k-nw-1)*TPI*i/nx+th[k]); flucz2[i]-=psi0*sin((k-nw-1)*TPI*i/nx+th[k]);}}
		for(k=0;k<2*nw+1;k++){th[k]=Noise(TPI);}
		for(i=0;i<nx;i++){for(k=1;k<nw+1;k++){
		flucy3[i]+=psi0*cos(k*TPI*i/nx+th[k]); flucz3[i]+=psi0*sin(k*TPI*i/nx+th[k]);}}	
		for(i=0;i<nx;i++){for(k=nw+1;k<2*nw+1;k++){
		flucy3[i]+=psi0*cos((k-nw-1)*TPI*i/nx+th[k]); flucz3[i]-=psi0*sin((k-nw-1)*TPI*i/nx+th[k]);}}
	#endif		
#endif		
		for(i=0;i<nr;i++){fluc[i]=1.;flucx[i]=0.;} 
	 } break;
	 case 4:{            //3D noise in density
		th= (double *) calloc(nw, sizeof(double));
		for(k=0;k<nw;k++){th[k]=Noise(TPI);}// printf("\n %1.3e", th[k]);}
		for(ix=0;ix<nx;ix++)
		for(iy=0;iy<ny;iy++)
		for(iz=0;iz<nz;iz++){ir=ix*ny*nz+iy*nz+iz;
			for(k=1;k<nw;k++){
			fluc[ir]+=psi0*sin(k*TPI*ix/nx+th[k])*sin(TPI*iy/ny)*sin(TPI*iz/nz);}}
		for(i=0;i<nr;i++){fluc[i]+=1.; flucy[i]=flucz[i]=0.;} //center around 1
	 }
}
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// PDF initialization /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if defined(_KAPPA)
    //kappa-distributions
	if(k1>1.5){norm=1./pow(.0000001,k1+.5)/gsl_sf_bessel_Knu(k1+.5,.0000001);
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv1;ivx++)
		for(ivy=0;ivy<nv1y;ivy++)	
		for(ivz=0;ivz<nv1z;ivz++){i=ir*nv13+ivx*nv12+ivy*nv1z+ivz;
			if(ivx==0 && ivy==nv1-1 && ivz==nv1-1){F1[i]=n1/TPI3;}
			else{
			F1[i]=n1*norm*pow(sqrt(2.*T1*(k1-1.5)
			*(Vx1[ivx]*Vx1[ivx]+(A1+1.)*(Vy1[ivy]*Vy1[ivy]+Vz1[ivz]*Vz1[ivz]))/M1),k1+.5)
			*gsl_sf_bessel_Knu(k1+.5,sqrt(2.*T1*(k1-1.5)*(Vx1[ivx]*Vx1[ivx]+(A1+1.)*(Vy1[ivy]*Vy1[ivy]+Vz1[ivz]*Vz1[ivz]))/M1))
			*(cos(ud1*Vx1[ivx])+I*sin(ud1*Vx1[ivx]))
			*(cos(flucy[ir]*Vy1[ivy])+I*sin(flucy[ir]*Vy1[ivy]))
			*(cos(flucz[ir]*Vz1[ivz])+I*sin(flucz[ir]*Vz1[ivz]))
			/TPI3;
			}
		}
	}
	}
	else{printf("\n WARNING: the code was compiled with flag -D_KAPPA for simulation of kappa-distributions, \n but kappa k1 <= 1.5, so Maxwellian is used! \n ");
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv1;ivx++)
		for(ivy=0;ivy<nv1y;ivy++)	
		for(ivz=0;ivz<nv1z;ivz++){i=ir*nv13+ivx*nv12+ivy*nv1z+ivz;
			F1[i]=n1*exp(-.5*T1*(Vx1[ivx]*Vx1[ivx]+(A1+1.)*(Vy1[ivy]*Vy1[ivy]+Vz1[ivz]*Vz1[ivz]))/M1)
			*(cos(ud1*Vx1[ivx])+I*sin(ud1*Vx1[ivx]))
			*(cos(flucy[ir]*Vy1[ivy])+I*sin(flucy[ir]*Vy1[ivy]))
			*(cos(flucz[ir]*Vz1[ivz])+I*sin(flucz[ir]*Vz1[ivz]))
			/TPI3;}
		}
	//printf("\n F1 initialization complete");
	}

	#if defined(_2POP) || defined(_3POP)
	if(k2>1.5){norm=1./pow(.0000001,k2+.5)/gsl_sf_bessel_Knu(k2+.5,.0000001);
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv2;ivx++)
		for(ivy=0;ivy<nv2y;ivy++)	
		for(ivz=0;ivz<nv2z;ivz++){i=ir*nv23+ivx*nv22+ivy*nv2z+ivz;
			if(ivx==0 && ivy==nv2-1 && ivz==nv2-1){F2[i]=n2/TPI3;}
			else{
			F2[i]=n2*norm*pow(sqrt(2.*T2*(k2-1.5)
			*(Vx2[ivx]*Vx2[ivx]+(A2+1.)*(Vy2[ivy]*Vy2[ivy]+Vz2[ivz]*Vz2[ivz]))/M2),k2+.5)
			*gsl_sf_bessel_Knu(k2+.5,sqrt(2.*T2*(k2-1.5)*(Vx2[ivx]*Vx2[ivx]+(A2+1.)*(Vy2[ivy]*Vy2[ivy]+Vz2[ivz]*Vz2[ivz]))/M2))
			*(cos(ud2*Vx2[ivx])+I*sin(ud2*Vx2[ivx]))
			*(cos(flucy2[ir]*Vy2[ivy])+I*sin(flucy2[ir]*Vy2[ivy]))
			*(cos(flucz2[ir]*Vz2[ivz])+I*sin(flucz2[ir]*Vz2[ivz]))
			/TPI3;
			}
		}
	}
	}
	else{printf("\n WARNING: the code was compiled with flag -D_KAPPA for simulation of kappa-distributions, \n but kappa k2 <= 1.5, so Maxwellian is used! \n ");
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv2;ivx++)
		for(ivy=0;ivy<nv2y;ivy++)	
		for(ivz=0;ivz<nv2z;ivz++){i=ir*nv23+ivx*nv22+ivy*nv2z+ivz;
			F2[i]=n2*exp(-.5*T2*(Vx2[ivx]*Vx2[ivx]+(A2+1.)*(Vy2[ivy]*Vy2[ivy]+Vz2[ivz]*Vz2[ivz]))/M2)
	#ifndef _NONINERT
			*(cos(ud2*Vx2[ivx])+I*sin(ud2*Vx2[ivx]))
	#endif
			*(cos(flucy2[ir]*Vy2[ivy])+I*sin(flucy2[ir]*Vy2[ivy]))
			*(cos(flucz2[ir]*Vz2[ivz])+I*sin(flucz2[ir]*Vz2[ivz]))
			/TPI3;}
		}
	}
	#endif
	#if defined(_3POP)
	if(k3>1.5){norm=1./pow(.0000001,k3+.5)/gsl_sf_bessel_Knu(k3+.5,.0000001);
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv3;ivx++)
		for(ivy=0;ivy<nv3y;ivy++)	
		for(ivz=0;ivz<nv3z;ivz++){i=ir*nv33+ivx*nv32+ivy*nv3z+ivz;
			if(ivx==0 && ivy==nv3-1 && ivz==nv3-1){F3[i]=n3/TPI3;}
			else{
			F3[i]=n3*norm*pow(sqrt(2.*T3*(k3-1.5)
			*(Vx3[ivx]*Vx3[ivx]+(A3+1.)*(Vy3[ivy]*Vy3[ivy]+Vz3[ivz]*Vz3[ivz]))/M3),k3+.5)
			*gsl_sf_bessel_Knu(k3+.5,sqrt(2.*T3*(k3-1.5)*(Vx3[ivx]*Vx3[ivx]+(A3+1.)*(Vy3[ivy]*Vy3[ivy]+Vz3[ivz]*Vz3[ivz]))/M3))
			*(cos(ud3*Vx3[ivx])+I*sin(ud3*Vx3[ivx]))
			*(cos(flucy3[ir]*Vy3[ivy])+I*sin(flucy3[ir]*Vy3[ivy]))
			*(cos(flucz3[ir]*Vz3[ivz])+I*sin(flucz3[ir]*Vz3[ivz]))
			/TPI3;
			}
		}
	}
	}
	else{printf("\n WARNING: the code was compiled with flag -D_KAPPA for simulation of kappa-distributions, \n but kappa k2 <= 1.5, so Maxwellian is used! \n ");
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv3;ivx++)
		for(ivy=0;ivy<nv3y;ivy++)	
		for(ivz=0;ivz<nv3z;ivz++){i=ir*nv33+ivx*nv32+ivy*nv3z+ivz;
			F3[i]=n3*exp(-.5*T3*(Vx3[ivx]*Vx3[ivx]+(A3+1.)*(Vy3[ivy]*Vy3[ivy]+Vz3[ivz]*Vz3[ivz]))/M3)
	#ifndef _NONINERT
			*(cos(ud3*Vx3[ivx])+I*sin(ud3*Vx3[ivx]))
	#endif
			*(cos(flucy3[ir]*Vy3[ivy])+I*sin(flucy3[ir]*Vy3[ivy]))
			*(cos(flucz3[ir]*Vz3[ivz])+I*sin(flucz3[ir]*Vz3[ivz]))
			/TPI3;}
		}
	}
	#endif

	
#else
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nv1;ivx++)
		for(ivy=0;ivy<nv1y;ivy++)	
		for(ivz=0;ivz<nv1z;ivz++){i=ir*nv13+ivx*nv12+ivy*nv1z+ivz;
			F1[i]=n1*exp(-.5*T1*(Vx1[ivx]*Vx1[ivx]+(A1+1.)*(Vy1[ivy]*Vy1[ivy]+Vz1[ivz]*Vz1[ivz]))/M1)
			*(cos(ud1*Vx1[ivx])+I*sin(ud1*Vx1[ivx]))
			*(cos(flucy[ir]*Vy1[ivy])+I*sin(flucy[ir]*Vy1[ivy]))
			*(cos(flucz[ir]*Vz1[ivz])+I*sin(flucz[ir]*Vz1[ivz]))
			/TPI3;}
#if defined(_2POP) || defined(_3POP)
		for(ivx=0;ivx<nv2;ivx++)
		for(ivy=0;ivy<nv2y;ivy++)	
		for(ivz=0;ivz<nv2z;ivz++){i=ir*nv23+ivx*nv22+ivy*nv2z+ivz;
			F2[i]=n2*exp(-.5*T2*(Vx2[ivx]*Vx2[ivx]+(A2+1.)*(Vy2[ivy]*Vy2[ivy]+Vz2[ivz]*Vz2[ivz]))/M2)
	#ifndef _NONINERT
			*(cos(ud2*Vx2[ivx])+I*sin(ud2*Vx2[ivx]))
	#endif
			*(cos(flucy2[ir]*Vy2[ivy])+I*sin(flucy2[ir]*Vy2[ivy]))
			*(cos(flucz2[ir]*Vz2[ivz])+I*sin(flucz2[ir]*Vz2[ivz]))
			/TPI3;}
#endif
#if defined(_3POP)	
		for(ivx=0;ivx<nv3;ivx++)
		for(ivy=0;ivy<nv3y;ivy++)	
		for(ivz=0;ivz<nv3z;ivz++){i=ir*nv33+ivx*nv32+ivy*nv3z+ivz;
			F3[i]=n3*exp(-.5*T3*(Vx3[ivx]*Vx3[ivx]+(A3+1.)*(Vy3[ivy]*Vy3[ivy]+Vz3[ivz]*Vz3[ivz]))/M3)
	#ifndef _NONINERT
			*(cos(ud3*Vx3[ivx])+I*sin(ud3*Vx3[ivx]))
	#endif
			*(cos(flucy3[ir]*Vy3[ivy])+I*sin(flucy3[ir]*Vy3[ivy]))
			*(cos(flucz3[ir]*Vz3[ivz])+I*sin(flucz3[ir]*Vz3[ivz]))
			/TPI3;}
#endif
	}



#endif

printf("\n PDFs initialized");

/*
	printf("\n localized anisotropy!!!");
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static) private(i,ir,ivx,ivy,ivz)
	#endif
	for(ir=0;ir<nr;ir++){
		for(ivx=0;ivx<nvi;ivx++)
		for(ivy=0;ivy<nviy;ivy++)	
		for(ivz=0;ivz<nviz;ivz++){i=ir*nvi3+ivx*nvi2+ivy*nviz+ivz;
#if defined(_3POP)	
			Fi[i]=(nb-n0-no)*exp(-.5*(tix[ir]*(Vxi[ivx]*Vxi[ivx])+tiy[ir]*(Vyi[ivy]*Vyi[ivy])+tiz[ir]*(Vzi[ivz]*Vzi[ivz]))/Mi)
			*(cos(flucy[ir]*Vyi[ivy])+I*sin(flucy[ir]*Vyi[ivy]))
			*(cos(flucz[ir]*Vzi[ivz])+I*sin(flucz[ir]*Vzi[ivz]))
			/TPI3;
		
			Fb[i]=n0*exp(-.5*Tb*(Vxe[ivx]*Vxe[ivx]+btper*(Vye[ivy]*Vye[ivy]+Vze[ivz]*Vze[ivz]))/Mi/mbratio)
			// *(cos(udi*Vxe[ivx])+I*sin(udi*Vxe[ivx]))
			*(cos(flucyb[ir]*Vye[ivy])+I*sin(flucyb[ir]*Vye[ivy]))
			*(cos(fluczb[ir]*Vze[ivz])+I*sin(fluczb[ir]*Vze[ivz]))
			/TPI3;

#elif defined(_2POP)
			Fi[i]=(nb-n0)*exp(-.5*(tix[ir]*(Vxi[ivx]*Vxi[ivx])+tiy[ir]*(Vyi[ivy]*Vyi[ivy])+tiz[ir]*(Vzi[ivz]*Vzi[ivz]))/Mi)
			*(cos(flucy[ir]*Vyi[ivy])+I*sin(flucy[ir]*Vyi[ivy]))
			*(cos(flucz[ir]*Vzi[ivz])+I*sin(flucz[ir]*Vzi[ivz]))
			/TPI3;


#endif
		}	
   }

*/

t=0.; step=0; printparam();
//if(2.*lD<dx){printf("\n WARNING: l_D not resolved... \n");}

save=fopen("timeave.dat","w"); fclose(save); saves = fopen("timeavesign.dat","w"); fclose(saves);
sEx=fopen("Ex.dat","wb"); sEy=fopen("Ey.dat","wb"); sEz=fopen("Ez.dat","wb"); 
sBy=fopen("By.dat","wb"); sBz=fopen("Bz.dat","wb"); //sBx=fopen("Bx.dat","wb"); 
//srho=fopen("rho.dat","wb");sjx=fopen("jx.dat","wb");sjy=fopen("jy.dat","wb");sjz=fopen("jz.dat","wb");
sf1=fopen("f1.dat","wb");  srho1=fopen("rho1.dat","wb"); 
sjx1=fopen("jx1.dat","wb");  sjy1=fopen("jy1.dat","wb");  sjz1=fopen("jz1.dat","wb");  

zeit(); runtime(0);
#if defined(_2POP) || defined(_3POP)
sf2=fopen("f2.dat","wb");  srho2=fopen("rho2.dat","wb"); 
sjx2=fopen("jx2.dat","wb");  sjy2=fopen("jy2.dat","wb");  sjz2=fopen("jz2.dat","wb");  
#if defined(_NONINERT)
suy=fopen("uy.dat","wb"); suz=fopen("uz.dat","wb"); 
udold=udave=udnew=ud2;
for(ix=0;ix<nx;ix++){
	uyave[ix]=uyold[ix]=uynew[ix]=0.; uzave[ix]=uzold[ix]=uznew[ix]=0.;}
#endif
#endif
#if defined(_3POP)
sf3=fopen("f3.dat","wb");  srho3=fopen("rho3.dat","wb"); 
sjx3=fopen("jx3.dat","wb");  sjy3=fopen("jy3.dat","wb");  sjz3=fopen("jz3.dat","wb");  
moments(F1, F2, F3); 
#if defined(_NONINERT)
udold=udave=udnew=ud3; 
#endif
#endif
#if defined(_2POP)
moments(F1, F2); 
#elif defined(_1POP)
moments(F1); 
#endif
fields0(); moments2(); histdata(); dumpdata(); dumpdist();
}
else{ //if starting from pre-saved data
	i=0;
	sprintf(tag,"restart0.dat");
	while((fopen(tag,"r")) != NULL){i++; //check for the last saved data set
	sprintf(tag,"restart%d.dat",i); } //i is the new index for the next check-point file
	dump=fopen(tag,"wb"); sprintf(tag,"restart%d.dat",i-1); restart=fopen(tag,"r");
	
zeit(); runtime(0);
reload(restart); //read data from the last check-point save
#if defined(_3POP)
moments(F1, F2, F3); 
#elif defined(_2POP)
moments(F1, F2); 
#else
moments(F1); 
#endif
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////// THIS IS THE MAIN PROGRAM CYCLE //////////////////////
////////////////////////////////////////////////////////////////////////////////
#if !defined(_AB4)
//The 4th order Runge-Kutta scheme
while(step<=endstep){
    t+=dt; step++; integrateRK4();	
	//printf("\n udold=%1.3e udave=%1.3e", udold, udave);
	//printf("\n udold=%1.3e udave=%1.3e udnew=%1.3e", udold, udave, udnew);
	//udold=udave; //udave=udnew;
	if(step%histstep==0){moments2(); histdata();}
    if(step%pdfstep==0){dumpdist();} //save the PDF
    if(step%savestep==0){dumpdata();} //save the fields
	
	if(ave(t1x)+ave(t1y)+ave(t1z)>9.*T1&&T1<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV1=.5*dV1; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv1;i++){Vx1[i]=dV1*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv1y;i++){Vy1[i]=Vz1[i]=dV1*(i-nv1+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F1, nv1); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T1 = %f before the refinement, T1 = %f, after the refinement", T1, 3.*T1);
		T1=T1*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#if defined(_2POP) || defined(_3POP)
	if(ave(t2x)+ave(t2y)+ave(t2z)>9.*T2&&T2<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV2=.5*dV2; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv2;i++){Vx2[i]=dV2*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv2y;i++){Vy2[i]=Vz2[i]=dV2*(i-nv2+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F2, nv2); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T2 = %f before the refinement, T2 = %f, after the refinement", T2, 3.*T2);
		T2=T2*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	#if defined(_3POP)
	if(ave(t3x)+ave(t3y)+ave(t3z)>9.*T3&&T3<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV3=.5*dV3; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv3;i++){Vx3[i]=dV3*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv3y;i++){Vy3[i]=Vz3[i]=dV3*(i-nv3+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F3, nv3); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T3 = %f before the refinement, T3 = %f, after the refinement", T3, 3.*T3);
		T3=T3*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	//printf("\n diagnostics for step %i complete in %1.3e sec", step, zeit());
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());

	#if defined(_1POP)
	if(fabs(ave(rhoe)-n1)>.01*n1){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_2POP)
	if(fabs(ave(rhoe)-n1-n2)>.01*(n1+n2)){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_3POP)
	if(fabs(ave(rhoe)-n1-n2-n3)>.01*(n1+n2+n3)){printf("\n density conservation violated! \n"); exit(1);}
	#endif
	
	//printf("\n total runtime to current step is %1.3e sec", runtime(1));
	#if defined(_OPENMP)
	if(runtime(1)>86100){ //86100 sec termination for 24-hour runs!!!
	zeit();
	printf("\n total runtime is %d sec, dumping the data for a restarting point \n", runtime(1));
	dump=fopen(tag,"wb"); dumpall(dump); //fclose(dump); 
	printf("\n dumping complete in %1.3e sec", zeit());
	exit(2);}
	#endif
}
#else
//the Adams-Bashforth scheme
	#if defined(_NONINERT)
		#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){F30[i]=F3[i];}	
	#endif
	ud0=udold;
	for(i=0;i<nr;i++){Ay0[i]=Ay[i]; Az0[i]=Az[i]; uy0[i]=uyold[i]; uz0[i]=uzold[i];}
	integrateRK4(); t+=dt; step++; 
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF24[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF24[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF34[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif
	dud4=ud-ud0; ud0=ud;
	for(i=0;i<nr;i++){dAy4[i]=Ay[i]-Ay0[i]; dAz4[i]=Az[i]-Az0[i]; 
	duy4[i]=uy[i]-uy0[i]; duz4[i]=uz[i]-uz0[i]; uy0[i]=uy[i]; uz0[i]=uz[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF23[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF23[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF33[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif
	dud3=ud-ud0; ud0=ud;
	for(i=0;i<nr;i++){dAy3[i]=Ay[i]-Ay0[i]; dAz3[i]=Az[i]-Az0[i];
	duy3[i]=uy[i]-uy0[i]; duz3[i]=uz[i]-uz0[i]; uy0[i]=uy[i]; uz0[i]=uz[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF22[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF22[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF32[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif
	dud2=ud-ud0; ud0=ud;
	for(i=0;i<nr;i++){dAy2[i]=Ay[i]-Ay0[i]; dAz2[i]=Az[i]-Az0[i];
	duy2[i]=uy[i]-uy0[i]; duz2[i]=uz[i]-uz0[i]; uy0[i]=uy[i]; uz0[i]=uz[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	//this is the first Adams-Bashforth time advancement

	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF21[i]=F2[i]-F20[i]; F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF21[i]=F2[i]-F20[i]; F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF31[i]=F3[i]-F30[i]; F3[i]+=B1*dF31[i]+B2*dF32[i]+B3*dF33[i]+B4*dF34[i];}
	#endif
	dud1=ud-ud0; udold=ud0+B1*dud1+B2*dud2+B3*dud3+B4*dud4;
	for(i=0;i<nr;i++){dAy1[i]=Ay[i]-Ay0[i]; dAz1[i]=Az[i]-Az0[i];  
		Ay[i]+=B1*dAy1[i]+B2*dAy2[i]+B3*dAy3[i]+B4*dAy4[i]; 
		Az[i]+=B1*dAz1[i]+B2*dAz2[i]+B3*dAz3[i]+B4*dAz4[i];//}
		duy1[i]=uy[i]-uy0[i]; 
		duz1[i]=uz[i]-uz0[i];
		uyold[i]+=B1*duy1[i]+B2*duy2[i]+B3*duy3[i]+B4*duy4[i]; 
		uzold[i]+=B1*duz1[i]+B2*duz2[i]+B3*duz3[i]+B4*duz4[i];}

	while(step<=endstep){
    t+=dt; step++; integrateAB4();

	if(step%histstep==0){moments2(); histdata();}
    if(step%pdfstep==0){dumpdist();} //save the PDF
    if(step%savestep==0){dumpdata();} //save the fields
	
	if(ave(t1x)+ave(t1y)+ave(t1z)>9.*T1&&T1<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV1=.5*dV1; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv1;i++){Vx1[i]=dV1*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv1y;i++){Vy1[i]=Vz1[i]=dV1*(i-nv1+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F1, nv1); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T1 = %f before the refinement, T1 = %f, after the refinement", T1, 3.*T1);
		T1=T1*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#if defined(_2POP) || defined(_3POP)
	if(ave(t2x)+ave(t2y)+ave(t2z)>9.*T2&&T2<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV2=.5*dV2; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv2;i++){Vx2[i]=dV2*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv2y;i++){Vy2[i]=Vz2[i]=dV2*(i-nv2+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F2, nv2); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T2 = %f before the refinement, T2 = %f, after the refinement", T2, 3.*T2);
		T2=T2*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	#if defined(_3POP)
	if(ave(t3x)+ave(t3y)+ave(t3z)>9.*T3&&T3<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV3=.5*dV3; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv3;i++){Vx3[i]=dV3*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv3y;i++){Vy3[i]=Vz3[i]=dV3*(i-nv3+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F3, nv3); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T3 = %f before the refinement, T3 = %f, after the refinement", T3, 3.*T3);
		T3=T3*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	printf("\n step %i AB4 complete in %1.3e sec", step, zeit());

	#if defined(_1POP)
	if(fabs(ave(rhoe)-n1)>.01*n1){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_2POP)
	if(fabs(ave(rhoe)-n1-n2)>.01*(n1+n2)){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_3POP)
	if(fabs(ave(rhoe)-n1-n2-n3)>.01*(n1+n2+n3)){printf("\n density conservation violated! \n"); exit(1);}
	#endif
	
	//printf("\n total runtime to current step is %1.3e sec", runtime(1));
	#if defined(_OPENMP)
	if(runtime(1)>86100){ //86100 sec termination for 24-hour runs!!!
	zeit();
	printf("\n total runtime is %d sec, dumping the data for a restarting point \n", runtime(1));
	dump=fopen(tag,"wb"); dumpall(dump); //fclose(dump); 
	printf("\n dumping complete in %1.3e sec", zeit());
	exit(2);}
	#endif
	}
	#else
	//Adams-Bashforth inegration with moving frame of reference
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){F0[i]=F1[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){F30[i]=F3[i];}	
	#endif
	for(i=0;i<nr;i++){Ay0[i]=Ay[i]; Az0[i]=Az[i];}
	integrateRK4(); t+=dt; step++; 
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF24[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF4[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF24[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF34[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif
	
	for(i=0;i<nr;i++){dAy4[i]=Ay[i]-Ay0[i]; dAz4[i]=Az[i]-Az0[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF23[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF3[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF23[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF33[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif

	for(i=0;i<nr;i++){dAy3[i]=Ay[i]-Ay0[i]; dAz3[i]=Az[i]-Az0[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	
	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF22[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF2[i]=F1[i]-F0[i]; F0[i]=F1[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF22[i]=F2[i]-F20[i]; F20[i]=F2[i];}	
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF32[i]=F3[i]-F30[i]; F30[i]=F3[i];}	
	#endif

	for(i=0;i<nr;i++){dAy2[i]=Ay[i]-Ay0[i]; dAz2[i]=Az[i]-Az0[i];}
	integrateRK4(); t+=dt; step++;
	if(step%histstep==0){moments2(); histdata();}
	if(step%pdfstep==0){dumpdist();} //save the PDF
	if(step%savestep==0){dumpdata();} //save the fields
	printf("\n step %i RK4 complete in %1.3e sec", step, zeit());
	//this is the first Adams-Bashforth time advancement

	#if defined(_1POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
	#elif defined(_2POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF21[i]=F2[i]-F20[i]; F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
	#elif defined(_3POP)
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n11;i++){dF1[i]=F1[i]-F0[i]; F1[i]+=B1*dF1[i]+B2*dF2[i]+B3*dF3[i]+B4*dF4[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n21;i++){dF21[i]=F2[i]-F20[i]; F2[i]+=B1*dF21[i]+B2*dF22[i]+B3*dF23[i]+B4*dF24[i];}
		#pragma omp parallel for schedule(static) private(i)
		for(i=0;i<n31;i++){dF31[i]=F3[i]-F30[i]; F3[i]+=B1*dF31[i]+B2*dF32[i]+B3*dF33[i]+B4*dF34[i];}
	#endif

	for(i=0;i<nr;i++){dAy1[i]=Ay[i]-Ay0[i]; dAz1[i]=Az[i]-Az0[i];  
		Ay[i]+=B1*dAy1[i]+B2*dAy2[i]+B3*dAy3[i]+B4*dAy4[i]; 
		Az[i]+=B1*dAz1[i]+B2*dAz2[i]+B3*dAz3[i]+B4*dAz4[i];}

	while(step<=endstep){
    t+=dt; step++; integrateAB4();

	if(step%histstep==0){moments2(); histdata();}
    if(step%pdfstep==0){dumpdist();} //save the PDF
    if(step%savestep==0){dumpdata();} //save the fields
	
	if(ave(t1x)+ave(t1y)+ave(t1z)>9.*T1&&T1<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV1=.5*dV1; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv1;i++){Vx1[i]=dV1*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv1y;i++){Vy1[i]=Vz1[i]=dV1*(i-nv1+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F1, nv1); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T1 = %f before the refinement, T1 = %f, after the refinement", T1, 3.*T1);
		T1=T1*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#if defined(_2POP) || defined(_3POP)
	if(ave(t2x)+ave(t2y)+ave(t2z)>9.*T2&&T2<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV2=.5*dV2; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv2;i++){Vx2[i]=dV2*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv2y;i++){Vy2[i]=Vz2[i]=dV2*(i-nv2+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F2, nv2); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T2 = %f before the refinement, T2 = %f, after the refinement", T2, 3.*T2);
		T2=T2*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	#if defined(_3POP)
	if(ave(t3x)+ave(t3y)+ave(t3z)>9.*T3&&T3<10.){ //if the average temperature increases by 3 
		dumpdist(); //save the old PDFs
		dV3=.5*dV3; //printf("\n dVe=%f",  dVe); //redefine the grid size
		for(i=0;i<nv3;i++){Vx3[i]=dV3*i;}//  printf("\n Vxe[%i]=%f", i, Vxe[i]);} //redefine the velocity grid
		for(i=0;i<nv3y;i++){Vy3[i]=Vz3[i]=dV3*(i-nv3+1);}// printf("\n Vye[%i]=%f", i, Vye[i]);}
		refine(F3, nv3); //refine the PDF on the new grid
		dumpdist(); //save the new PDF
		ref=fopen("refinement.dat","a"); //report the refinement 
		fprintf(ref,"\n Adaptive mesh refinement invoked at t = %f, step = %i", t, step);
		fprintf(ref,"\n T3 = %f before the refinement, T3 = %f, after the refinement", T3, 3.*T3);
		T3=T3*3.; //redefine the new average temperature
		fclose(ref); //release the handle of the report file
	}
	#endif
	
	printf("\n step %i AB4 complete in %1.3e sec", step, zeit());

	#if defined(_1POP)
	if(fabs(ave(rhoe)-n1)>.01*n1){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_2POP)
	if(fabs(ave(rhoe)-n1-n2)>.01*(n1+n2)){printf("\n density conservation violated! \n"); exit(1);}
	#elif defined(_3POP)
	if(fabs(ave(rhoe)-n1-n2-n3)>.01*(n1+n2+n3)){printf("\n density conservation violated! \n"); exit(1);}
	#endif
	
	//printf("\n total runtime to current step is %1.3e sec", runtime(1));
	#if defined(_OPENMP)
	if(runtime(1)>86100){ //86100 sec termination for 24-hour runs!!!
	zeit();
	printf("\n total runtime is %d sec, dumping the data for a restarting point \n", runtime(1));
	dump=fopen(tag,"wb"); dumpall(dump); //fclose(dump); 
	printf("\n dumping complete in %1.3e sec", zeit());
	exit(2);}
	#endif
	}
	#endif
#endif

#if !defined(_OPENMP)
dump=fopen(tag,"wb"); dumpall(dump);
#endif
printf("\n\n Simulation complete \n\n");
}
