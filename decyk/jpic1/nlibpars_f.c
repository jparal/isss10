/* parsing library */
/* Wrappers for calling the Fortran routines from a C main program */

#include <string.h>
#include "nlibpars_f.h"

void vinput_(int *iunit, int *junit, int *kunit, int *icmpl, 
             int *ircopy);

void parsvf_(char *page, int *ip, char *group, int *ig, char *code,
             char *helpv, int *icode, double *value, double *range,
             char *cvalue, int *lpmax, int *lgmax, int *nvmax,
             int *ncmax, int *itwo, int *npage, int *ngroup, int *nvars,
             int *iunit, int l1, int l2, int l3, int l4, int l5);

void rdnlst_(char *code, int *icode, double *value, double *range,
             char *cvalue, int *nvars, int *ncmax, int *itwo,
             int *iunit, int l1, int l2);

void wrnlst_(char *code, int *icode, double *value, char *cvalue,
             int *nvars, int *ncmax, int *iunit, int l1, int l2);

void wrvarf_(char *code, int *icode, double *value, char *cvalue, 
             int *nvars, int *ncmax, int *iunit, int *kstrt, int *knum,
             int *isel, int *iquote, int *lvm, int l1, int l2);

void writvf_(char *code, int *icode, double *value, char *cvalue,
             int *nvars, int *ncmax, int *iunit, int l1, int l2);

void menuvfd_(char *page, int *ip, char *group, int *ig, char *code,
              char *helpv, int *icode, double *value, double *range,
              char *cvalue, int *npage, int *ngroup, int *nvars,
              int *ncmax, int *itwo, int *icmpl, int *ircopy, 
              int *iunit, int *irc, int l1, int l2, int l3, int l4,
              int l5);

void vpars_(const char *input, char *code, int *ic, double *value,
            int *id, char *cval, int *irc, int l1, int l2, int l3);

void rpars_(const char *input, int *icode, double *range, int *itwo,
            int l1);

void findvn_(char *code, const char *codev, int *nvars, int *nvar,
             int l1, int l2);

void findv_(const char *input, char *code, int *ltou, int *ic, int *irc,
            int l1, int l2);

void evalc_(const char *cval, int *ival, double *val, int *id, int l1);

void frmtv_(const char *code, int *icode, double *value, char *cvalue, 
            int *iquote, char *chrv, int *lv, int l1, int l2, int l3);

void frmtn_(const char *code, char *chrv, int *lv, int l1, int l2);

void dsvarf_(char *code, int *icode, double *value, char *cvalue,
             int *nvars, int *ncmax, int *kstrt, int *knum, int *isel,
             int *iquote, int *lvm, double *ax, double *ay,
             double *space, int l1, int l2);

void opsysc_(const char *input, int l1);

void helper_(double *ax, double *ay, double *space);

void helpvf_(const char *code, double *ax, double *ay, double *space,
             int *iunit, int l1);

void wparms_(int *iunit, int *irc);

void menups_(int *iunit, double *bvl, double *bvr, int *ndp, int *ndv,
             int *ndd, int *jnmv, int *jpro, int *jps, int *nplot,
             double *amodex, double *freq, double *trmp, double *toff,
             int *ntr, double *vtest, int *ibcs, double *anx,
             double *edge1, double *edge2);

void menup2_(int *iunit, int *ndp, int *ndv, int *ndd, int *jnmv,
             int *jpro, int *jps, int *nplot, int *ntr, int *ibcs,
             double *anx, double *edge1, double *edge2);


/* Interfaces to C */

/*--------------------------------------------------------------------*/
void vinput(int iunit, int junit, int kunit, int *icmpl, int *ircopy) {
   vinput_(&iunit,&junit,&kunit,icmpl,ircopy);
   return;
}

/*--------------------------------------------------------------------*/
void parsvf(char *page, int *ip, char *group, int *ig, char *code,
            char *helpv, int *icode, double *value, double *range,
            char *cvalue, int lpmax, int lgmax, int nvmax, int ncmax,
            int itwo, int *npage, int *ngroup, int *nvars, int iunit) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   parsvf_(page,ip,group,ig,code,helpv,icode,value,range,cvalue,&lpmax,
           &lgmax,&nvmax,&ncmax,&itwo,npage,ngroup,nvars,&iunit,
           strlen(page),strlen(group),strlen(code),strlen(helpv),
           strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void rdnlst(char *code, int *icode, double *value, double *range,
            char *cvalue, int nvars, int ncmax, int itwo,
            int iunit) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   rdnlst_(code,icode,value,range,cvalue,&nvars,&ncmax,&itwo,&iunit,
           strlen(code),strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void wrnlst(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int iunit) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   wrnlst_(code,icode,value,cvalue,&nvars,&ncmax,&iunit,strlen(code),
           strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void wrvarf(char *code, int *icode, double *value, char *cvalue, 
            int nvars, int ncmax, int iunit, int kstrt, int knum,
            int isel, int iquote, int lvm) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   wrvarf_(code,icode,value,cvalue, &nvars,&ncmax,&iunit,&kstrt,&knum,
           &isel,&iquote,&lvm,strlen(code),strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void writvf(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int iunit) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   writvf_(code,icode,value,cvalue,&nvars,&ncmax,&iunit,strlen(code),
           strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void menuvfd(char *page, int *ip, char *group, int *ig, char *code,
             char *helpv, int *icode, double *value, double *range,
             char *cvalue, int npage, int ngroup, int nvars, int ncmax,
             int itwo, int *icmpl, int *ircopy, int iunit, int *irc) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   menuvfd_(page,ip,group,ig,code,helpv,icode,value,range,cvalue,&npage,
            &ngroup,&nvars,&ncmax,&itwo,icmpl,ircopy,&iunit,irc,
            strlen(page),strlen(group),strlen(code),strlen(helpv),
            strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void vpars(const char *input, char *code, int *ic, double *value,
           int *id, char *cval, int *irc) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   vpars_(input,code,ic,value,id,cval,irc,strlen(input),strlen(code),
          strlen(cval));
   return;
}

/*--------------------------------------------------------------------*/
void rpars(const char *input, int icode, double *range, int itwo) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   rpars_(input,&icode,range,&itwo,strlen(input));
   return;
}

/*--------------------------------------------------------------------*/
void findvn(char *code, const char *codev, int nvars, int *nvar) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   findvn_(code,codev,&nvars,nvar,strlen(code),strlen(codev));
   return;
}

/*--------------------------------------------------------------------*/
void findv(const char *input, char *code, int ltou, int *ic, int *irc) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   findv_(input,code,&ltou,ic,irc,strlen(input),strlen(code));
   return;
}

/*--------------------------------------------------------------------*/
void evalc(const char *cval, int *ival, double *val, int *id) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   evalc_(cval,ival,val,id,strlen(cval));
   return;
}

/*--------------------------------------------------------------------*/
void frmtv(const char *code, int icode, double *value, char *cvalue, 
           int iquote, char *chrv, int *lv) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   frmtv_(code,&icode,value,cvalue,&iquote,chrv,lv,strlen(code),
          strlen(cvalue),strlen(chrv));
   return;
}

/*--------------------------------------------------------------------*/
void frmtn(const char *code, char *chrv, int *lv) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   frmtn_(code,chrv,lv,strlen(code),strlen(chrv));
   return;
}

/*--------------------------------------------------------------------*/
void dsvarf(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int kstrt, int knum, int isel,
            int iquote, int lvm, double ax, double ay, double space) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   dsvarf_(code,icode,value,cvalue,&nvars,&ncmax,&kstrt,&knum,&isel,
           &iquote,&lvm,&ax,&ay,&space,strlen(code),strlen(cvalue));
   return;
}

/*--------------------------------------------------------------------*/
void opsysc(const char *input) {
   opsysc_(input,strlen(input));
   return;
}

/*--------------------------------------------------------------------*/
void helper(double ax, double *ay, double space) {
   helper_(&ax,ay,&space);
   return;
}

/*--------------------------------------------------------------------*/
void helpvf(const char *code, double ax, double *ay, double space,
            int iunit) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   helpvf_(code,&ax,ay,&space,&iunit,strlen(code));
   return;
}

/*--------------------------------------------------------------------*/
void wparms(int iunit, int *irc) {
   wparms_(&iunit,irc);
   return;
}

/*--------------------------------------------------------------------*/
void menups(int iunit, double *bvl, double *bvr, int *ndp, int *ndv,
            int *ndd, int *jnmv, int *jpro, int *jps, int *nplot,
            double *amodex, double *freq, double *trmp, double *toff,
            int *ntr, double *vtest, int *ibcs, double *anx,
            double *edge1, double *edge2) {
   menups_(&iunit,bvl,bvr,ndp,ndv,ndd,jnmv,jpro,jps,nplot,amodex,freq,
           trmp,toff,ntr,vtest,ibcs,anx,edge1,edge2);
   return;
}

/*--------------------------------------------------------------------*/
void menup2(int iunit, int *ndp, int *ndv, int *ndd, int *jnmv,
            int *jpro, int *jps, int *nplot, int *ntr, int *ibcs,
            double *anx, double *edge1, double *edge2) {
   menup2_(&iunit,ndp,ndv,ndd,jnmv,jpro,jps,nplot,ntr,ibcs,anx,edge1,
           edge2);
   return;
}