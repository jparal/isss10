/* header file for nlibpars_f.c */

void vinput(int iunit, int junit, int kunit, int *icmpl, int *ircopy);

void parsvf(char *page, int *ip, char *group, int *ig, char *code,
            char *helpv, int *icode, double *value, double *range,
            char *cvalue, int lpmax, int lgmax, int nvmax, int ncmax,
            int itwo, int *npage, int *ngroup, int *nvars, int iunit);

void rdnlst(char *code, int *icode, double *value, double *range,
            char *cvalue, int nvars, int ncmax, int itwo,
            int iunit);

void wrnlst(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int iunit);

void wrvarf(char *code, int *icode, double *value, char *cvalue, 
            int nvars, int ncmax, int iunit, int kstrt, int knum,
            int isel, int iquote, int lvm);

void writvf(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int iunit);

void menuvfd(char *page, int *ip, char *group, int *ig, char *code,
             char *helpv, int *icode, double *value, double *range,
             char *cvalue, int npage, int ngroup, int nvars, int ncmax,
             int itwo, int *icmpl, int *ircopy, int iunit, int *irc);

void vpars(const char *input, char *code, int *ic, double *value,
           int *id, char *cval, int *irc);

void rpars(const char *input, int icode, double *range, int itwo);

void findvn(char *code, const char *codev, int nvars, int *nvar);

void findv(const char *input, char *code, int ltou, int *ic, int *irc);

void evalc(const char *cval, int *ival, double *val, int *id);

void frmtv(const char *code, int icode, double *value, char *cvalue, 
           int iquote, char *chrv, int *lv);

void frmtn(const char *code, char *chrv, int *lv);

void dsvarf(char *code, int *icode, double *value, char *cvalue,
            int nvars, int ncmax, int kstrt, int knum, int isel,
            int iquote, int lvm, double ax, double ay, double space);

void opsysc(const char *input);

void helper(double ax, double *ay, double space);

void helpvf(const char *code, double ax, double *ay, double space,
            int iunit);

void wparms(int iunit, int *irc);

void menups(int iunit, double *bvl, double *bvr, int *ndp, int *ndv,
            int *ndd, int *jnmv, int *jpro, int *jps, int *nplot,
            double *amodex, double *freq, double *trmp, double *toff,
            int *ntr, double *vtest, int *ibcs, double *anx,
            double *edge1, double *edge2);

void menup2(int iunit, int *ndp, int *ndv, int *ndd, int *jnmv,
            int *jpro, int *jps, int *nplot, int *ntr, int *ibcs,
            double *anx, double *edge1, double *edge2);
