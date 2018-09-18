#ifndef __RADAR5_H__
#define __RADAR5_H__

#ifdef __cplusplus
extern "C" {
#endif


/* Header for FORTRAN functions and subroutines of RADAR5 */    

typedef double (*PHI_t)(int *I, double *X, double *RPAR, int *IPAR);
typedef double (*ARGLAG_t)(int *IL, double *X, double *Y,double *RPAR, int *IPAR, PHI_t PHI, double *PAST, int *IPAST, int *NRDS);
typedef void (*FCN_t)(int *N, double *X, double *Y, double *F, ARGLAG_t ARGLAG, PHI_t PHI, double *RPAR, int *IPAR, double *PAST,int *IPAST, int *NRDS);
typedef void (*JFCN_t)(int *N, double *X, double *Y, double *DFY, int LDFY, ARGLAG_t ARGLAG, PHI_t PHI, double *RPAR, int *IPAR, double *PAST,int *IPAST, int *NRDS);
typedef void (*JACLAG_t)(int *N, double *X, double *Y, int *DFYL, ARGLAG_t ARGLAG, PHI_t PHI, int *IVE, int *IVC, int *IVL, double *RPAR, int *IPAR, double *PAST, int *IPAST, int *NRDS);
typedef void (*SOLOUT_t)(int *NR, double *XOLD, double *X, double *HSOL, double *Y, double *CONT, int *LRC, int *N, double *RPAR, int *IPAR, int *IRTRN);


extern int radar5_(int *N, FCN_t FCN, PHI_t PHI, ARGLAG_t ARGLAG, double *X, double *Y, double *XEND, double *H, \
                       double *RTOL, double *ATOL, double *ITOL, \
                       JFCN_t JFCN, int *IJAC, int *MLJAC, int *MUJAC, \
                       JACLAG_t JACLAG, int *NLAGS, int *NJACL, \
                       int *IMAS, SOLOUT_t SOLOUT, int *IOUT, \
                       double *WORK, int *IWORK, double *RPAR, int *IPAR, int *IDID, \
                       double *GRID, int *IPAST, int *DUMMY, int *MLMAS, int *MUMAS);

extern void lagr5_(int *IL, double *X,double *Y,ARGLAG_t ARGLAG,double *PAST,double *ALPHA1,int *IPOS1,double *RPAR,int *IPAR,PHI_t PHI, int *IPAST, int *NRDS);
extern double ylagr5_(int *IP, double *ALPHA1,int *IPOS1,PHI_t PHI,double *RPAR,int *IPAR,double *PAST,int *IPAST, int *NRDS);

extern double contr5_(int *I, int *N, double *XOUT, double *CONT, double *X, double *HSOL);


#ifdef __cplusplus
}
#endif

#endif/*__RADAR5_H__*/
