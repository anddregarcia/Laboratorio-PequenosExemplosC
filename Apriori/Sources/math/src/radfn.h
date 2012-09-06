
#ifndef __RADFN__
#define __RADFN__

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef double RADFN     (double d2, const double *a);
typedef double DERIVFN   (double d2, const double *a, double f);

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern double rf_cauchy  (double d2, const double *a);
extern double rf_gauss   (double d2, const double *a);

extern double rfi_cauchy (double d2, const double *a);
extern double rfi_gauss  (double d2, const double *a);

extern double rfd_cauchy (double d2, const double *a, double f);
extern double rfd_gauss  (double d2, const double *a, double f);

#endif
