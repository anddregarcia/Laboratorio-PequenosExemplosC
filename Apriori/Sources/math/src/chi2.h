
#ifndef __CHI2__
#define __CHI2__

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern double chi2pdf  (double x,    double df);
extern double chi2cdf  (double x,    double df);
extern double chi2cdfP (double x,    double df);
extern double chi2cdfQ (double x,    double df);
extern double chi2qtl  (double prob, double df);
extern double chi2qtlP (double prob, double df);
extern double chi2qtlQ (double prob, double df);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define chi2cdf(x,d)   chi2cdfP(x,d)
#define chi2qtl(p,d)   chi2qtlP(p,d)

#endif
