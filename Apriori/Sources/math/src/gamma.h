
#ifndef __GAMMA__
#define __GAMMA__

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern double logGamma   (double n);
extern double Gamma      (double n);
extern double lowerGamma (double n, double x);
extern double upperGamma (double n, double x);
extern double GammaP     (double n, double x);
extern double GammaQ     (double n, double x);

extern double Gammapdf   (double x,    double k, double theta);
extern double Gammacdf   (double x,    double k, double theta);
extern double GammacdfP  (double x,    double k, double theta);
extern double GammacdfQ  (double x,    double k, double theta);
extern double Gammaqtl   (double prob, double k, double theta);
extern double GammaqtlP  (double prob, double k, double theta);
extern double GammaqtlQ  (double prob, double k, double theta);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define Gammacdf(x,k,t)  GammaP(k,(x)/(t))
#define GammacdfP(x,k,t) GammaP(k,(x)/(t))
#define GammacdfQ(x,k,t) GammaQ(k,(x)/(t))

#endif
