
#ifndef __PARAMS__
#define __PARAMS__

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int getints   (char *s, char **end, int n, ...);
extern int getdbls   (char *s, char **end, int n, ...);

extern int getintvec (char *s, char **end, int n, int    *p);
extern int getdblvec (char *s, char **end, int n, double *p);

#endif
