
#ifndef __ARRAYS__
#define __ARRAYS__
#include "fntypes.h"

/*----------------------------------------------------------------------
  Functions for Pointer Arrays
----------------------------------------------------------------------*/
extern void ptr_clear    (void *array, int n);
extern void ptr_copy     (void *src, void *dst, int n);
extern void ptr_move     (void *array, int off, int n, int pos);
extern void ptr_select   (void *array, int n, int k, RANDFN *randfn);
extern void ptr_shuffle  (void *array, int n,        RANDFN *randfn);
extern void ptr_reverse  (void *array, int n);
extern void ptr_qsort    (void *array, int n, CMPFN *cmpfn, void *data);
extern void ptr_heapsort (void *array, int n, CMPFN *cmpfn, void *data);
extern int  ptr_bsearch  (const void *key,
                          void *array, int n, CMPFN *cmpfn, void *data);
extern int  ptr_unique   (void *array, int n, CMPFN *cmpfn, void *data,
                          OBJFN *delfn);

/*----------------------------------------------------------------------
  Functions for Arrays of Basic Data Types
----------------------------------------------------------------------*/
extern void sht_clear    (short  *array, int n);
extern void sht_copy     (short  *src,   short *dst, int n);
extern void sht_move     (short  *array, int off, int n, int pos);
extern void sht_select   (short  *array, int n, int k, RANDFN *randfn);
extern void sht_shuffle  (short  *array, int n,        RANDFN *randfn);
extern void sht_reverse  (short  *array, int n);
extern void sht_qsort    (short  *array, int n);
extern void sht_heapsort (short  *array, int n);
extern int  sht_unique   (short  *array, int n);
extern int  sht_bsearch  (short  key, short *array, int n);

/*--------------------------------------------------------------------*/

extern void int_clear    (int    *array, int n);
extern void int_copy     (int    *src,   int *dst, int n);
extern void int_move     (int    *array, int off, int n, int pos);
extern void int_select   (int    *array, int n, int k, RANDFN *randfn);
extern void int_shuffle  (int    *array, int n,        RANDFN *randfn);
extern void int_reverse  (int    *array, int n);
extern void int_qsort    (int    *array, int n);
extern void int_heapsort (int    *array, int n);
extern int  int_unique   (int    *array, int n);
extern int  int_bsearch  (int    key, int *array, int n);

/*--------------------------------------------------------------------*/

extern void flt_clear    (float  *array, int n);
extern void flt_copy     (float  *src,   float *dst, int n);
extern void flt_move     (float  *array, int off, int n, int pos);
extern void flt_select   (float  *array, int n, int k, RANDFN *randfn);
extern void flt_shuffle  (float  *array, int n,        RANDFN *randfn);
extern void flt_reverse  (float  *array, int n);
extern void flt_qsort    (float  *array, int n);
extern void flt_heapsort (float  *array, int n);
extern int  flt_unique   (float  *array, int n);
extern int  flt_bsearch  (float  key, float *array, int n);

/*--------------------------------------------------------------------*/

extern void dbl_clear    (double *array, int n);
extern void dbl_copy     (double *src,   double *dst, int n);
extern void dbl_move     (double *array, int off, int n, int pos);
extern void dbl_select   (double *array, int n, int k, RANDFN *randfn);
extern void dbl_shuffle  (double *array, int n,        RANDFN *randfn);
extern void dbl_reverse  (double *array, int n);
extern void dbl_qsort    (double *array, int n);
extern void dbl_heapsort (double *array, int n);
extern int  dbl_unique   (double *array, int n);
extern int  dbl_bsearch  (double key, double *array, int n);

#endif
