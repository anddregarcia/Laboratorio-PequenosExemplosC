
#ifndef __REPORT__
#define __REPORT__
#include "tract.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- report modes --- */
#define ISR_ALL        0        /* report all item sets */
#define ISR_CLOSED     1        /* report only closed item sets */
#define ISR_SCAN    0x10        /* report in scanable form */
#define ISR_LOGS    0x20        /* compute sums of logarithms */
#define ISR_DOUBLE  0x40        /* double precision support values */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
struct _isreport;               /* --- an item set eval. function --- */
typedef double ISEVALFN (struct _isreport *rep, void *data);

typedef struct _isreport {      /* --- an item set reporter --- */
  ITEMBASE   *base;             /* underlying item base */
  FILE       *file;             /* output file to write to */
  int        closed;            /* whether to report only closed sets */
  int        min;               /* minimum number of items in set */
  int        max;               /* maximum number of items in set */
  int        cnt;               /* current number of items in set */
  int        pfx;               /* number of items in valid prefix */
  int        *items;            /* current item set (array of items) */
  int        *supps;            /* (prefix) item sets support values */
  double     *sdbls;            /* ditto, as double precision values */
  int        *pexs;             /* perfect extension items */
  int        *pxpp;             /* number of perfect exts. per prefix */
  FILE       *ftid;             /* output file for transaction ids */
  int        *tids;             /* array  of transaction ids */
  int        tidcnt;            /* number of transaction ids */
  double     logwgt;            /* logarithm of total trans. weight */
  double     *logs;             /* logarithms of item frequencies */
  double     *sums;             /* sums of logarithms for prefixes */
  ISEVALFN   *eval;             /* additional evaluation function */
  void       *data;             /* additional evaluation data */
  double     minval;            /* minimal value of evaluation */
  const char *isep;             /* item separator for output */
  const char *impl;             /* implication sign for rule output */
  const char *format;           /* format for information output */
  const char **names;           /* (formatted) item names */
  char       *apos[1];          /* append positions in output buffer */
} ISREPORT;                     /* (item set reporter) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern ISREPORT*   isr_create  (ITEMBASE *base, FILE *file, int mode,
                                const char *sep, const char *impl);
extern void        isr_delete  (ISREPORT *rep, int delis);
extern ITEMBASE*   isr_base    (ISREPORT *rep);
extern FILE*       isr_file    (ISREPORT *rep);

extern void        isr_setfmt  (ISREPORT *rep, const char *format);
extern void        isr_setsize (ISREPORT *rep, int min, int max);
extern void        isr_seteval (ISREPORT *rep, ISEVALFN eval,
                                void *data, double minval);
extern void        isr_setftid (ISREPORT *rep, FILE *ftid);
extern FILE*       isr_getftid (ISREPORT *rep);

extern int         isr_add     (ISREPORT *rep, int item, int supp);
extern int         isr_addx    (ISREPORT *rep, int item, double supp);
extern int         isr_addpex  (ISREPORT *rep, int item);
extern int         isr_uses    (ISREPORT *rep, int item);
extern int         isr_remove  (ISREPORT *rep, int cnt);
extern int         isr_xable   (ISREPORT *rep);

extern int         isr_cnt     (ISREPORT *rep);
extern int         isr_item    (ISREPORT *rep);
extern int         isr_itemx   (ISREPORT *rep, int index);
extern int         isr_supp    (ISREPORT *rep);
extern int         isr_suppx   (ISREPORT *rep, int index);
extern double      isr_logsum  (ISREPORT *rep);
extern double      isr_logsumx (ISREPORT *rep, int index);

extern int         isr_pexcnt  (ISREPORT *rep);
extern int         isr_pex     (ISREPORT *rep, int index);

extern const char* isr_name    (ISREPORT *rep, int item);
extern double      isr_log     (ISREPORT *rep, int item);
extern double      isr_logwgt  (ISREPORT *rep);
extern double      isr_logq    (ISREPORT *rep, void *data);

extern int         isr_report  (ISREPORT *rep);
extern int         isr_reportx (ISREPORT *rep, int *tids, int n);
extern int         isr_sinfo   (ISREPORT *rep, int supp, double eval);
extern int         isr_sinfox  (ISREPORT *rep, double supp,double eval);
extern int         isr_rinfo   (ISREPORT *rep, int supp,
                                int body, int head, double eval);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define isr_base(r)       ((r)->base)
#define isr_file(r)       ((r)->file)
#define isr_setfmt(r,f)   ((r)->format = (f))
#define isr_setftid(r,f)  ((r)->ftid = (f))
#define isr_getftid(r)    ((r)->ftid)

#define isr_uses(r,i)     ((r)->pxpp[i] < 0)
#define isr_xable(r)      ((r)->cnt < (r)->max)

#define isr_cnt(r)        ((r)->cnt)
#define isr_item(r)       ((r)->items[(r)->cnt -1])
#define isr_itemx(r,i)    ((r)->items[i])
#define isr_supp(r)       ((r)->supps[(r)->cnt])
#define isr_suppx(r,i)    ((r)->supps[i])
#define isr_logsum(r)     ((r)->sums [(r)->cnt])
#define isr_logsumx(r,i)  ((r)->sums [i])

#define isr_pexcnt(r)     ((int)((r)->pxpp -(r)->pexs))
#define isr_pex(r,t)      ((r)->pexs [i])

#define isr_name(r,i)     ((r)->names[i])
#define isr_log(r,i)      ((r)->logs [i])
#define isr_logwgt(r)     ((r)->logwgt)

#endif
