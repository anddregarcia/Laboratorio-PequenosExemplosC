
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include "report.h"
#include "scan.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define BS_INT       32         /* buffer size for integer output */
#define BS_DBL      400         /* buffer size for double  output */
#define LN_2        0.69314718055994530942  /* ln(2) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/

ISREPORT* isr_create (ITEMBASE *base, FILE *file, int mode,
                      const char *sep, const char *impl)
{                               /* --- create an item set reporter */
  int        i, k, n;           /* loop variables, buffers */
  ISREPORT   *rep;              /* created item set reporter */
  int        len, sum;          /* length of an item name and sum */
  char       buf[4*TS_SIZE+4];  /* buffer for formatting */
  const char *name;             /* to traverse the item names */

  assert(base && file);         /* check the function arguments */
  n   = ib_cnt(base);           /* get the number of items */
  rep = (ISREPORT*)malloc(sizeof(ISREPORT) +(n+n+1) *sizeof(char*));
  if (!rep) return NULL;        /* allocate the base structure */
  rep->base   = base;           /* store the item base */
  rep->file   = file;           /* and the output file */
  rep->closed = mode & ISR_CLOSED;
  rep->min    = 1;              /* init. the range of item set sizes */
  rep->max    = ib_cnt(base);   /* (minimum and maximum size) */
  rep->cnt    = rep->pfx = 0;   /* init. the number of items */
  rep->eval   = (ISEVALFN*)0;   /* clear add. evaluation function */
  rep->data   = NULL;           /* and the corresponding data */
  rep->names  = (const char**)(rep->apos +n+1);
  *rep->apos  = NULL;           /* organize the pointer arrays */
  rep->logs   = rep->sums = NULL;
  k = (!(mode & ISR_DOUBLE)) ? n+1 : 0;
  rep->items  = (int*)malloc((n+n+n+1+k) *sizeof(int));
  if (!rep->items) { isr_delete(rep, 0); return NULL; }
  rep->pexs   = rep->items + n; /* allocate memory for the arrays */
  rep->pxpp   = rep->pexs += n; /* and organize it (split it) */
  for (i = n; --i >= 0; )       /* clear the item usage flags */
    rep->pxpp[i] = 0;           /* (sign bits) for all items */
  if (k) {                      /* if to use integer support values */
    rep->sdbls = NULL;          /* clear the double prec. array */
    rep->supps = rep->pxpp +k;  /* set pointer to integer array */
    rep->supps[0] = base->wgt;} /* and store the empty set support */
  else {                        /* if to use double precision support */
    rep->supps = NULL;          /* clear the integer array */
    rep->sdbls = (double*)malloc((n+1) *sizeof(double));
    if (!rep->sdbls) { isr_delete(rep, 0); return NULL; }
    rep->sdbls[0] = base->wgt;  /* create a double prec. array */
  }                             /* and store the empty set support */
  rep->ftid   = NULL;           /* clear the transaction id file, */
  rep->tids   = NULL;           /* the array  of transaction ids and */
  rep->tidcnt = 0;              /* the number of transaction ids */
  if (mode & ISR_LOGS) {        /* if to compute logarithms of freqs. */
    rep->logs = (double*)malloc((n+n+1) *sizeof(double));
    if (!rep->logs) { isr_delete(rep, 0); return NULL; }
    rep->sums = rep->logs +n;   /* allocate the needed arrays */
    for (i = 0; i < n; i++)     /* compute logarithms of item freqs. */
      rep->logs[i] = log(ib_getfrq(base, i));
    rep->logwgt  = log(base->wgt);
    rep->sums[0] = 0;           /* store the log of the total weight */
  }                             /* and init. the sum of logarithms */
  for (i = sum = 0; i < n; i++) {
    name = ib_name(base, i);    /* traverse the items and their names */
    len  = strlen(name);        /* get length and format item name */
    sum += k = (mode & ISR_SCAN) ? sc_format(buf, name, 0) : len;
    if (k > len) {              /* if name formatting was needed */
      name = (const char*)malloc((k+1) *sizeof(char));
      if (name) name = strcpy((char*)name, buf);
    }                           /* clone the formatted item name */
    rep->names[i] = name;       /* store the (formatted) item name */
    if (!name) { isr_delete(rep, 0); return NULL; }
  }                             /* check for proper name copying */
  rep->names[n] = NULL;         /* store a sentinel after the names */
  if (!sep) sep = " ";          /* allocate the output buffer */
  *rep->apos = (char*)malloc((sum +(n-1)*strlen(sep) +1) *sizeof(char));
  if (!*rep->apos) { isr_delete(rep, 0); return NULL; }
  rep->isep   = sep;            /* item separator, implication sign */
  rep->impl   = (impl) ? impl : " <- ";
  rep->format = "  (%s)";       /* format for absolute support */
  return rep;                   /* return created item set reporter */
}  /* isr_create() */

/*--------------------------------------------------------------------*/

void isr_delete (ISREPORT *rep, int delis)
{                               /* --- delete an item set reporter */
  int i;                        /* loop variable */

  assert(rep);                  /* check the function argument */
  if (*rep->apos) free(*rep->apos);  /* delete the output buffer */
  for (i = 0; rep->names[i]; i++)    /* traverse the item names */
    if (rep->names[i] != ib_name(rep->base, i))
      free((void*)rep->names[i]);    /* delete all cloned names */
  if (rep->logs)  free(rep->logs);   /* delete the arrays */
  if (rep->items) free(rep->items);  /* (if they are present) */
  if (delis) ib_delete(rep->base);   /* delete the item base */
  free(rep);                    /* delete the base structure */
}  /* isr_delete() */

/*--------------------------------------------------------------------*/

void isr_setsize (ISREPORT *rep, int min, int max)
{                               /* --- set size range for item set */
  assert(rep                    /* check the function arguments */
     && (min >= 0) && (max >= 0));
  rep->min = min;               /* store the minimum and maximum */
  rep->max = max;               /* size of an item set to report */
}  /* isr_setsize() */

/*--------------------------------------------------------------------*/

void isr_seteval (ISREPORT *rep, ISEVALFN eval,
                  void *data, double minval)
{                               /* --- set evaluation function */
  assert(rep);                  /* check the function argument */
  rep->eval   = eval;           /* store the evaluation function, */
  rep->data   = data;           /* the corresponding user data */
  rep->minval = minval;         /* and the minimum value */
}  /* isr_seteval() */

/*--------------------------------------------------------------------*/

int isr_add (ISREPORT *rep, int item, int supp)
{                               /* --- add an item (integer support) */
  assert(rep && (item >= 0)     /* check the function arguments */
             && (item < (int)(rep->names -(const char**)rep->apos)));
  if (isr_uses(rep, item))      /* if the item is already in use, */
    return -1;                  /* abort the function */
  rep->pxpp [item] |= INT_MIN;  /* mark the item as used */
  rep->items[  rep->cnt] = item;/* store the item and its support */
  rep->supps[++rep->cnt] = supp;/* clear the perfect ext. counter */
  rep->pxpp [  rep->cnt] &= INT_MIN;
  return rep->cnt;              /* return the new number of items */
}  /* isr_add() */

/*--------------------------------------------------------------------*/

int isr_addx (ISREPORT *rep, int item, double supp)
{                               /* --- add an item (float support) */
  assert(rep && (item >= 0)     /* check the function arguments */
             && (item < (int)(rep->names -(const char**)rep->apos)));
  if (isr_uses(rep, item))      /* if the item is already in use, */
    return -1;                  /* abort the function */
  rep->pxpp [item] |= INT_MIN;  /* mark the item as used */
  rep->items[  rep->cnt] = item;/* store the item and its support */
  rep->sdbls[++rep->cnt] = supp;/* clear the perfect ext. counter */
  rep->pxpp [  rep->cnt] &= INT_MIN;
  return rep->cnt;              /* return the new number of items */
}  /* isr_addx() */

/*--------------------------------------------------------------------*/

int isr_addpex (ISREPORT *rep, int item)
{                               /* --- add a perfect extension */
  assert(rep && (item >= 0)     /* check the function arguments */
             && (item < (int)(rep->names -(const char**)rep->apos)));
  if (isr_uses(rep, item))      /* if the item is already in use, */
    return -1;                  /* abort the function */
  rep->pxpp[item] |= INT_MIN;   /* mark the item as used */
  *--rep->pexs = item;          /* store the item and */
  rep->pxpp[rep->cnt]++;        /* count it for the current prefix */
  return (int)(rep->pxpp -rep->pexs);
}  /* isr_addpex() */           /* return the number of perf. exts. */

/*--------------------------------------------------------------------*/

int isr_remove (ISREPORT *rep, int cnt)
{                               /* --- remove one or more items */
  int k;                        /* loop variable */

  assert(rep && (rep->cnt >= 0));  /* check the function argument */
  while (--cnt >= 0) {          /* traverse the added items */
    if (rep->cnt <= 0) break;   /* check for items to remove */
    for (k = rep->pxpp[rep->cnt] & ~INT_MIN; --k >= 0; )
      rep->pxpp[*rep->pexs++] &= ~INT_MIN;
    rep->pxpp[rep->items[--rep->cnt]] &= ~INT_MIN;
  }                             /* remove the item markers */
  if (rep->cnt < rep->pfx)      /* remove the item and its perf. exts. */
    rep->pfx = rep->cnt;        /* and adapt the prefix size if nec. */
  return rep->cnt;              /* return the new number of items */
}  /* isr_remove() */

/*--------------------------------------------------------------------*/

int isr_intout (ISREPORT *rep, int num)
{                               /* --- print an integer number */
  int  i, s;                    /* loop variable, sign flag */
  char buf[BS_INT];             /* output buffer */

  assert(rep && rep->file);     /* check the function arguments */
  if (num == 0) {               /* treat zero as a special case */
    fputc('0', rep->file); return 1; }
  if (num <= INT_MIN)           /* treat INT_MIN as a special case */
    return fwrite("-2147483648", sizeof(char), 11, rep->file);
  s = 0;                        /* default: no sign printed */
  if (num < 0) {                /* if the number is negative, */
    fputc('-', rep->file); s = 1; }   /* print a leading sign */
  num = abs(num);               /* remove the sign (just printed) */
  i   = BS_INT;                 /* start at the end of the buffer */
  do {                          /* digit output loop */
    buf[--i] = (num % 10) +'0'; /* store the next digit and */
    num /= 10;                  /* remove it from the number */
  } while (num > 0);            /* while there are more digits */
  fwrite(buf+i, sizeof(char), BS_INT-i, rep->file);
  return BS_INT -i +s;          /* print the digits and */
}  /* isr_intout() */           /* return the number of characters */

/*--------------------------------------------------------------------*/

int isr_dblout (ISREPORT *rep, double num, int dec)
{                               /* --- print a double prec. number */
  int    i, s;                  /* loop variable, sign flag */
  double x, r;                  /* integer value and fraction */
  char   buf[BS_DBL];           /* output buffer */

  assert(rep && rep->file);     /* check the function arguments */
  s = 0;                        /* default: no sign printed */
  if (num < 0) {                /* if the number is negative, */
    fputc('-', rep->file); s = 1; }   /* print a leading sign */
  num = fabs(num);              /* remove the sign (just printed) */
  if (dec > 32) dec = 32;       /* limit the number of decimals */
  r = (dec & 1) ? 0.05 : 0.5; x = 0.1;
  for (i = dec >> 1; i > 0; i >>= 1) {
    x *= x; if (i & 1) r *= x;} /* compute value for rounding and */
  num += r;                     /* round number for given decimals */
  num -= x = floor(num);        /* get integer part and fraction */
  i = BS_DBL;                   /* start at the end of the buffer */
  do {                          /* digit output loop */
    if (--i < 0) break;         /* prevent a buffer overflow */
    buf[i] = (char)fmod(x, 10) +'0';
    x = floor(x/10);            /* compute and store next digit */
  } while (x > 0);              /* while there are more digits */
  fwrite(buf+i, sizeof(char), BS_DBL-i, rep->file);
  i = BS_DBL -i +s;             /* print the stored characters */
  if (dec <= 0)                 /* if to print no decimals, */
    return i;                   /* return the number of characters */
  fputc('.', rep->file); i++;   /* print and count the decimal point */
  while (--dec >= 0) {          /* while to print more decimals */
    num *= 10;                  /* compute the next decimal */
    fputc((int)num +'0', rep->file); i++;
    num -= floor(num);          /* print and count the next decimal */
  }                             /* and remove it from the number */
  return i;                     /* return the number of characters */
}  /* isr_dblout() */

/*--------------------------------------------------------------------*/

double isr_logq (ISREPORT *rep, void *data)
{                               /* --- logarithm of support quotient */
  double x;                     /* buffer for computation */
  assert(rep);                  /* check the function arguments */
  x = (rep->supps) ? log(rep->supps[rep->cnt])
                   : log(rep->sdbls[rep->cnt]);
  return (x -rep->sums[rep->cnt] +(rep->cnt-1) *rep->logwgt) *(1/LN_2);
}  /* isr_logq() */

/* Evaluate an itemset by the logarithm of the quotient of the actual */
/* support of an item set and the support that is expected under full */
/* independence of the items (product of item probabilities times the */
/* total transaction weight). 'data' is needed for the interface.     */

/*--------------------------------------------------------------------*/

static int _output (ISREPORT *rep)
{                               /* --- output an item set */
  int        i;                 /* loop variable */
  char       *s;                /* to traverse the output buffer */
  const char *name;             /* to traverse the item names */
  double     eval = 0;          /* additional evaluation */

  assert(rep                    /* check the function argument */
     && (rep->cnt >= rep->min)
     && (rep->cnt <= rep->max));
  if (rep->eval) {              /* if an evaluation function is given */
    if (rep->logs) {            /* if to compute sums of logarithms */
      eval = rep->sums[rep->pfx]; /* get the valid sum for a prefix */
      for (i = rep->pfx; i < rep->cnt; ) {
        eval += rep->logs[rep->items[i]];
        rep->sums[++i] = eval;  /* traverse the additional items */
      }                         /* and add the logarithms of */
    }                           /* their individual frequencies */
    eval = rep->eval(rep, rep->data);
    if (eval < rep->minval) return 0;
  }                             /* evaluate the current item set */
  s = rep->apos[rep->pfx];      /* get the position for appending */
  if (rep->pfx < rep->cnt) {    /* if items have been added */
    if (rep->pfx <= 0) {        /* if this is the first item */
      for (name = rep->names[rep->items[rep->pfx]]; *name; )
        *s++ = *name++;         /* copy the item name to the buffer */
      rep->apos[++rep->pfx] = s;/* record the position for */
    }                           /* appending the second item */
    while (rep->pfx < rep->cnt){/* traverse the additional items */
      for (name = rep->isep; *name; )
        *s++ = *name++;         /* copy the item separator */
      for (name = rep->names[rep->items[rep->pfx]]; *name; )
        *s++ = *name++;         /* copy the item name to the buffer */
      rep->apos[++rep->pfx] = s;/* record the new position */
    }                           /* for appending the next item */
  }
  fwrite(*rep->apos, sizeof(char), s -*rep->apos, rep->file);
  /* Writing the formatted item set with fwrite seems to be slightly */
  /* faster than terminating the string and writing it with fputs.   */
  if (rep->format) {            /* check for information output */
    if (rep->supps) isr_sinfo (rep, rep->supps[rep->cnt], eval);
    else            isr_sinfox(rep, rep->sdbls[rep->cnt], eval);
  }                             /* call appropriate info function */
  fputc('\n', rep->file);       /* terminate the item set output */
  if (!rep->tids || !rep->ftid) /* check whether to report */
    return 1;                   /* a list of transaction ids */
  if (rep->tidcnt > 0) {        /* if tids are in ascending order */
    for (i = 0; i < rep->tidcnt; i++) {
      if (i > 0) fputs(rep->isep, rep->ftid);
      fprintf(rep->ftid, "%d", rep->tids[i]+1);
    } }                         /* report the transaction ids */
  else if (rep->tidcnt < 0) {   /* if tids are in descending order */
    for (i = -rep->tidcnt; --i >= 0; ) {
      fprintf(rep->ftid, "%d", rep->tids[i]+1);
      if (i > 0) fputs(rep->isep, rep->ftid);
    }                           /* report the transaction ids */
  }                             /* (use the item separator) */
  fputc('\n', rep->ftid);       /* terminate the transaction id list */
  return 1;                     /* return 'item set reported' */
}  /* _output() */

/*--------------------------------------------------------------------*/

static int _report (ISREPORT *rep, int k)
{                               /* --- recursively report item sets */
  int r = 0;                    /* counter for reported item sets */

  assert(rep && (k > 0));       /* check the function arguments */
  while (--k >= 0) {            /* traverse the perfect extensions */
    rep->items[rep->cnt++] = rep->pexs[k];
    if (rep->supps) rep->supps[rep->cnt] = rep->supps[rep->cnt-1];
    else            rep->sdbls[rep->cnt] = rep->sdbls[rep->cnt-1];
    if (rep->cnt >= rep->min)   /* if it has the req. min. size, */
      r += _output(rep);        /* report and count the item set */
    if ((k > 0)                 /* if another item can be added */
    &&  (rep->cnt +k >= rep->min)
    &&  (rep->cnt    <  rep->max))
      r += _report(rep, k);     /* recurse for remaining item sets */
    if (--rep->cnt < rep->pfx)  /* remove the current item again */
      rep->pfx = rep->cnt;      /* and adapt the valid prefix */
  }                             /* if necessary */
  return r;                     /* return number of rep. item sets */
}  /* _report() */

/*--------------------------------------------------------------------*/

int isr_report (ISREPORT *rep)
{                               /* --- report the current item set */
  int r = 0;                    /* number of reported item sets */
  int k;                        /* number of perfect extensions */
  int min, max;                 /* min. and max. item set size */

  assert(rep);                  /* check the function argument */
  if (rep->cnt > rep->max)      /* if the item set is too large, */
    return 0;                   /* abort the function */
  k   = (int)(rep->pxpp -rep->pexs);
  max = rep->cnt +k;            /* get number of perfect extensions */
  min = rep->min;               /* and check whether the minimum size */
  if (max < min) return 0;      /* can be reached with perfect exts. */
  if (rep->closed) {            /* if to report only closed item sets */
    if (max > rep->max) max = rep->max;
    if (max > rep->min) rep->min = max;
  }                             /* adapt the minimal item set size */
  if (rep->cnt >= rep->min)     /* if the item set is large enough, */
    r += _output(rep);          /* report and count it */
  if ((k > 0)                   /* if there are perfect extensions */
  &&  (rep->cnt < rep->max))    /* and maximum size not yet reached, */
    r += _report(rep, k);       /* recursively add and report them */
  rep->min  = min;              /* restore the minimum item set size */
  return r;                     /* return number of rep. item sets */
}  /* isr_report() */

/*--------------------------------------------------------------------*/

int isr_reportx (ISREPORT *rep, int *tids, int n)
{                               /* --- report the current item set */
  assert(rep);                  /* check the function arguments */
  rep->tids   = tids;           /* store the transaction id array */
  rep->tidcnt = n;              /* and the number of transaction ids */
  n = isr_report(rep);          /* report the current item set */
  rep->tids   = NULL;           /* clear the transaction id array */
  return n;                     /* return number of rep. item sets */
}  /* isr_reportx() */

/*--------------------------------------------------------------------*/

static int _getdec (const char *s, const char **end)
{                               /* --- get number of decimal places */
  int k = 0;                    /* number of decimal places */

  assert(s && end);             /* check the function arguments */
  if ((*s >= '0') && (*s <= '9')) {
    k = *s++ -'0';              /* get the first digit */
    if ((*s >= '0') && (*s <= '9'))
      k = 10 *k +*s++ -'0';     /* get a possible second digit and */
  }                             /* compute the number of decimals */
  *end = s; return k;           /* return  the number of decimals */
}  /* _getdec() */

/*--------------------------------------------------------------------*/

int isr_sinfo (ISREPORT *rep, int supp, double eval)
{                               /* --- print item set information */
  int        k, n = 0;          /* number of decimals, char. counter */
  double     wgt;               /* total transaction weight */
  const char *s, *t;            /* to traverse the format */

  assert(rep);                  /* check the function arguments */
  if (!rep->format) return 0;   /* check for a given format */
  wgt = rep->supps[0];          /* get the total transaction weight */
  for (s = rep->format; *s; ) { /* traverse the output format */
    if (*s != '%') {            /* copy everything except '%' */
      fputc(*s++, rep->file); n++; continue; }
    t = s++; k = _getdec(s,&s); /* get the number of decimal places */
    switch (*s++) {             /* evaluate the indicator character */
      case '%': fputc('%', rep->file); n++;                    break;
      case 'a': n += isr_intout(rep,      supp);               break;
      case 's': n += isr_dblout(rep,      supp/wgt,  k);       break;
      case 'S': n += isr_dblout(rep, 100*(supp/wgt), k);       break;
      case 'e': n += isr_dblout(rep,      eval,      k);       break;
      case 'E': n += isr_dblout(rep, 100* eval,      k);       break;
      case  0 : --s;            /* print the requested quantity */
      default : while (t < s) { fputc(*t++, rep->file); n++; } break;
    }                           /* otherwise copy characters */
  }
  return n;                     /* return the number of characters */
}  /* isr_sinfo() */

/*--------------------------------------------------------------------*/

int isr_sinfox (ISREPORT *rep, double supp, double eval)
{                               /* --- print item set information */
  int        k, n = 0;          /* number of decimals, char. counter */
  double     wgt;               /* total transaction weight */
  const char *s, *t;            /* to traverse the format */

  assert(rep);                  /* check the function arguments */
  if (!rep->format) return 0;   /* check for a given format */
  wgt = rep->sdbls[0];          /* get the total transaction weight */
  for (s = rep->format; *s; ) { /* traverse the output format */
    if (*s != '%') {            /* copy everything except '%' */
      fputc(*s++, rep->file); n++; continue; }
    t = s++; k = _getdec(s,&s); /* get the number of decimal places */
    switch (*s++) {             /* evaluate the indicator character */
      case '%': fputc('%', rep->file); n++;                    break;
      case 'a': n += isr_dblout(rep,      supp,      k);       break;
      case 's': n += isr_dblout(rep,      supp/wgt,  k);       break;
      case 'S': n += isr_dblout(rep, 100*(supp/wgt), k);       break;
      case 'e': n += isr_dblout(rep,      eval,      k);       break;
      case 'E': n += isr_dblout(rep, 100* eval,      k);       break;
      case  0 : --s;            /* print the requested quantity */
      default : while (t < s) { fputc(*t++, rep->file); n++; } break;
    }                           /* otherwise copy characters */
  }
  return n;                     /* return the number of characters */
}  /* isr_sinfox() */

/*--------------------------------------------------------------------*/

int isr_rinfo (ISREPORT *rep, int supp, int body, int head, double eval)
{                               /* --- print ass. rule information */
  int        k, n = 0;          /* number of decimals, char. counter */
  double     wgt;               /* total transaction weight */
  double     conf, lift;        /* buffers for computations */
  const char *s, *t;            /* to traverse the format */

  assert(rep);                  /* check the function arguments */
  if (!rep->format) return 0;   /* check for a given format */
  wgt = rep->supps[0];          /* get the total transaction weight */
  if (wgt <= 0) wgt = 1;        /* avoid divisions by zero */
  for (s = rep->format; *s; ) { /* traverse the output format */
    if (*s != '%') {            /* copy everything except '%' */
      fputc(*s++, rep->file); n++; continue; }
    t = s++; k = _getdec(s,&s); /* get the number of decimal places */
    switch (*s++) {             /* evaluate the indicator character */
      case '%': fputc('%', rep->file); n++;                    break;
      case 'a': n += isr_intout(rep,      supp);               break;
      case 's': n += isr_dblout(rep,      supp/wgt,  k);       break;
      case 'S': n += isr_dblout(rep, 100*(supp/wgt), k);       break;
      case 'b': n += isr_intout(rep,      body);               break;
      case 'x': n += isr_dblout(rep,      body/wgt,  k);       break;
      case 'X': n += isr_dblout(rep, 100*(body/wgt), k);       break;
      case 'h': n += isr_intout(rep,      head);               break;
      case 'y': n += isr_dblout(rep,      head/wgt,  k);       break;
      case 'Y': n += isr_dblout(rep, 100*(head/wgt), k);       break;
      case 'c': conf = (body > 0) ? supp/(double)body : 0;
                n += isr_dblout(rep,      conf,      k);       break;
      case 'C': conf = (body > 0) ? supp/(double)body : 0;
                n += isr_dblout(rep, 100* conf,      k);       break;
      case 'l': lift = ((body > 0) && (head > 0))
                     ? (supp*wgt) /(body*(double)head) : 0;
                n += isr_dblout(rep,      lift,      k);       break;
      case 'L': lift = ((body > 0) && (head > 0))
                     ? (supp*wgt) /(body*head) : 0;
                n += isr_dblout(rep, 100* lift,      k);       break;
      case 'e': n += isr_dblout(rep,      eval,      k);       break;
      case 'E': n += isr_dblout(rep, 100* eval,      k);       break;
      case  0 : --s;            /* print the requested quantity */
      default : while (t < s) { fputc(*t++, rep->file); n++; } break;
    }                           /* otherwise copy characters */
  }
  return n;                     /* return the number of characters */
}  /* isr_rinfo() */

/*--------------------------------------------------------------------*/

int isr_rinfox (ISREPORT *rep, double supp,
                double body, double head, double eval)
{                               /* --- print ass. rule information */
  int        k, n = 0;          /* number of decimals, char. counter */
  double     wgt;               /* total transaction weight */
  double     conf, lift;        /* buffers for computations */
  const char *s, *t;            /* to traverse the format */

  assert(rep);                  /* check the function arguments */
  if (!rep->format) return 0;   /* check for a given format */
  wgt = rep->supps[0];          /* get the total transaction weight */
  if (wgt <= 0) wgt = 1;        /* avoid divisions by zero */
  for (s = rep->format; *s; ) { /* traverse the output format */
    if (*s != '%') {            /* copy everything except '%' */
      fputc(*s++, rep->file); n++; continue; }
    t = s++; k = _getdec(s,&s); /* get the number of decimal places */
    switch (*s++) {             /* evaluate the indicator character */
      case '%': fputc('%', rep->file); n++;                    break;
      case 'a': n += isr_dblout(rep,      supp,      k);       break;
      case 's': n += isr_dblout(rep,      supp/wgt,  k);       break;
      case 'S': n += isr_dblout(rep, 100*(supp/wgt), k);       break;
      case 'b': n += isr_dblout(rep,      body,      k);       break;
      case 'x': n += isr_dblout(rep,      body/wgt,  k);       break;
      case 'X': n += isr_dblout(rep, 100*(body/wgt), k);       break;
      case 'h': n += isr_dblout(rep,      head,      k);       break;
      case 'y': n += isr_dblout(rep,      head/wgt,  k);       break;
      case 'Y': n += isr_dblout(rep, 100*(head/wgt), k);       break;
      case 'c': conf = (body > 0) ? supp/body : 0;
                n += isr_dblout(rep,      conf,      k);       break;
      case 'C': conf = (body > 0) ? supp/body : 0;
                n += isr_dblout(rep, 100* conf,      k);       break;
      case 'l': lift = ((body > 0) && (head > 0))
                     ? (supp*wgt) /(body*head) : 0;
      case 'L': lift = ((body > 0) && (head > 0))
                     ? (supp*wgt) /(body*head) : 0;
                n += isr_dblout(rep, 100* lift,      k);       break;
      case 'e': n += isr_dblout(rep,      eval,      k);       break;
      case 'E': n += isr_dblout(rep, 100* eval,      k);       break;
      case  0 : --s;            /* print the requested quantity */
      default : while (t < s) { fputc(*t++, rep->file); n++; } break;
    }                           /* otherwise copy characters */
  }
  return n;                     /* return the number of characters */
}  /* isr_rinfox() */
