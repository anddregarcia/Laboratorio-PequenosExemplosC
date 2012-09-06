/*
Projeto de Organização e recuperação da Informação (ORI)
Escola de Engenharia de Piracicaba
Authors: Adilson Perecin(a.perecin(at)hotmail.com) & Cristiano Benato(benato(at)hst.com.br)

Tema: Organização e verificação de itemsets e ordenação de coincidencias através
do algoritmo de ordenação Apriori

Instrutor: Luiz Camolesi
 


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "istree.h"
#include "chi2.h"
#ifdef STORAGE
#include "storage.h"
#endif


#define LN_2        0.69314718055994530942  /* ln(2) */
#define BLKSIZE     32          /* block size for level array */
#define F_HDONLY    INT_MIN     /* flag for head only item in path */
#define F_SKIP      INT_MIN     /* flag for subtree skipping */
#define ID(n)       ((int)((n)->id & ~F_HDONLY))
#define HDONLY(n)   ((int)((n)->id &  F_HDONLY))
#define COUNT(n)    ((n) & ~F_SKIP)
#define CHCNT(n)    ((n)->chcnt & ~F_SKIP)

#ifdef ARCH64                   /* if 64 bit architecture, */
#define PAD(x)      ((x) & 1)   /* padding to an even number */
#else                           /* if sizeof(int) == sizeof(void*), */
#define PAD(x)      0           /* no padding is needed */
#endif


typedef double EVALFN (int supp, int body, int head, int n);

typedef double AGGRFN (double aggr, double val);




static double _conf (int supp, int body, int head, int base)
{                               /* --- rule confidence */
  return (body > 0) ? supp/(double)body : 0;
}  /* _conf() */



static double _diff (int supp, int body, int head, int base)
{                               /* --- absolute confidence difference */
  if ((body <= 0) || (base <= 0)) return 0;
  return fabs(supp/(double)body -head/(double)base);
}  /* _diff() */



static double _lift (int supp, int body, int head, int base)
{                               /* --- lift value */
  if ((body <= 0) || (head <= 0)) return 0;
  return (supp*(double)base) /(body*(double)head);
  /* =   (supp/(double)body) /(head/(double)base) */
}  /* _lift() */



static double _ld21 (int supp, int body, int head, int base)
{                               /* --- abs. diff. of lift value to 1 */
  if ((body <= 0) || (head <= 0)) return 0;
  return fabs((supp*(double)base) /(body*(double)head) -1);
}  /* _ld21() */

/*--------------------------------------------------------------------*/

static double _quot (int supp, int body, int head, int base)
{                               /* --- diff. of lift quotient to 1 */
  double t;                     /* temporary buffer */
  if ((body <= 0) || (head <= 0)) return 0;
  t = (supp*(double)base) /(body*(double)head);
  return 1 -((t > 1) ? 1/t : t);
}  /* _quot() */



static double _chi2 (int supp, int body, int head, int base)
{                               /* --- normalized chi^2 measure */
  double t;                     /* temporary buffer */

  if ((head <= 0) || (head >= base)
  ||  (body <= 0) || (body >= base))
    return 0;                   /* check for strict positivity */
  t = head *(double)body -supp *(double)base;
  return (t*t) /(((double)head)*(base-head)*((double)body)*(base-body));
}  /* _chi2() */                /* compute and return chi^2 measure */


static double _pval (int supp, int body, int head, int base)
{                               /* --- p-value from chi^2 measure */
  return chi2cdf(base *_chi2(supp, body, head, base), 1);
}  /* _pval() */


static double _info (int supp, int body, int head, int base)
{                               /* --- information diff. to prior */
  double sum, t;                /* result, temporary buffer */

  if ((head <= 0) || (head >= base)
  ||  (body <= 0) || (body >= base))
    return 0;                   /* check for strict positivity */
  t = supp; sum = 0;            /* support of     head and     body */
  if (t > 0) sum += t *log(t /(      head  *(double)      body ));
  t = body -supp;               /* support of not head and     body */
  if (t > 0) sum += t *log(t /((base-head) *(double)      body ));
  t = head -supp;               /* support of     head and not body */
  if (t > 0) sum += t *log(t /(      head  *(double)(base-body)));
  t = base -head -body +supp;   /* support of not head and not body */
  if (t > 0) sum += t *log(t /((base-head) *(double)(base-body)));
  return (log(base) +sum/base) /LN_2;
}  /* _info() */                /* return information gain in bits */


static double _pgst (int supp, int body, int head, int base)
{                               /* --- p-value from G statistic */
  return chi2cdf(base *2*LN_2 *_info(supp, body, head, base), 1);
}  /* _pgst() */

static EVALFN *_evalfns[] = {   /* --- add. evaluation measures */
  /* IST_NONE  0 */ (EVALFN*)0, /* no additional evaluation measure */
  /* IST_CONF  1 */ _conf,      /* rule confidence */
  /* IST_DIFF  2 */ _diff,      /* abs. confidence diff. to prior */
  /* IST_LIFT  3 */ _lift,      /* lift value (confidence/prior) */
  /* IST_LD21  4 */ _ld21,      /* abs. diff. of lift value to 1 */
  /* IST_QUOT  5 */ _quot,      /* diff. of lift quotient to 1 */
  /* IST_CHI2  6 */ _chi2,      /* normalized chi^2 measure */
  /* IST_PVAL  7 */ _pval,      /* p-value from chi^2 measure */
  /* IST_INFO  8 */ _info,      /* information difference to prior */
  /* IST_PGST  9 */ _pgst,      /* p-value from G statistic */
  /* IST_LOGQ 10 */ (EVALFN*)0, /* binary log. of support quotient */
};                              /* table of evaluation functions */



static int _search (int id, ISNODE **chn, int n)
{                               /* --- find a child node (index) */
  int i, k, x;                  /* left and middle index */

  assert(chn && (n > 0));       /* check the function arguments */
  for (i = 0; i < n; ) {        /* while the range is not empty */
    k = (i+n) >> 1;             /* get index of the middle element */
    x = ID(chn[k]);             /* compare the item identifier */
    if      (id > x) i = k+1;   /* to the middle element and */
    else if (id < x) n = k;     /* adapt the range boundaries */
    else return k;              /* if there is an exact match, */
  }                             /* return the child node index */
  return -1-i;                  /* return the insertion position */
}  /* _search() */



static int _getsupp (ISNODE *node, int *items, int n)
{                               /* --- get support of an item set */
  int    i, k;                  /* array indices, number of children */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* child node array */

  assert(node                   /* check the function arguments */
     && (n >= 0) && (items || (n <= 0)));
  while (--n > 0) {             /* follow the set/path from the node */
    k = CHCNT(node);            /* if there are no children, */
    if (k <= 0) return F_SKIP;  /* the support is less than minsupp */
    if (node->offset >= 0) {    /* if a pure array is used */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      i = *items++ -ID(chn[0]); /* compute the child array index and */
      if (i >= k) return F_SKIP; }  /* abort if child does not exist */
    else {                      /* if an identifier map is used */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      i = _search(*items++, chn, k);
    }                           /* search for the proper index */
    if (i < 0) return F_SKIP;   /* abort if index is out of range */
    node = chn[i];              /* go to the corresponding child */
    if (!node) return F_SKIP;   /* if the child does not exists, */
  }                             /* the support is less than minsupp */
  if (node->offset >= 0) {      /* if a pure array is used, */
    i = *items -node->offset;   /* compute the counter index */
    if (i >= node->size) return F_SKIP; }
  else {                        /* if an identifier map is used */
    map = node->cnts +(k = node->size);
    i   = int_bsearch(*items, map, k);
  }                             /* search for the proper index */
  if (i < 0) return F_SKIP;     /* abort if index is out of range */
  return node->cnts[i];         /* return the item set support */
}  /* _getsupp() */



static double _min (double aggr, double val)
{ return (val < aggr) ? val : aggr; }

static double _max (double aggr, double val)
{ return (val > aggr) ? val : aggr; }

static double _sum (double aggr, double val)
{ return aggr +val; }


static AGGRFN *_aggrfns[] = {   /* --- measure aggregation functions */
  /* IST_NONE 0 */  (AGGRFN*)0, /* no aggregation (use first value) */
  /* IST_MIN  1 */  _min,       /* minimum */
  /* IST_MAX  2 */  _max,       /* maximum */
  /* IST_AVG  3 */  _sum,       /* average (needs to sum first) */
};

static double _aggregate (ISTREE *ist)
{                               /* --- average rule confidence */
  int    n, item;               /* loop variable, current (head) item */
  int    supp;                  /* support of item set */
  int    body, head;            /* support of rule body and head */
  int    base;                  /* total transaction weight */
  int    *path;                 /* path to follow for body support */
  ISNODE *node, *curr;          /* to traverse the nodes */
  EVALFN *eval;                 /* add. evaluation function */
  AGGRFN *aggr;                 /* aggregation function */
  double val;                   /* (aggregated) value of measure */

  assert(ist);                  /* check the function argument */
  node = ist->node;             /* get the item set support */
  if (node->offset >= 0) item = ist->index +node->offset;
  else                   item = node->cnts[node->size +ist->index];
  supp = COUNT(node->cnts[ist->index]);
  head = COUNT(ist->lvls[0]->cnts[item]);
  base = COUNT(ist->wgt);       /* get head and empty set support */
  eval = _evalfns[ist->eval];   /* get the evaluation function */
  curr = node->parent;          /* get subset support from parent */
  if (!curr)                    /* if there is no parent (root node), */
    return eval(supp, base, head, base); /* evaluate the set directly */
  if (curr->offset >= 0)        /* if a pure array is used */
    body = COUNT(curr->cnts[ID(node) -curr->offset]);
  else {                        /* if an identifier map is used */
    path = curr->cnts +(n = curr->size);
    body = COUNT(curr->cnts[int_bsearch(ID(node), path, n)]);
  }                             /* find index and get subset support */
  val = eval(supp, body, head, base);
  if (ist->agg <= IST_NONE)     /* compute the first measure value */
    return val;                 /* and check whether to return it */
  path = ist->buf +ist->maxht;  /* initialize the path */
  *--path = item;               /* for the support retrieval */
  item = ID(node);              /* get the next head item */
  aggr = _aggrfns[ist->agg];    /* and the aggregation function */
  for (n = 0; curr; curr = curr->parent) {
    body = COUNT(_getsupp(curr, path, ++n));
    val  = aggr(val, eval(supp, body, head, base));
    *--path = item;             /* get the support of the rule body */
    item = ID(curr);            /* and sum the rule confidences, */
  }                             /* then extend the path (store head) */
  return (ist->agg >= IST_AVG) ? val /(n+1) : val;
}  /* _aggregate() */           /* return the measure aggregate */

static double _logq (ISTREE *ist)
{                               /* --- logarithm of support quotient */
  int    n, item;               /* loop variable, current item */
  int    *cnts;                 /* array of item frequencies */
  ISNODE *node;                 /* to traverse the nodes */
  double sum;                   /* sum of logs. of item frequencies */

  assert(ist);                  /* check the function argument */
  cnts = ist->lvls[0]->cnts;    /* get the item frequencies */
  node = ist->node;             /* and the item set support */
  if (node->offset >= 0) item = ist->index +node->offset;
  else                   item = node->cnts[node->size +ist->index];
  sum = log(COUNT(node->cnts[ist->index])) -log(COUNT(cnts[item]));
  for (n = 0; node->parent; node = node->parent) {
    sum -= log(COUNT(cnts[ID(node)])); n++; }
  if (n > 0) sum += n *log(COUNT(ist->wgt));
  return sum /LN_2;             /* sum logs. of item frequencies and */
}  /* _logq() */                /* subtract from log. of set freq., */


static void _count (ISNODE *node,
                    const int *items, int n, int wgt, int min)
{                               /* --- count transaction recursively */
  int    i, k, o;               /* array index, offset, map size */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* array of child nodes */

  assert(node                   /* check the function arguments */
     && (n >= 0) && (items || (n <= 0)));
  if (node->offset >= 0) {      /* if a pure array is used */
    if (node->chcnt == 0) {     /* if this is a new node (leaf) */
      o = node->offset;         /* get the index offset */
      while ((n > 0) && (*items < o)) {
        n--; items++; }         /* skip items before first counter */
      while (--n >= 0) {        /* traverse the transaction's items */
        i = *items++ -o;        /* compute the counter array index */
        if (i >= node->size) return;
        node->cnts[i] += wgt;   /* if the corresp. counter exists, */
      } }                       /* add the transaction weight to it */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      o   = ID(chn[0]);         /* get the child node array */
      while ((n >= min) && (*items < o)) {
        n--; items++; }         /* skip items before the first child */
      for (--min; --n >= min;){ /* traverse the transaction's items */
        i = *items++ -o;        /* compute the child array index */
        if (i >= node->chcnt) return;
        if (chn[i]) _count(chn[i], items, n, wgt, min);
      }                         /* if the corresp. child node exists, */
    } }                         /* count the transaction recursively */
  else {                        /* if an identifer map is used */
    if (node->chcnt == 0) {     /* if this is a new node (leaf) */
      map = node->cnts +(k = node->size);
      o   = map[0];             /* get the identifier map */
      while ((n > 0) && (*items < o)) {
        n--; items++; }         /* skip items before first counter */
      o   = map[k-1];           /* get the last item with a counter */
      for (i = 0; --n >= 0; ) { /* traverse the transaction's items */
        if (*items > o) return; /* if beyond last item, abort */
        #ifdef IST_BSEARCH      /* if to use a binary search */
        i = int_bsearch(*items++, map, k);
        if (i >= 0) node->cnts[i] += wgt;
        #else                   /* if to use a linear search */
        while (*items > map[i]) i++;
        if (*items++ == map[i]) node->cnts[i] += wgt;
        #endif                  /* if the corresp. counter exists, */
      } }                       /* add the transaction weight to it */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      o   = ID(chn[0]);         /* get the child node array */
      while ((n >= min) && (*items < o)) {
        n--; items++; }         /* skip items before first child */
      k   = node->chcnt;        /* get the number of children and */
      o   = ID(chn[k-1]);       /* the index of the last item */
      for (--min; --n >= min; ) {
        if (*items > o) return; /* traverse the transaction */
        #ifdef IST_BSEARCH      /* if to use a binary search */
        i = _search(*items++, chn, k);
        if (i >= 0) _count(chn[i], items, n, wgt, min);
        else        i = -1-i;   /* count the transaction recursively */
        chn += i; k -= i;       /* and adapt the child node range */
        #else                   /* if to use a linear search */
        while (*items > ID(*chn)) chn++;
        if (*items++ == ID(*chn)) _count(*chn, items, n, wgt, min);
        #endif                  /* find the proper child node index */
      }                         /* if the corresp. child node exists, */
    }                           /* count the transaction recursively */
  }
}  /* _count() */


static void _countx (ISNODE *node, const TTNODE *tree, int min)
{                               /* --- count trans. tree recursively */
  int    i, k, o, n;            /* array indices, loop variables */
  int    item;                  /* buffer for an item */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* child node array */

  assert(node && tree);         /* check the function arguments */
  if (ttn_max(tree) < min)      /* if the transactions are too short, */
    return;                     /* abort the recursion */
  n = ttn_size(tree);           /* get the number of children */
  if (n <= 0) {                 /* if there are no children */
    if (n < 0) _count(node, ttn_items(tree), -n, ttn_wgt(tree), min);
    return;                     /* count the normal transaction */
  }                             /* and abort the function */
  while (--n >= 0)              /* count the transactions recursively */
    _countx(node, ttn_child(tree, n), min);
  if (node->offset >= 0) {      /* if a pure array is used */
    if (node->chcnt == 0) {     /* if this is a new node (leaf) */
      o = node->offset;         /* get the index offset */
      for (n = ttn_size(tree); --n >= 0; ) {
        i = ttn_item(tree,n)-o; /* traverse the node's items */
        if (i < 0) return;      /* if before the first item, abort */
        if (i < node->size)     /* if the corresp. counter exists */
          node->cnts[i] += ttn_wgt(ttn_child(tree, n));
      } }                       /* add the transaction weight to it */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      o   = ID(chn[0]);         /* get the child node array */
      for (--min, n = ttn_size(tree); --n >= 0; ) {
        i = ttn_item(tree,n)-o; /* traverse the node's items */
        if (i < 0) return;      /* if before the first item, abort */
        if ((i < node->chcnt) && chn[i])
          _countx(chn[i], ttn_child(tree, n), min);
      }                         /* if the corresp. child node exists, */
    } }                         /* count the trans. tree recursively */
  else {                        /* if an identifer map is used */
    if (node->chcnt == 0) {     /* if this is a new node (leaf) */
      map = node->cnts +(k = node->size);
      o   = map[0];             /* get the item identifier map */
      for (n = ttn_size(tree); n >= 0; --n) {
        item = ttn_item(tree,n);/* traverse the node's items */
        if (item < o) return;   /* if before the first item, abort */
        #ifdef IST_BSEARCH      /* if to use a binary search */
        i = int_bsearch(item, map, k);
        if (i >= 0) node->cnts[k = i] += ttn_wgt(ttn_child(tree, n));
        else        k = -1-i;   /* add trans. weight to the counter */
        #else                   /* if to use a linear search */
        while (item < map[--k]);
        if (item == map[k]) node->cnts[k] += ttn_wgt(ttn_child(tree,n));
        else k++;               /* if the corresp. counter exists, */
        #endif                  /* add the transaction weight to it, */
      } }                       /* otherwise adapt the map index */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      k   = node->chcnt;        /* get the child node array and */
      o   = ID(chn[0]);         /* the last item with a child */
      for (--min, n = ttn_size(tree); --n >= 0; ) {
        item = ttn_item(tree,n);/* traverse the node's items */
        if (item < o) return;   /* if before the first item, abort */
        #ifdef IST_BSEARCH      /* if to use a binary search */
        i = _search(item, chn, k);
        if (i >= 0) _countx(chn[i], ttn_child(tree, n), min);
        else        k = -1-i;   /* add trans. weight to the counter */
        #else                   /* if to use a linear search */
        while (item < ID(chn[--k]));
        if (item == ID(chn[k])) _countx(chn[k], ttn_child(tree,n), min);
        else k++;               /* if the corresp. counter exists, */
        #endif                  /* count the transaction recursively, */
      }                         /* otherwise adapt the child index */
    }                           /* into the child node array */
  }
}  /* _countx() */



static int _needed (ISNODE *node)
{                               /* --- recursively check nodes */
  int    i, r;                  /* array index, check result */
  ISNODE **chn;                 /* child node array */

  assert(node);                 /* check the function argument */
  if (node->chcnt == 0) return -1; /* do not skip new leaves, */
  if (node->chcnt <= 0) return  0; /* but skip marked subtrees */
  i   = (node->offset < 0) ? node->size : PAD(node->size);
  chn = (ISNODE**)(node->cnts +node->size +i);
  for (r = 0, i = node->chcnt; --i >= 0; )
    if (chn[i]) r |= _needed(chn[i]);
  if (r) return -1;             /* recursively check all children */
  node->chcnt |= F_SKIP;        /* set the skip flag if possible */
  return 0;                     /* return 'subtree can be skipped' */
}  /* _needed() */


static int _used (ISNODE *node, int *marks, int supp)
{                               /* --- recursively check item usage */
  int    i, k, r = 0;           /* array index, map size, result */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* child node array */

  assert(node && marks);        /* check the function arguments */
  if (node->offset >= 0) {      /* if a pure array is used */
    if (node->chcnt == 0) {     /* if this is a new node (leaf) */
      k = node->offset;         /* get the index offset */
      for (i = node->size; --i >= 0; ) {
        if (node->cnts[i] >= supp)
          marks[k+i] = r = 1;   /* mark items in set that satisfy */
      } }                       /* the minimum support criterion */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      for (i = node->chcnt; --i >= 0; )
        if (chn[i]) r |= _used(chn[i], marks, supp);
    } }                         /* recursively process all children */
  else {                        /* if an identifer map is used */
    if (node->chcnt == 0) {     /* if this is a new node */
      map = node->cnts +node->size;
      for (i = node->size; --i >= 0; ) {
        if (node->cnts[i] >= supp)
          marks[map[i]] = r = 1;/* mark items in set that satisfies */
      } }                       /* the minimum support criterion */
    else if (node->chcnt > 0) { /* if there are child nodes */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      for (i = node->chcnt; --i >= 0; )
        r |= _used(chn[i], marks, supp);
    }                           /* get the child node array and */
  }                             /* recursively process all children */
  if ((r != 0) && node->parent) /* if the check succeeded, mark */
    marks[ID(node)] = 1;        /* the item associated with the node */
  return r;                     /* return the check result */
}  /* _used() */



static void _mark (ISNODE *node, int *items, int n, int supp)
{                               /* --- mark an item set */
  int    i, k;                  /* array index, map size */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* child node array */

  assert(node                   /* check the function arguments */
     && (n >= 0) && (items || (n <= 0)));
  while (--n > 0) {             /* follow the set/path from the node */
    if (node->offset >= 0) {    /* if a pure array is used */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      i   = *items++ -ID(chn[0]); }
    else {                      /* if an identifier map is used */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      i   = _search(*items++, chn, CHCNT(node));
    }                           /* get the proper child array index */
    node = chn[i];              /* go to the corresponding child */
  }
  if (node->offset >= 0)        /* if a pure array is used, */
    i   = *items -node->offset; /* compute the counter index */
  else {                        /* if an identifier map is used */
    map = node->cnts +(k = node->size);
    i   = int_bsearch(*items, map, k);
  }                             /* search for the proper index */
  if ((supp < 0)                /* if to clear unconditionally */
  ||  (node->cnts[i] == supp))  /* or the support is the same */
    node->cnts[i] |= F_SKIP;    /* mark the item set by the sign bit */
}  /* _mark() */



static void _marksub (ISTREE *ist, ISNODE *node, int index, int supp)
{                               /* --- mark all n-1 subsets */
  int i;                        /* next item, loop variable */
  int *items;                   /* (partial) item set */

  if (node->offset >= 0) i = node->offset +index;
  else                   i = node->cnts[node->size +index];
  items = ist->buf +ist->maxht; /* get and store the first two items */
  *--items = i;        _mark(node->parent, items, 1, supp);
  *--items = ID(node); _mark(node->parent, items, 1, supp);
  i = 2;                        /* mark counters in parent node */
  for (node = node->parent; node->parent; node = node->parent) {
    _mark(node->parent, items, i, supp);
    *--items = ID(node); i++;   /* climb up the tree and mark */
  }                             /* counters for all n-1 subsets */
}  /* _marksub() */



static ISNODE* _child (ISTREE *ist, ISNODE *node, int index, int spx)
{                               /* --- create child node (extend set) */
  int    i, k, n;               /* loop variables, counters */
  ISNODE *curr;                 /* to traverse the path to the root */
  int    item, cnt;             /* item identifier, number of items */
  int    *set;                  /* next (partial) item set to check */
  int    body;                  /* enough support for a rule body */
  int    hdonly;                /* whether head only item on path */
  int    app;                   /* appearance flags of an item */
  int    s_set;                 /* support of an item set */

  assert(ist && node            /* check the function arguments */
     && (index >= 0) && (index < node->size));

  /* --- initialize --- */
  s_set = node->cnts[index];    /* get support of item set to extend */
  if ((s_set <  ist->supp)      /* if the support is insufficient */
  ||  (s_set >= spx))           /* or item is a perfect extension, */
    return NULL;                /* abort (do not create a child) */
  if (node->offset >= 0) item = node->offset +index; 
  else                   item = node->cnts[node->size +index];
  app = ib_getapp(ist->base, item);  /* get item id. and app. flag */
  if ((app == APP_NONE)         /* do not extend an item to ignore */
  || ((app == APP_HEAD) && (HDONLY(node))))
    return NULL;                /* do not combine two head only items */
  hdonly = (app == APP_HEAD) || HDONLY(node);
  if ((ist->eval   >  IST_NONE) /* if to prune with evaluation */
  &&  (ist->height >= ist->prune)) {
    ist->index = index;         /* note index for aggregation */
    if (_aggregate(ist) < ist->minval) {
      node->cnts[index] |= F_SKIP; return NULL; }
  }                             /* check whether item set qualifies */
  body = (s_set >= ist->rule)   /* if the set has enough support for */
       ? 1 : 0;                 /* a rule body, set the body flag */
  ist->buf[ist->maxht-2] = item;/* init. set for support checks */

  /* --- check candidates --- */
  for (n = 0, i = index; ++i < node->size; ) {
    if (node->offset >= 0) k = node->offset +i;
    else                   k = node->cnts[node->size +i];
    app = ib_getapp(ist->base, k); /* traverse the candidate items */
    if ((app == APP_NONE) || (hdonly && (app == APP_HEAD)))
      continue;                 /* skip sets with two head only items */
    s_set = node->cnts[i];      /* traverse the candidate items */
    if ((s_set <  ist->supp)    /* if set support is insufficient */
    ||  (s_set >= spx))         /* or item is a perfect extension, */
      continue;                 /* ignore the corresponding candidate */
    body &= 1;                  /* restrict body flags to set support */
    if (s_set >= ist->rule)     /* if set support is sufficient for */
      body |= 2;                /* a rule body, set the body flag */ 
    set    = ist->buf +ist->maxht -(cnt = 2);
    set[1] = k;                 /* add the candidate item to the set */
    for (curr = node; curr->parent; curr = curr->parent) {
      s_set = _getsupp(curr->parent, set, cnt);
      if (s_set <  ist->supp)   /* get the item set support and */
        break;                  /* if it is too low, abort the loop */
      if (s_set >= ist->rule)   /* if some subset has enough support */
        body |= 4;              /* for a rule body, set the body flag */
      *--set = ID(curr); cnt++; /* add id of current node to the set */
    }                           /* and adapt the number of items */
    if (!curr->parent && body)  /* if subset support is high enough */
      ist->map[n++] = k;        /* for a full rule and a rule body, */
  }                             /* note the item identifier */
  if (n <= 0) return NULL;      /* if no child is needed, abort */
  #ifdef BENCH                  /* if benchmark version, */
  ist->scnec += n;              /* sum the necessary counters */
  #endif

  /* --- decide on node structure --- */
  k = ist->map[n-1] -ist->map[0] +1;
  if (n+n >= k) n = k;          /* use a pure array if it is small, */
  else          k = n+n;        /* otherwise use an identifier map */
  #ifdef ARCH64                 /* if 64 bit architecture, */
  if (n == k) n = k += PAD(k);  /* pad to an even number of counters */
  #endif
  #ifdef BENCH                  /* if benchmark version, */
  ist->sccnt += n;              /* sum the number of counters */
  if (n != k) ist->mapsz += n;  /* sum the size of the maps */
  ist->ndcnt++;                 /* count the node to be created */
  #endif

  /* --- create child --- */
  curr = (ISNODE*)malloc(sizeof(ISNODE) +(k-1) *sizeof(int));
  if (!curr) return (void*)-1;  /* create a child node */
  if (hdonly) item |= F_HDONLY; /* set the head only flag and */
  curr->id    = item;           /* initialize the item identifier */
  curr->chcnt = 0;              /* there are no children yet */
  curr->size  = n;              /* set size of counter array */
  if (n == k)                   /* if to use a pure array, */
    curr->offset = ist->map[0]; /* note the first item as an offset */
  else {                        /* if to use an identifier map, */
    curr->offset = -1;          /* use the offset as an indicator */
    for (set = curr->cnts +n +(i = n); --i >= 0; )
      *--set = ist->map[i];     /* copy the identifier map */
  }                             /* from the buffer to the node */
  for (set = curr->cnts +(i = n); --i >= 0; )
    *--set = 0;                 /* clear all counters of the node */
  return curr;                  /* return pointer to created child */
}  /* _child() */

static void _cleanup (ISTREE *ist)
{                               /* --- clean up on error */
  ISNODE *node, *t;             /* to traverse the nodes */

  assert(ist);                  /* check the function argument */
  for (node = ist->lvls[ist->height]; node; ) {
    t = node; node = node->succ; free(t); }
  ist->lvls[ist->height] = NULL;/* delete all created nodes */
  for (node = ist->lvls[ist->height -1]; node; node = node->succ)
    node->chcnt = 0;            /* clear the child node counters */
}  /* _cleanup() */             /* of the deepest nodes in the tree */



ISTREE* ist_create (ITEMBASE *base, int mode,
                    int supp, int smax, double conf)
{                               /* --- create an item set tree */
  int    cnt, n;                /* number of items, buffer */
  ISTREE *ist;                  /* created item set tree */
  ISNODE *root;                 /* root node of the tree */

  assert(base                   /* check the function arguments */
     && (supp >= 0) && (conf >= 0) && (conf <= 1));

  /* --- allocate memory --- */ 
  cnt = ib_cnt(base);           /* get the number of items */
  ist = (ISTREE*)malloc(sizeof(ISTREE));
  if (!ist) return NULL;        /* allocate the tree body */
  ist->lvls = (ISNODE**)malloc(BLKSIZE *sizeof(ISNODE*));
  if (!ist->lvls) {                  free(ist); return NULL; }
  ist->buf  = (int*)    malloc(BLKSIZE *sizeof(int));
  if (!ist->buf)  { free(ist->lvls); free(ist); return NULL; }
  ist->map  = (int*)    malloc(cnt *sizeof(int));
  if (!ist->map)  { free(ist->buf);
                    free(ist->lvls); free(ist); return NULL; }
  n = cnt +PAD(cnt);            /* compute the array size */
  ist->lvls[0] = ist->curr =    /* allocate a root node */
  root = (ISNODE*)calloc(1, sizeof(ISNODE) +(n-1) *sizeof(int));
  if (!root)      { free(ist->map);  free(ist->buf);
                    free(ist->lvls); free(ist); return NULL; }

  /* --- initialize structures --- */
  ist->base   = base;           /* copy parameters to the structure */
  ist->mode   = mode;
  ist->wgt    = ib_getwgt(base);
  ist->maxht  = BLKSIZE;
  ist->height = 1;
  ist->rule   = (supp > 0)         ? supp : 1;
  ist->smax   = (smax > ist->rule) ? smax : ist->rule;
  if (!(mode & APP_HEAD)) supp = (int)ceil(conf *supp);
  ist->supp   = (supp > 0)         ? supp : 1;
  ist->conf   = conf *(1.0-DBL_EPSILON);
  /* Multiplying the minimum confidence with (1.0-DBL_EPSILON) takes */
  /* care of rounding errors. For example, a minimum confidence of   */
  /* 0.8 (or 80%) cannot be coded accurately with a double precision */
  /* floating point number. It is rather stored as a value slightly  */
  /* larger number, which can cause missing rules. To prevent this,  */
  /* the confidence is made smaller by the largest possible factor.  */
  #ifdef BENCH                  /* if benchmark version */
  ist->ndcnt  = 1;   ist->ndprn = ist->mapsz = 0;
  ist->sccnt  = ist->scnec = cnt; ist->scprn = 0;
  ist->cpcnt  = ist->cpnec =      ist->cpprn = 0;
  #endif                        /* initialize the benchmark variables */
  ist_setsize(ist, 1, 1, 1);    /* init. the extraction variables */
  ist_seteval(ist, IST_NONE, IST_NONE, 1, INT_MAX);
  ist_init(ist);
  root->parent = root->succ  = NULL;
  root->offset = root->chcnt = root->id = 0;
  root->size   = n;             /* initialize the root node */
  while (--cnt >= 0)            /* copy the item frequencies */
    root->cnts[cnt] = ib_getfrq(base, cnt); 
  return ist;                   /* return created item set tree */
}  /* ist_create() */



void ist_delete (ISTREE *ist)
{                               /* --- delete an item set tree */
  int    i;                     /* loop variables */
  ISNODE *node, *t;             /* to traverse the nodes */

  assert(ist);                  /* check the function argument */
  for (i = ist->height; --i >= 0; ) {
    for (node = ist->lvls[i]; node; ) {
      t = node; node = node->succ; free(t); }
  }                             /* delete all nodes, */
  free(ist->lvls);              /* the level array, */
  free(ist->map);               /* the identifier map, */
  free(ist->buf);               /* the path buffer, */
  free(ist);                    /* and the tree body */
}  /* ist_delete() */



void ist_count (ISTREE *ist, const int *items, int n, int wgt)
{                               /* --- count a transaction */
  assert(ist                    /* check the function arguments */
     && (n >= 0) && (items || (n <= 0)));
  if (n >= ist->height)         /* recursively count the transaction */
    _count(ist->lvls[0], items, n, wgt, ist->height);
}  /* ist_count() */



void ist_countt (ISTREE *ist, const TRACT *t)
{                               /* --- count a transaction */
  int k;                        /* number of items */

  assert(ist && t);             /* check the function arguments */
  k = t_size(t);                /* get the transaction size and */
  if (k >= ist->height)         /* count the transaction recursively */
    _count(ist->lvls[0], t_items(t), k, t_wgt(t), ist->height);
}  /* ist_countt() */



void ist_countb (ISTREE *ist, const TABAG *bag)
{                               /* --- count a transaction bag */
  int   i, k;                   /* loop variable, number of items */
  TRACT *t;                     /* to traverse the transactions */

  assert(ist && bag);           /* check the function arguments */
  if (!tb_max(bag) >= ist->height)
    return;                     /* check for suff. long transactions */
  for (i = tb_cnt(bag); --i >= 0; ) {
    t = tb_tract(bag, i);       /* traverse the transactions */
    k = t_size(t);              /* get the transaction size and */
    if (k >= ist->height)       /* count the transaction recursively */
      _count(ist->lvls[0], t_items(t), k, t_wgt(t), ist->height);
  }
}  /* ist_countt() */


void ist_countx (ISTREE *ist, const TATREE *tree)
{                               /* --- count transaction in tree */
  assert(ist && tree);          /* check the function arguments */
  _countx(ist->lvls[0], tt_root(tree), ist->height);
}  /* ist_countx() */           /* recursively count the trans. tree */

/*--------------------------------------------------------------------*/

void ist_prune (ISTREE *ist)
{                               /* --- prune counters and pointers */
  int    i, k, n;               /* loop variables */
  int    *c, *map;              /* counter array, item identifier map */
  ISNODE **np, *node;           /* to traverse the nodes */
  ISNODE **chn;                 /* child node array */

  assert(ist);                  /* check the function argument */
  if (ist->height <= 1)         /* if there is only the root node, */
    return;                     /* there is nothing to prune */

  /* -- prune counters for infrequent items -- */
  for (node = ist->lvls[ist->height-1]; node; node = node->succ) {
    c = node->cnts;             /* traverse the deepest level */
    if (node->offset >= 0) {    /* if a pure array is used */
      for (n = node->size; --n >= 0; ) /* find the last */
        if (c[n] >= ist->supp) break;  /* frequent item */
      for (i = 0; i < n; i++)          /* find the first */
        if (c[i] >= ist->supp) break;  /* frequent item  */
      node->size = ++n-i;       /* set the new node size */
      #ifdef BENCH              /* if benchmark version */
      k = node->size -(n-i);    /* get the number of pruned counters */
      ist->sccnt -= k;          /* update the number of counters */
      ist->scprn += k;          /* and of pruned counters */
      #endif                    /* update the memory usage */
      if (i > 0) {              /* if there are leading infreq. items */
        node->offset += i;      /* set the new item offset */
        for (k = 0; i < n; i++) /* copy the frequent items */
          c[k++] = c[i];        /* to the front of the array */
      } }
    else {                      /* if an identifier map is used */
      map = c +node->size;      /* get the item identifier map */
      for (i = n = 0; i < node->size; i++) {
        if (c[i] >= ist->supp) {
          c[n] = c[i]; map[n++] = map[i]; }
      }                         /* remove infrequent items */
      k = node->size -n;        /* get the number of pruned counters */
      if (k <= 0) continue;     /* if no items were pruned, continue */
      #ifdef BENCH              /* if benchmark version, */
      ist->sccnt -= k;          /* update the number of counters */
      ist->scprn += k;          /* and of pruned counters */
      ist->mapsz -= k;          /* update the total item map size */
      #endif
      node->size = n; c += n;   /* set the new node size */
      for (i = 0; i < n; i++)   /* move the item identifier map */
        c[i] = map[i];          /* so that it starts directly */
    }                           /* after the support counters */
  }

  /* -- prune pointers to empty children -- */
  for (node = ist->lvls[ist->height-2]; node; node = node->succ) {
    n = CHCNT(node);            /* traverse the parent nodes */
    if (n <= 0) continue;       /* skip childless nodes */
    if (node->offset >= 0) {    /* if a pure array is used */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      while (--n >= 0)          /* find the last  non-empty child */
        if (chn[n] && (chn[n]->size > 0)) break;                
      for (i = 0; i < n; i++)   /* find the first non-empty child */
        if (chn[i] && (chn[i]->size > 0)) break;             
      node->chcnt = ++n-i;      /* set the new number of children */
      #ifdef BENCH              /* if benchmark version, */
      k = node->chcnt -(n-i);   /* get the number of pruned pointers */
      ist->cpcnt -= k;          /* update the number of pointers */
      ist->cpprn += k;          /* and of pruned pointers */
      #endif
      for (k = 0; i < n; i++)   /* remove all empty children */
        chn[k++] = (chn[i] && (chn[i]->size > 0)) ? chn[i] : NULL; }
    else {                      /* if an item identifier map is used */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      for (i = k = 0; i < n; i++)
        if (chn[i]->size > 0)   /* collect the child nodes */
          chn[k++] = chn[i];    /* that are not empty */
      node->chcnt = k;          /* set the new number of children */
      #ifdef BENCH              /* if benchmark version, */
      n -= k;                   /* get the number of pruned pointers */
      ist->cpcnt -= n;          /* update the number of pointers */
      ist->cpprn += n;          /* and of pruned pointers */
      #endif
    }                 
    if (node->chcnt <= 0)       /* if all children were removed, */
      node->chcnt |= F_SKIP;    /* set the skip flag, so that */
  }                             /* no recounting takes place */

  /* -- remove empty children -- */
  for (np = ist->lvls +ist->height-1; *np; ) {
    node = *np;                 /* traverse the deepest level again */
    if (node->size > 0) { np = &node->succ; continue; }
    *np = node->succ; free(node); /* remove empty nodes */
    #ifdef BENCH                /* if benchmark version */
    ist->ndcnt--; ist->ndprn++; /* update the number nodes */
    #endif                      /* and of pruned nodes */
  }
}  /* ist_prune() */

/*--------------------------------------------------------------------*/

int ist_check (ISTREE *ist, int *marks)
{                               /* --- check item usage */
  int i, n;                     /* loop variable, number of items */

  assert(ist);                  /* check the function argument */
  for (i = ist->lvls[0]->size; --i >= 0; )
    marks[i] = 0;               /* clear the marker array */
  _used(ist->lvls[0], marks, ist->supp);
  for (n = 0, i = ist->lvls[0]->size; --i >= 0; )
    if (marks[i]) n++;          /* count used items */
  return n;                     /* and return this number */
}  /* ist_check() */

/*--------------------------------------------------------------------*/

int ist_addlvl (ISTREE *ist)
{                               /* --- add a level to item set tree */
  int    i, n;                  /* loop variable, node counter */
  int    spx;                   /* support for a perfect extension */
  ISNODE **np;                  /* to traverse the nodes */
  ISNODE *node;                 /* current node in deepest level */
  ISNODE *par;                  /* parent of current node */
  ISNODE *cur;                  /* current node in new level (child) */
  ISNODE **frst;                /* first child of current node */
  ISNODE *last;                 /* last  child of current node */
  ISNODE **end;                 /* end of node list of new level */
  ISNODE **chn;                 /* child node array */
  void   *t;                    /* temporary buffer for reallocation */

  assert(ist);                  /* check the function arguments */

  /* --- enlarge level array --- */
  if (ist->height >= ist->maxht) {
    n = ist->maxht +BLKSIZE;    /* if the level array is full */
    t = realloc(ist->lvls, n *sizeof(ISNODE*));
    if (!t) return -1;          /* enlarge the level array */
    ist->lvls = (ISNODE**)t;    /* and set the new array */
    t = realloc(ist->buf,  n *sizeof(int));
    if (!t) return -1;          /* enlarge the buffer array */
    ist->buf   = (int*)t;       /* and set the new array */
    ist->maxht = n;             /* set the new array size */
  }                             /* (applies to buf and levels) */
  end  = ist->lvls +ist->height;
  *end = NULL;                  /* start a new tree level */

  /* --- add tree level --- */
  for (np = ist->lvls +ist->height -1; *np; np = &(*np)->succ) {
    node = ist->node = *np;     /* traverse the deepest nodes */
    frst = end; last = NULL;    /* note start of the child node list */
    if (!(ist->mode & IST_PERFECT)) spx = INT_MAX;
    else if (!node->parent)         spx = ist->wgt;
    else spx = _getsupp(node->parent, &node->id, 1);
    spx = COUNT(spx);           /* get support for perfect extension */
    for (i = n = 0; i < node->size; i++) {
      cur = _child(ist,node,i,spx);   /* traverse the counter array */
      if (!cur) continue;       /* create a child node if necessary */
      if (cur == (void*)-1) { *end = NULL; _cleanup(ist); return -1; }
      *end = last = cur;        /* add node at the end of the list */
      end  = &cur->succ; n++;   /* that contains the new level */
    }                           /* and advance the end pointer */
    if (n <= 0) {               /* if no child node was created, */
      node->chcnt = F_SKIP; continue; }         /* skip the node */
    *end = NULL;                /* terminate the child node list */
    #ifdef BENCH                /* if benchmark version, */
    ist->cpnec += n;            /* sum the number of */
    #endif                      /* necessary child pointers */
    chn = np; par = node->parent;
    if (par) {                  /* if there is a parent node */
      if (par->offset >= 0) {   /* if a pure array is used */
        chn = (ISNODE**)(par->cnts +par->size +PAD(par->size));
        chn += ID(node) -ID(chn[0]); }
      else {                    /* if an identifier map is used */
        chn = (ISNODE**)(par->cnts +par->size +    par->size);
        chn += _search(ID(node), chn, CHCNT(par));
      }                         /* find the child node pointer */
    }                           /* in the parent node */
    if (node->offset >= 0) {    /* if a pure counter array is used, */
      i = (node->size +PAD(node->size) -1) *sizeof(int);
      n = ID(last)-ID(*frst)+1;}/* always add a pure child array */
    else {                      /* if an identifier map is used */
      i = (node->size +    node->size  -1) *sizeof(int);
    }                           /* add a compact child array */
    node = (ISNODE*)realloc(node, sizeof(ISNODE) +i +n *sizeof(ISNODE));
    if (!node) { _cleanup(ist); return -1; }
    node->chcnt = n;            /* add a child array to the node */
    #ifdef BENCH                /* if benchmark version, */
    ist->cpcnt += n;            /* sum the number of child pointers */
    #endif
    *chn = *np = node;          /* set the new (reallocated) node */
    if (node->offset >= 0) {    /* if a pure array is used */
      chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
      while (--n >= 0) chn[n] = NULL;
      i   = ID(*frst);          /* get the child node array */
      for (cur = *frst; cur; cur = cur->succ) {
        chn[ID(cur)-i] = cur;   /* set the child node pointer */
        cur->parent    = node;  /* and the parent pointer */
      } }                       /* in the new node */
    else {                      /* if an identifier map is used */
      chn = (ISNODE**)(node->cnts +node->size +node->size);
      i   = 0;                  /* get the child node array */
      for (cur = *frst; cur; cur = cur->succ) {
        chn[i++]    = cur;      /* set the child node pointer */
        cur->parent = node;     /* and the parent pointer */
      }                         /* in the new node */
    }                           /* (store pointers to children */
  }                             /*  in the current node) */
  if (!ist->lvls[ist->height])  /* if no child has been added, */
    return 1;                   /* abort the function, otherwise */
  ist->height++;                /* increment the level counter */
  _needed(ist->lvls[0]);        /* mark unnecessary subtrees */
  return 0;                     /* return 'ok' */
}  /* ist_addlvl() */

/*--------------------------------------------------------------------*/

void ist_up (ISTREE *ist, int root)
{                               /* --- go up in item set tree */
  assert(ist && ist->curr);     /* check the function argument */
  if      (root)                /* if root flag set, */
    ist->curr = ist->lvls[0];   /* go to the root node */
  else if (ist->curr->parent)   /* if it exists, go to the parent */
    ist->curr = ist->curr->parent;
}  /* ist_up() */

/*--------------------------------------------------------------------*/

int ist_down (ISTREE *ist, int item)
{                               /* --- go down in item set tree */
  ISNODE *node;                 /* current node */
  ISNODE **chn;                 /* child node array */
  int    cnt;                   /* number of children */

  assert(ist && ist->curr);     /* check the function argument */
  node = ist->curr;             /* get the current node */
  cnt  = CHCNT(node);           /* if there are no child nodes, */
  if (cnt <= 0) return -1;      /* abort the function */
  if (node->offset >= 0) {      /* if a pure array is used */
    chn   = (ISNODE**)(node->cnts +node->size +PAD(node->size));
    item -= ID(chn[0]);         /* compute index in child node array */
    if ((item >= cnt) || !chn[item]) return -1; }
  else {                        /* if an identifier map is used */
    chn   = (ISNODE**)(node->cnts +node->size +    node->size);
    item  = _search(item, chn, cnt);
  }                             /* search for the proper index */
  if (item < 0) return -1;      /* if index is out of range, abort */
  ist->curr = chn[item];        /* otherwise go to the child node */
  return 0;                     /* return 'ok' */
}  /* ist_down() */

/*--------------------------------------------------------------------*/

int ist_next (ISTREE *ist, int item)
{                               /* --- get next item with a counter */
  int    i, n;                  /* array index, map size */
  int    *map;                  /* item identifier map */
  ISNODE *node;                 /* current node in tree */

  assert(ist && ist->curr);     /* check the function argument */
  node = ist->curr;             /* get the current node */
  if (node->offset >= 0) {      /* if a pure array is used, */
    i = item -node->offset;     /* compute the array index */
    if (i <  0) return node->offset;
    if (i >= node->size) return -1;
    return item +1; }           /* return the next item identifier */
  else {                        /* if an identifier map is used */
    map = node->cnts +(n = node->size);
    i = int_bsearch(item, map, n);
    i = (i < 0) ? -1-i : i+1;   /* try to find the item in the map */
    return (i < n) ? map[i] : -1;
  }                             /* return the following item */
}  /* ist_next() */

/*--------------------------------------------------------------------*/

int ist_supp (ISTREE *ist, int item)
{                               /* --- get support for an item */
  ISNODE *node;                 /* current node in tree */

  assert(ist && ist->curr);     /* check the function argument */
  node = ist->curr;             /* get the current node */
  if (node->offset >= 0) {      /* if pure arrays are used, */
    item -= node->offset;       /* get index in counter array */
    if (item >= node->size) return 0; }
  else                          /* if an identifier map is used */
    item = int_bsearch(item, node->cnts +node->size, node->size);
  if (item < 0) return 0;       /* abort if index is out of range */
  return COUNT(node->cnts[item]);
}  /* ist_supp() */             /* return the item set support */

/*--------------------------------------------------------------------*/

int ist_suppx (ISTREE *ist, int *items, int n)
{                               /* --- get support of an item set */
  assert(ist                    /* check the function arguments */
     && (n >= 0) && (items || (n <= 0)));
  if (n <= 0)                   /* if the item set is empty, */
    return COUNT(ist->wgt);     /* return the total trans. weight */
  return COUNT(_getsupp(ist->lvls[0], items, n));
}  /* ist_suppx() */            /* return the item set support */

/*--------------------------------------------------------------------*/

void ist_mark (ISTREE *ist, int mode)
{                               /* --- mark frequent item sets */
  int    i, k;                  /* loop variables */
  ISNODE *node;                 /* to traverse the nodes */
  int    supp;                  /* support of an item set */

  assert(ist);                  /* check the function argument */
  if ((mode & ~IST_EVAL) == IST_CLEAR) {
    ist->wgt &= ~F_SKIP;        /* if to clear all skip flags */
    for (k = 0; ++k < ist->height; )
      for (node = ist->lvls[k]; node; node = node->succ)
        for (i = node->size; --i >= 0; )
          node->cnts[i] &= ~F_SKIP;
    return;                     /* clear skip flags of all sets */
  }                             /* and abort the function */
  if ((mode      & IST_EVAL)    /* if maximal sets w.r.t. evaluation */
  &&  (ist->eval > IST_NONE)) { /* and evaluation measure is given */
    supp = -1;                  /* set default support filter */
    for (k = ist->height; --k > 0; ) {
      for (node = ist->lvls[k]; node; node = node->succ) {
        ist->node = node;       /* traverse the node of each level */
        for (i = node->size; --i >= 0; ) {
          ist->index = i;       /* traverse the items in each node */
          if ((node->cnts[i]   >=  0)
          && ((node->cnts[i]   <  ist->supp)
          ||  (_aggregate(ist) <  ist->minval))) {
            node->cnts[i] |= F_SKIP; continue; }
          if (mode & IST_CLOSED) supp = node->cnts[i];
          _marksub(ist, node, i, supp);
	}                       /* traverse the tree bottom up */
      }                         /* (that is, from leaves to root) */
    }                           /* and mark relevant n-1 subsets */
    supp = (mode & IST_CLOSED) ? ist->wgt : ist->supp;
    node = ist->lvls[0];        /* traverse the root node elements */
    for (i = node->size; --i >= 0; ) {
      if ((node->cnts[i] >= supp) || (node->cnts[i] < 0)) {
        ist->wgt |= F_SKIP; break; }
    } }                         /* mark the empty set if necessary */
  else {                        /* if to use only the support */
    supp = (mode & IST_CLOSED) ? ist->wgt : ist->supp;
    node = ist->lvls[0];        /* traverse the root node elements */
    for (i = node->size; --i >= 0; ) {
      if (node->cnts[i] >= supp) {
        ist->wgt |= F_SKIP; break; }
    }                           /* mark the empty set if necessary */
    supp = -1;                  /* set default support filter */
    for (k = 0; ++k < ist->height; ) {
      for (node = ist->lvls[k]; node; node = node->succ) {
        for (i = node->size; --i >= 0; ) {
          if (node->cnts[i] < ist->supp) {
            node->cnts[i] |= F_SKIP; continue; }
          if (mode == IST_CLOSED) supp = node->cnts[i];
          _marksub(ist, node, i, supp);
        }                       /* traverse the tree top down */
      }                         /* (that is, from root to leaves) */
    }                           /* and mark all n-1 subsets */
  }                             /* of the current item set */
}  /* ist_mark() */

/*--------------------------------------------------------------------*/

void ist_setsize (ISTREE *ist, int min, int max, int dir)
{                               /* --- set the set/rule size range */
  assert(ist);                  /* check the function arguments */
  ist->maxsz = (max < 0) ? -1 : max;  /* store the size range */
  ist->minsz = (min < 0) ?  0 : min;  /* (min. and max. size) and */
  ist->dir   = (dir < 0) ? -1 : 1;    /* the traversal direction */
}  /* ist_setsize() */

/*--------------------------------------------------------------------*/

void ist_seteval (ISTREE *ist, int eval, int agg, double minval,
                  int prune)
{                               /* --- set additional evaluation */
  assert(ist);                  /* check the function arguments */
  ist->eval   = ((eval > IST_NONE) && (eval < sizeof(_evalfns)))
              ? eval : IST_NONE;/* check and note the eval. measure */
  ist->agg    = ((agg  > IST_NONE) && (agg  < sizeof(_aggrfns)))
              ? agg  : IST_NONE;/* check and note the agg. mode */
  ist->minval = minval *(1.0-DBL_EPSILON);
  ist->prune  = (prune > 0) ? prune : INT_MAX;
}  /* ist_seteval() */          /* note the parameters */

/*--------------------------------------------------------------------*/

void ist_init (ISTREE *ist)
{                               /* --- initialize (rule) extraction */
  assert(ist);                  /* check the function argument */
  if (ist->maxsz > ist->height)  ist->maxsz = ist->height;
  ist->size  = (ist->dir >= 0) ? ist->minsz : ist->maxsz;
  ist->node  = ist->lvls[(ist->size > 0) ? ist->size -1 : 0];
  ist->index = ist->item = -1;  /* initialize the */
  ist->head  = NULL;            /* extraction variables */
}  /* ist_init() */

/*--------------------------------------------------------------------*/

static int _emptyset (ISTREE *ist, int *supp, double *eval)
{                               /* --- whether to report empty set */
  assert(ist);                  /* check the function argument */
  ist->size += ist->dir;        /* immediately go the next level */
  if ((ist->wgt >= ist->supp)   /* if the empty set qualifies */
  &&  (ist->wgt <= ist->smax)   /* (w.r.t. support and evaluation) */
  && ((ist->eval == IST_NONE) || (ist->minval <= 0))) {
    if (supp) *supp = COUNT(ist->wgt);
    if (eval) *eval = 0;        /* store the empty item set support */
    return -1;                  /* the add. evaluation is always 0 */
  }                             /* return 'report the empty set' */
  return 0;                     /* return 'do not report' */
}  /* _empty() */

/*--------------------------------------------------------------------*/

int ist_set (ISTREE *ist, int *set, int *supp, double *eval)
{                               /* --- extract next frequent item set */
  int    i;                     /* loop variable */
  int    item;                  /* an item identifier */
  ISNODE *node;                 /* current item set node */
  int    s_set;                 /* support of the current set */
  double val;                   /* value of evaluation measure */

  assert(ist && set);           /* check the function arguments */
  if ((ist->size < ist->minsz)  /* if below the minimal size */
  ||  (ist->size > ist->maxsz)) /* or above the maximal size, */
    return -1;                  /* abort the function */
  if ((ist->size == 0)          /* if to report the empty item set */
  &&  _emptyset(ist, supp, eval))
    return  0;                  /* check whether it qualifies */

  /* --- find frequent item set --- */
  node = ist->node;             /* get the current item set node */
  while (1) {                   /* search for a frequent item set */
    if (++ist->index >= node->size) { /* if all subsets have been */
      node = node->succ;        /* processed, go to the successor */
      if (!node) {              /* if at the end of a level, */
        ist->size += ist->dir;  /* go to the next level */
        if ((ist->size < ist->minsz)
        ||  (ist->size > ist->maxsz))
          return -1;            /* if outside size range, abort */
        if ((ist->size == 0)    /* if to report the empty item set */
        &&  _emptyset(ist, supp, eval))
          return  0;            /* check whether it qualifies */
        node = ist->lvls[ist->size -1];
      }                         /* get the 1st node of the new level */
      ist->node  = node;        /* note the new item set node */
      ist->index = 0;           /* start with the first item set */
    }                           /* of the new item set node */
    if (node->offset >= 0) item = node->offset +ist->index;
    else                   item = node->cnts[node->size +ist->index];
    if (ib_getapp(ist->base, item) == APP_NONE)
      continue;                 /* skip items to ignore */
    s_set = node->cnts[ist->index];
    if ((s_set < ist->supp)     /* if the support is not sufficient */
    ||  (s_set > ist->smax))    /* or larger than the maximum, */
      continue;                 /* go to the next item set */
    /* Note that this check automatically skips all item sets that */
    /* are marked with the flag F_SKIP, because s_set is negative  */
    /* with this flag and thus necessarily smaller than ist->supp. */
    if (ist->eval <= IST_NONE){ /* if no add. eval. measure given */
      val = 0; break; }         /* abort the loop (select the set) */
    val = (ist->eval >= IST_LOGQ) ? _logq(ist) : _aggregate(ist);
    if (val >= ist->minval)     /* if the evaluation is high enough, */
      break;                    /* abort the loop (select the set) */
  }  /* while (1) */
  if (supp) *supp = s_set;      /* store the item set support and */
  if (eval) *eval = val;        /* the value of the add. measure */

  /* --- build frequent item set --- */
  i        = ist->size;         /* get the current item set size */
  set[--i] = item;              /* and store the first item */
  while (node->parent) {        /* while not at the root node */
    set[--i] = ID(node);        /* add item to the item set */
    node = node->parent;        /* and go to the parent node */
  }
  return ist->size;             /* return the item set size */
}  /* ist_set() */

/*--------------------------------------------------------------------*/

int ist_rule (ISTREE *ist, int *rule,
              int *supp, int *body, int *head, double *eval)
{                               /* --- extract next association rule */
  int    i;                     /* loop variable */
  int    item;                  /* an item identifier */
  ISNODE *node;                 /* current item set node */
  ISNODE *parent;               /* parent of the item set node */
  int    *map, n;               /* identifier map and its size */
  int    s_set;                 /* support of set  (body & head) */
  int    s_body;                /* support of body (antecedent) */
  int    s_head;                /* support of head (consequent) */
  double val;                   /* value of evaluation measure */
  int    app;                   /* appearance flag of head item */

  assert(ist && rule);          /* check the function arguments */
  if (ist->size == 0)           /* if at the empty item set, */
    ist->size += ist->dir;      /* go to the next item set size */
  if ((ist->size < ist->minsz)  /* if the item set is too small */
  ||  (ist->size > ist->maxsz)) /* or too large (number of items), */
    return -1;                  /* abort the function */

  /* --- find rule --- */
  node = ist->node;             /* get the current item set node */
  while (1) {                   /* search for a rule */
    if (ist->item >= 0) {       /* --- select next item subset */
      *--ist->path = ist->item; /* add previous head to the path */
      ist->item = ID(ist->head);/* and get the next head item */
      ist->head = ist->head->parent;
      if (!ist->head)           /* if all subsets have been processed */
        ist->item = -1;         /* clear the head item to trigger the */
    }                           /* selection of a new item set */
    if (ist->item < 0) {        /* --- select next item set */
      if (++ist->index >= node->size){/* if all subsets have been */
        node = node->succ;      /* processed, go to the successor */
        if (!node) {            /* if at the end of a level, go down */
          ist->size += ist->dir;/* go to the next level */
          if ((ist->size < ist->minsz) || (ist->size <= 0)
          ||  (ist->size > ist->maxsz))
            return -1;          /* if outside the size range, abort */
          node = ist->lvls[ist->size -1];
        }                       /* get the 1st node of the new level */
        ist->node = node;       /* note the new item set node and */
        ist->index  = 0;        /* start with the first item set */
      }                         /* of the new item set node */
      if (node->offset >= 0) item = node->offset +ist->index;
      else                   item = node->cnts[node->size +ist->index];
      app = ib_getapp(ist->base, item);
      if ((app == APP_NONE) || ((app == APP_HEAD) && HDONLY(node)))
        continue;               /* skip sets with two head only items */
      ist->item   = item;       /* set the head item identifier */
      ist->hdonly = (app == APP_HEAD) || HDONLY(node);
      ist->head   = node;       /* set the new head item node */
      ist->path   = ist->buf +ist->maxht;
    }                           /* clear the path (reinitialize it) */
    app = ib_getapp(ist->base, ist->item);
    if (!(app &  APP_HEAD)      /* get head item appearance indicator */
    ||  ((app != APP_HEAD) && ist->hdonly))
      continue;                 /* if rule is not allowed, skip it */
    s_set = COUNT(node->cnts[ist->index]);
    if ((s_set < ist->supp)     /* if the support is not sufficient */
    ||  (s_set > ist->smax)) {  /* or larger than the maximum, */
      ist->item = -1; continue; }   /* go to the next item set */
    parent = node->parent;      /* get the parent node */
    n = (int)(ist->buf +ist->maxht -ist->path);
    if (n > 0)                  /* if there is a path, use it */
      s_body = COUNT(_getsupp(ist->head, ist->path, n));
    else if (!parent)           /* if there is no parent (root node), */
      s_body = COUNT(ist->wgt); /* get the total trans. weight */
    else if (parent->offset >= 0)  /* if a pure array is used */
      s_body = COUNT(parent->cnts[ID(node) -parent->offset]);
    else {                      /* if an identifier map is used */
      map = parent->cnts +(n = parent->size);
      s_body = COUNT(parent->cnts[int_bsearch(ID(node), map, n)]);
    }                           /* find array index and get support */
    if ((s_body < ist->rule)    /* if the body support is too low */
    ||  (s_set  < s_body *ist->conf))       /* or the confidence, */
      continue;                 /* go to the next item (sub)set */
    s_head = COUNT(ist->lvls[0]->cnts[ist->item]);
    if ((ist->eval <= IST_NONE) /* if no add. eval. measure given */
    ||  !_evalfns[ist->eval]) { /* or the measure does not exist, */
      val = 0; break; }         /* abort the loop (select the rule) */
    val = _evalfns[ist->eval](s_set, s_body, s_head, COUNT(ist->wgt));
    if (val >= ist->minval)     /* if the evaluation is high enough, */
      break;                    /* abort the loop (select the rule) */
  }  /* while (1) */
  if (supp) *supp = s_set;      /* store the rule support values */
  if (body) *body = s_body;     /* (whole rule and only body) */
  if (head) *head = s_head;     /* store the head item support, */
  if (eval) *eval = val;        /* the value of the add. measure */

  /* --- build rule --- */
  if (node->offset >= 0) item = node->offset +ist->index;
  else                   item = node->cnts[node->size +ist->index];
  i = ist->size;                /* get the current item and */
  if (item != ist->item)        /* if this item is not the head, */
    rule[--i] = item;           /* add it to the rule body */
  while (node->parent) {        /* traverse the path to the root */
    if (ID(node) != ist->item)  /* and add all items on this */
      rule[--i] = ID(node);     /* path to the rule body */
    node = node->parent;        /* (except the head of the rule) */
  }
  rule[0] = ist->item;          /* set the head of the rule, */
  return ist->size;             /* return the rule size */
}  /* ist_rule() */

/*--------------------------------------------------------------------*/

static int _report (ISTREE *ist, ISREPORT *rep, ISNODE *node, int supp)
{                               /* --- recursive item set reporting */
  int    i, k, c, n = 0;        /* loop variable, set counter */
  int    spx;                   /* support for perfect extension */
  int    off;                   /* item offset */
  int    *map;                  /* item identifier map */
  ISNODE **chn;                 /* child node array */

  assert(ist && rep);           /* check the function arguments */
  if (!(ist->mode & IST_PERFECT))  /* if no perfext extension pruning */
    spx = INT_MAX;              /* clear perfect extension support */
  else {                        /* if perfect extensions pruning */
    spx = supp;                 /* note the parent set support */
    for (i = 0; i < node->size; i++) {
      if (COUNT(node->cnts[i]) < spx) continue;
      if (node->offset >= 0) k = i +node->offset;
      else k = node->cnts[node->size +i];
      isr_addpex(rep, k);       /* traverse the node's items and */
    }                           /* collect the perfect extensions */
  }                             /* in the item set reporter */
  if ((supp >= 0)               /* if current item set is not marked */
  &&  (supp <= ist->smax))      /* and does not exceed max. support, */
    n += isr_report(rep);       /* report the current item set */
  if (node->offset >= 0) {      /* if a pure array is used */
    chn = (ISNODE**)(node->cnts +node->size +PAD(node->size));
    c   = CHCNT(node);          /* get the child node array */
    off = (c > 0) ? ID(chn[0]) : 0;
    for (i = 0; i < node->size; i++) {
      supp = COUNT(node->cnts[i]);
      if ((supp <  ist->supp)   /* traverse the node's items and */
      ||  (supp >= spx))        /* check against minimum support */
        continue;               /* and the parent set support */
      ist->node  = node;        /* store the node and the index */
      ist->index = i;           /* in the node for evaluation */
      k = node->offset +i;      /* compute the item identifier */
      isr_add(rep, k, supp);    /* add the item to the reporter */
      supp = node->cnts[i];     /* get the item support (with flag) */
      k -= off;                 /* compute the child node index */
      if ((k >= 0)              /* if the corresp. child node exists, */
      &&  (k <  c) && chn[k])   /* recursively report the subtree */
        n += _report(ist, rep, chn[k], supp);
      else if ((supp >= 0)      /* if the item set is not marked */
      &&       (supp <= ist->smax))
        n += isr_report(rep);   /* report the current item set */
      isr_remove(rep, 1);       /* remove the last item */
    } }                         /* from the current item set */
  else {                        /* if an identifier map is used */
    map = node->cnts +(k = node->size);
    chn = (ISNODE**)(map +k +PAD(k)); /* get the item id map */
    c   = CHCNT(node);          /* and the child node array  */
    c   = (c > 0) ? ID(chn[c-1]) : -1;
    for (i = 0; i < node->size; i++) {
      supp = COUNT(node->cnts[i]);
      if ((supp <  ist->supp)   /* traverse the node's items and */
      ||  (supp >= spx))        /* check against minimum support */
        continue;               /* and the parent set support */
      ist->node  = node;        /* store the node and the index */
      ist->index = i;           /* in the node for evaluation */
      k = map[i];               /* retrieve the item identifier */
      isr_add(rep, k, supp);    /* add the item to the reporter */
      supp = node->cnts[i];     /* get the item support (with flag) */
      if (k <= c)               /* if there may be a child node, */
        while (k > ID(*chn)) chn++;  /* skip the preceding items */
      if ((k <= c)              /* if the corresp. child node exists, */
      &&  (k == ID(*chn)))      /* recursively report the subtree */
        n += _report(ist, rep, *chn, supp);
      else if ((supp >= 0)      /* if the item set is not marked */
      &&       (supp <= ist->smax))
        n += isr_report(rep);   /* report the current item set */
      isr_remove(rep, 1);       /* remove the last item */
    }                           /* from the current item set */
  }
  return n;                     /* return the number of item sets */
}  /* _report() */

/*--------------------------------------------------------------------*/

int ist_report (ISTREE *ist, ISREPORT *rep)
{                               /* --- recursive item set reporting */
  assert(ist && rep);           /* check the function arguments */
  return _report(ist, rep, ist->lvls[0], ist->wgt);
}  /* ist_report() */           /* recursively report item sets */

/*--------------------------------------------------------------------*/

double ist_eval (ISTREE *ist)
{                               /* --- evaluate current item set */
  assert(ist);                  /* check the function argument */
  if (ist->eval <= IST_NONE) return 0;
  if (ist->eval >= IST_LOGQ) return _logq(ist);
  return _aggregate(ist);       /* compute the evaluation measure */
}  /* ist_eval() */

/*--------------------------------------------------------------------*/

double ist_evalx (ISREPORT *rep, void *data)
{                               /* --- evaluate current item set */
  ISTREE *ist;                  /* item set tree to work on */

  assert(rep && data);          /* check the function arguments */
  ist = (ISTREE*)data;          /* type the user data */
  if (ist->eval <= IST_NONE) return 0;
  if (ist->eval >= IST_LOGQ) return _logq(ist);
  return _aggregate(ist);       /* compute the evaluation measure */
}  /* ist_evalx() */

/*--------------------------------------------------------------------*/
#ifndef NDEBUG

static void _showtree (ISNODE *node, ITEMBASE *base, int level)
{                               /* --- show subtree */
  int    i, k, cnt;             /* loop variables, number of children */
  ISNODE **chn;                 /* child node array */

  assert(node && (level >= 0)); /* check the function arguments */
  i   = (node->offset < 0) ? node->size : PAD(node->size);
  chn = (ISNODE**)(node->cnts +node->size +i);
  cnt = CHCNT(node);            /* get the child node array */
  for (i = 0; i < node->size; i++) {
    for (k = level; --k >= 0; ) /* indent and print */
      printf("   ");            /* item identifier and counter */
    if (node->offset >= 0) k = node->offset +i;
    else                   k = node->cnts[node->size +i];
    printf("%s: %d\n", ib_name(base, k), COUNT(node->cnts[i]));
    if (cnt <= 0) continue;     /* check whether there are children */
    if (node->offset >= 0) k -= ID(chn[0]);
    else                   k  = _search(k, chn, cnt);
    if ((k >= 0) && (k < cnt) && chn[k])
      _showtree(chn[k], base, level +1);
  }                             /* show subtree recursively */
}  /* _showtree() */

/*--------------------------------------------------------------------*/

void ist_show (ISTREE *ist)
{                               /* --- show an item set tree */
  assert(ist);                  /* check the function argument */
  _showtree(ist->lvls[0], ist->base, 0);
  printf("total: %d\n", COUNT(ist->wgt));
}  /* ist_show() */             /* show the nodes recursively */

#endif
