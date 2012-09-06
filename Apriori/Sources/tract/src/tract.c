#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "tract.h"
#include "scan.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define BLKSIZE     256         /* block size for enlarging arrays */

#ifdef ARCH64
#define CHOFF(n)    ((n) +1 -((n) & 1))
#else
#define CHOFF(n)    (n)
#endif
/* The item array offset to the child node array must be odd for */
/* a 64 bit architecture, because this offset +3 must be an even */
/* number, as there are three int fields before the items array. */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
static int   sorted   = 0;      /* flag for sorted indicators */
static char* appmap[] = {       /* item appearance indicators */
  "0:-",  "0:none",  "0:neither", "0:ignore",
  "1:i",  "1:in",    "1:a",  "1:antecedent", "1:b", "1:body",
  "2:o",  "2:out",   "2:c",  "2:consequent", "2:h", "2:head",
  "3:io", "3:inout", "3:ac", "3:a&c", "3:both", "3:bh", "3:b&h",
};                              /* (code:name) */

/*----------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------*/

static int _appcmp (const void *p1, const void *p2, void *data)
{                               /* --- compare appearance indicators */
  const char *s1 = (const char*)p1 +2; /* type the two pointers */
  const char *s2 = (const char*)p2 +2; /* to strings (skip code) */
  register int d;               /* difference of character values */

  for ( ; 1; s1++, s2++) {      /* traverse the strings */
    d = *s1 -*s2;               /* if one string is smaller than */
    if (d   != 0) return d;     /* the other, return the result */
    if (*s1 == 0) return 0;     /* if at the end of the string, */
  }                             /* return 'strings are equal' */
}  /* _appcmp() */

/*--------------------------------------------------------------------*/

static int _appcode (const char *s)
{                               /* --- get appearance indicator code */
  int  i, n;                    /* index, number of app. indicators */
  char *t;                      /* to traverse the indicator name */

  assert(s);                    /* check the function argument */
  n = sizeof(appmap)/sizeof(char*);
  if (!sorted) {                /* if app. indicators not sorted yet */
    ptr_qsort(appmap, n, _appcmp, NULL);
    sorted = -1;                /* sort the appearance indicators */
  }                             /* and set the sorted flag */
  i = ptr_bsearch(s-2, appmap, n, _appcmp, NULL);
  if (i >= n) return -1;        /* try to find appearance indicator */
  if (i >= 0)                   /* if an exact match was found, */
    return appmap[i][0] -'0';   /* return the app. indicator code */
  t = appmap[i = -1-i] +2;      /* otherwise check for a prefix match */
  while (*s && (*s == *t)) s++, t++;
  if (*s) return -1;            /* if no match, abort the function */
  return appmap[i][0] -'0';     /* return the app. indicator code */
}  /* _appcode() */

/*--------------------------------------------------------------------*/

static int _read (ITEMBASE *base, FILE *file)
{                               /* --- read an item */
  int   d, n;                   /* delimiter type, array size */
  char  *buf;                   /* read buffer */
  ITEM  *item;                  /* item corresponding to read name */
  TRACT *t;                     /* to access the transaction buffer */

  assert(base && file);         /* check the function arguments */
  d = ts_next(base->tscan, file, NULL, 0);
  if (d == TS_ERR)    return d; /* read the next field (item name) */
  buf = ts_buf(base->tscan);    /* get the read buffer and */
  if (buf[0] == '\0') return d; /* check whether it is empty */
  item = nim_byname(base->nimap, buf);
  if (!item) {                  /* look up the name in name/id map */
    if (base->app == APP_NONE)  /* if new items are to be ignored, */
      return d;                 /* do not register the item */
    item = nim_add(base->nimap, buf, sizeof(ITEM));
    if (!item) return E_NOMEM;  /* add the new item to the map, */
    item->frq = item->xfq = 0;  /* initialize the frequency counters */
    item->app = base->app;      /* (occurrence and sum of t.a. sizes) */
    item->pen = base->pen;      /* set the appearance indicator */
  }                             /* and the insertion penalty */
  t = base->tract;              /* get the transaction buffer */
  n = base->size;               /* and its current size */
  if (t->size >= n) {           /* if the transaction buffer is full */
    n += (n > BLKSIZE) ? (n >> 1) : BLKSIZE;
    t  = (TRACT*)realloc(t, sizeof(TRACT) +n *sizeof(int));
    if (!t) return E_NOMEM;     /* enlarge the transaction buffer */
    base->tract = t; base->size = n;
  }                             /* set the new buffer and its size */
  t->items[t->size++] = item->id;
  return d;                     /* add the item to the transaction */
}  /* _read() */                /* and return the delimiter type */

/*--------------------------------------------------------------------*/

static int _nocmp (const void *p1, const void *p2, void *data)
{                               /* --- compare item frequencies */
  const ITEM *a = p1, *b = p2;  /* type the item pointers */
  if (a->app == APP_NONE)   return (b->app == APP_NONE)   ? 0 : 1;
  if (b->app == APP_NONE)   return -1;
  if (a->frq < *(int*)data) return (b->frq < *(int*)data) ? 0 : 1;
  if (b->frq < *(int*)data) return -1;
  return a->id -b->id;          /* return sign of identifier diff. */
}  /* _nocmp() */

/*--------------------------------------------------------------------*/

static int _asccmp (const void *p1, const void *p2, void *data)
{                               /* --- compare item frequencies */
  const ITEM *a = p1, *b = p2;  /* type the item pointers */
  if (a->app == APP_NONE)   return (b->app == APP_NONE)   ? 0 : 1;
  if (b->app == APP_NONE)   return -1;
  if (a->frq < *(int*)data) return (b->frq < *(int*)data) ? 0 : 1;
  if (b->frq < *(int*)data) return -1;
  return a->frq -b->frq;        /* return sign of frequency diff. */
}  /* _asccmp() */

/*--------------------------------------------------------------------*/

static int _descmp (const void *p1, const void *p2, void *data)
{                               /* --- compare item frequencies */
  const ITEM *a = p1, *b = p2;  /* type the item pointers */
  if (a->app == APP_NONE)   return (b->app == APP_NONE)   ? 0 : 1;
  if (b->app == APP_NONE)   return -1;
  if (a->frq < *(int*)data) return (b->frq < *(int*)data) ? 0 : 1;
  if (b->frq < *(int*)data) return -1;
  return b->frq -a->frq;        /* return sign of frequency diff. */
}  /* _descmp() */

/*--------------------------------------------------------------------*/

static int _asccmpx (const void *p1, const void *p2, void *data)
{                               /* --- compare item frequencies */
  const ITEM *a = p1, *b = p2;  /* type the item pointers */
  if (a->app == APP_NONE)   return (b->app == APP_NONE)   ? 0 : 1;
  if (b->app == APP_NONE)   return -1;
  if (a->frq < *(int*)data) return (b->frq < *(int*)data) ? 0 : 1;
  if (b->frq < *(int*)data) return -1;
  return a->xfq -b->xfq;        /* return sign of frequency diff. */
}  /* _asccmpx() */

/*--------------------------------------------------------------------*/

static int _descmpx (const void *p1, const void *p2, void *data)
{                               /* --- compare item frequencies */
  const ITEM *a = p1, *b = p2;  /* type the pointers */
  if (a->app == APP_NONE)   return (b->app == APP_NONE)   ? 0 : 1;
  if (b->app == APP_NONE)   return -1;
  if (a->frq < *(int*)data) return (b->frq < *(int*)data) ? 0 : 1;
  if (b->frq < *(int*)data) return -1;
  return b->xfq -a->xfq;        /* return sign of frequency diff. */
}  /* _descmpx() */

/*----------------------------------------------------------------------
  Item Base Functions
----------------------------------------------------------------------*/

ITEMBASE* ib_create (int size)
{                               /* --- create an item base */
  ITEMBASE *base;               /* created item base */

  if (size <= 0) size = BLKSIZE;/* check and adapt number of items */
  base = malloc(sizeof(ITEMBASE));
  if (!base) return NULL;       /* create an item base */
  base->tscan = ts_create();    /* and its components */
  base->tract = (TRACT*)malloc(sizeof(TRACT) +size *sizeof(int));
  base->nimap = nim_create(0, 0, (HASHFN*)0, (OBJFN*)0);
  if (!base->tscan || !base->tract || !base->nimap) {
    ib_delete(base); return NULL; }
  base->wgt   = 0;              /* initialize the fields */
  base->app   = APP_BOTH;
  base->pen   = 0.0;
  base->size  = size;
  ts_chars(base->tscan, TS_NULL, "");
  base->chars[0] = ' ';         /* blank */
  base->chars[1] = ' ';         /* field  separator */
  base->chars[2] = '\n';        /* record separator */
  base->chars[3] = '\0';        /* dummy to terminate the string */
  base->tract->size     =  0;   /* clear the transaction buffer */
  base->tract->wgt      =  0;
  base->tract->items[0] = -1;
  return base;                  /* return the created item set */
}  /* ib_create() */

/*--------------------------------------------------------------------*/

void ib_delete (ITEMBASE *base)
{                               /* --- delete an item set */
  assert(base);                 /* check the function argument */
  if (base->nimap) nim_delete(base->nimap);
  if (base->tract) t_delete  (base->tract);
  if (base->tscan) ts_delete (base->tscan);
  free(base);                   /* delete the components */
}  /* ib_delete() */            /* and the item set body */

/*--------------------------------------------------------------------*/

void ib_chars (ITEMBASE *base,
               const char *blanks,  const char *fldseps,
               const char *recseps, const char *comment)
{                               /* --- set special characters */
  assert(base);                 /* check the function argument */
  if (blanks)                   /* set blank characters */
    base->chars[0] = ts_chars(base->tscan, TS_BLANK,  blanks);
  if (fldseps)                  /* set field separators */
    base->chars[1] = ts_chars(base->tscan, TS_FLDSEP, fldseps);
  if (recseps)                  /* set record separators */
    base->chars[2] = ts_chars(base->tscan, TS_RECSEP, recseps);
  if (comment)                  /* set comment indicators */
    ts_chars(base->tscan, TS_COMMENT, comment);
}  /* ib_chars() */

/*--------------------------------------------------------------------*/

int ib_add (ITEMBASE *base, const char *name)
{                               /* --- add an item to the set */
  ITEM *item;                   /* to access the item data */

  assert(base && name);         /* check the function arguments */
  item = nim_add(base->nimap, name, sizeof(ITEM));
  if (item == NULL)   return -1;/* add the new item */
  if (item == EXISTS) return -2;/* to the name/id map */
  return item->id;              /* return the item identifier */
}  /* ib_add() */

/*--------------------------------------------------------------------*/

int ib_readapp (ITEMBASE *base, FILE *file)
{                               /* --- read appearance indicators */
  int  d;                       /* delimiter type */
  char *buf;                    /* read buffer */
  ITEM *item;                   /* to access the item data */

  assert(base && file);         /* check the function arguments */
  buf = ts_buf(base->tscan);    /* read the first record (one field) */
  d   = ts_next(base->tscan, file, NULL, 0);
  if (d == TS_ERR)      return E_FREAD;
  if (d != TS_REC)      return E_FLDCNT;
  base->app = _appcode(buf);    /* get default appearance code */
  if (base->app < 0)    return E_UNKAPP;
  while (1) {                   /* read item/indicator pairs */
    d = ts_next(base->tscan, file, NULL, 0);
    if (d <= TS_EOF)            /* read the next item */
      return (d == TS_ERR) ? E_FREAD : 0;
    if (buf[0] == '\0')         /* check for end of file */
      return E_ITEMEXP;         /* and for a missing item */
    item = nim_add(base->nimap, buf, sizeof(ITEM));
    if (item == NULL)   return E_NOMEM;    /* add the new item */
    if (item == EXISTS) return E_DUPITEM;  /* to the name/id map */
    item->frq = item->xfq = 0;  /* clear the frequency counters */
    item->pen = base->pen;      /* and the insertion penalty */
    if (d != TS_FLD)    return E_APPEXP;
    d = ts_next(base->tscan, file, NULL, 0);
    if (d == TS_ERR)    return E_FREAD;
    if (d == TS_FLD)    return E_FLDCNT;
    item->app = _appcode(buf);  /* get the appearance indicator */
    if (item->app <  0) return E_UNKAPP;
  }                             /* (abort with error code on failure) */
  return 0;                     /* return 'ok' */
}  /* ib_readapp() */

/*--------------------------------------------------------------------*/

int ib_readpen (ITEMBASE *base, FILE *file)
{                               /* --- read insertion penalties */
  int  d;                       /* delimiter type */
  char *buf, *s;                /* read buffer, end pointer */
  ITEM *item;                   /* to access the item data */

  assert(base && file);         /* check the function arguments */
  buf = ts_buf(base->tscan);    /* read the first record (one field) */
  d   = ts_next(base->tscan, file, NULL, 0);
  if (d == TS_ERR)      return E_FREAD;
  if (d != TS_REC)      return E_FLDCNT;
  base->pen = strtod(buf, &s);  /* get default appearance code */
  if (*s || (base->pen < 0) || (base->pen > 1)) return E_PENALTY;
  while (1) {                   /* read item/indicator pairs */
    d = ts_next(base->tscan, file, NULL, 0);
    if (d <= TS_EOF)            /* read the next item */
      return (d == TS_ERR) ? E_FREAD : 0;
    if (buf[0] == '\0')         /* check for end of file */
      return E_ITEMEXP;         /* and for a missing item */
    item = nim_add(base->nimap, buf, sizeof(ITEM));
    if (item == NULL)   return E_NOMEM;    /* add the new item */
    if (item == EXISTS) return E_DUPITEM;  /* to the name/id map */
    item->frq = item->xfq = 0;  /* clear the frequency counters */
    item->app = base->app;      /* and the appearance indicator */
    if (d != TS_FLD)    return E_PENEXP;
    d = ts_next(base->tscan, file, NULL, 0);
    if (d == TS_ERR)    return E_FREAD;
    if (d == TS_FLD)    return E_FLDCNT;
    item->pen = strtod(buf,&s); /* get the insertion penalty */
    if (*s || (item->pen < 0) || (item->pen > 1)) return E_PENALTY;
  }                             /* (abort with error code on failure) */
  return 0;                     /* return 'ok' */
}  /* ib_readpen() */

/*--------------------------------------------------------------------*/

int ib_read (ITEMBASE *base, FILE *file)
{                               /* --- read a transaction */
  int   i, d, x;                /* loop variable, delimiter, buffer */
  char  *buf;                   /* read buffer */
  ITEM  *item;                  /* pointer to an item */
  TRACT *t;                     /* to access the transaction buffer */

  assert(base && file);         /* check the function arguments */
  base->tract->size = 0;        /* initialize the item counter */
  d   = _read(base, file);      /* read the first item and */
  buf = ts_buf(base->tscan);    /* get the read buffer */
  if ((d      == TS_EOF)        /* if at the end of the file */
  &&  (buf[0] == '\0'))         /* and no item has been read, */
    return 1;                   /* return 'end of file' */
  while ((d      == TS_FLD)     /* read the other items */
  &&     (buf[0] != '\0'))      /* of the transaction */
    d = _read(base, file);      /* up to the end of the record */
  if (d == TS_ERR) return d;    /* check for a read error */
  t = base->tract;              /* get the transaction buffer */
  if ((buf[0] == '\0') && (d == TS_FLD) && (t->size > 0))
    return E_ITEMEXP;           /* check for an empty field */
  int_qsort(t->items, t->size); /* prepare the read transaction */
  t->size = int_unique(t->items, t->size);
  t->items[t->size] = -1;       /* store a sentinel after the items */
  base->wgt += t->wgt = 1;      /* sum the transaction weight */
  x = t->size *t->wgt;          /* compute extended frequency weight */
  for (i = t->size; --i >= 0;){ /* traverse the items */
    item = nim_byid(base->nimap, t->items[i]);
    item->frq += t->wgt;        /* sum the transaction weights */
    item->xfq += x;             /* and the transaction sizes */
  }                             /* as importance indicators */
  return 0;                     /* return 'ok' */
}  /* ib_read() */

/*--------------------------------------------------------------------*/

void ib_penfrq (ITEMBASE *base)
{                               /* --- include insertion penalties */
  int  i;                       /* loop variable */
  ITEM *item;                   /* to traverse the items */

  assert(base);                 /* check the function argument */
  for (i = nim_cnt(base->nimap); --i >= 0; ) {
    item = (ITEM*)nim_byid(base->nimap, i);
    item->frq += (int)ceil(item->pen *base->wgt);
  }                             /* recompute item frequencies */
}  /* ib_penfrq() */

/*--------------------------------------------------------------------*/

int ib_recode (ITEMBASE *base, int minfrq, int dir, int *map)
{                               /* --- recode items w.r.t. frequency */
  int   i, n, x;                /* loop variables, item buffer */
  ITEM  *item;                  /* to traverse the items */
  TRACT *t;                     /* to access the transaction buffer */
  int   *s, *d;                 /* to traverse the items */
  CMPFN *cmp;                   /* comparison function */

  assert(base);                 /* check the function arguments */
  if      (dir >  1) cmp = _asccmpx;  /* get the appropriate */
  else if (dir >  0) cmp = _asccmp;   /* comparison function */
  else if (dir >= 0) cmp = _nocmp;    /* (ascending/descending) */
  else if (dir > -2) cmp = _descmp;   /* and sort the items */
  else               cmp = _descmpx;  /* w.r.t. their frequency */
  nim_sort(base->nimap, cmp, &minfrq, map, 1);
  for (i = n = nim_cnt(base->nimap); --n >= 0; ) {
    item = (ITEM*)nim_byid(base->nimap, n);
    if (item->frq < minfrq)     /* determine frequent items and */
      item->app = APP_NONE;     /* set all others to 'ignore' */
    else if (item->app != APP_NONE)
      break;                    /* in addition, skip all items */
  }                             /* that have been set to 'ignore' */
  nim_trunc(base->nimap, ++n);  /* remove all items to ignore */
  if (!map) return n;           /* if no map is provided, abort */
  while (--i >= 0)              /* mark all removed items */
    if (map[i] >= n) map[i] = -1;
  t = base->tract;              /* traverse the buffered transaction */
  for (s = d = t->items; *s >= 0; s++) {
    x = map[*s];                /* recode all items and */
    if (x >= 0) *d++ = x;       /* remove all items to ignore */
  }                             /* from the buffered transaction */
  t->size = (int)(d -t->items); /* compute the new number of items */
  t->items[t->size] = -1;       /* store a sentinel after the items */
  int_qsort(t->items, t->size); /* and resort the remaining items */
  return n;                     /* return number of frequent items */
}  /* ib_recode() */

/*--------------------------------------------------------------------*/

void ib_trunc (ITEMBASE *base, int cnt)
{                               /* --- truncate an item base */
  TRACT *t;                     /* to access the transaction buffer */
  int   *s, *d;                 /* to traverse the items */

  assert(base && (cnt >= 0));   /* check the function arguments */
  nim_trunc(base->nimap, cnt);  /* truncate the item base */
  t = base->tract;              /* traverse the buffered transaction */
  for (s = d = t->items; *s >= 0; s++)
    if (*s < cnt) *d++ = *s;    /* remove all deleted items */
  t->size = (int)(d -t->items); /* compute the new number of items */
  t->items[t->size] = -1;       /* store a sentinel after the items */
}  /* ib_trunc() */

/*----------------------------------------------------------------------
  Transaction Functions
----------------------------------------------------------------------*/

TRACT* t_create (const int *items, int n, int wgt)
{                               /* --- create a transaction */
  TRACT *t;                     /* created transaction */

  assert(items || (n <= 0));    /* check the function arguments */
  t = (TRACT*)malloc(sizeof(TRACT) +n *sizeof(int));
  if (!t) return NULL;          /* allocate a new transaction */
  t->wgt      = wgt;            /* copy the transaction weight */
  t->items[n] = -1;             /* store a sentinel after the items */
  for (t->size = n; --n >= 0; ) /* copy all items into */
    t->items[n] = items[n];     /* the created transaction */
  return t;                     /* return the created transaction */
}  /* t_create() */

/*--------------------------------------------------------------------*/

TRACT* t_clone (const TRACT *t)
{ return t_create(t->items, t->size, t->wgt); }

/*--------------------------------------------------------------------*/

int t_cmp (const void *p1, const void *p2, void *data)
{                               /* --- compare transactions */
  const    int *i1, *i2;        /* to traverse the items */
  register int d;               /* difference of item identifiers */

  assert(p1 && p2);             /* check the function arguments */
  i1 = ((const TRACT*)p1)->items;
  i2 = ((const TRACT*)p2)->items;
  for ( ; 1; i1++, i2++) {      /* traverse the item arrays */
    d = *i1 -*i2;               /* compare corresponding items */
    if (d  != 0) return d;      /* and if one is greater, abort */
    if (*i1 < 0) return 0;      /* otherwise check for the sentinel */
  }                             /* and abort if it is reached */
}  /* t_cmp() */

/*--------------------------------------------------------------------*/

int t_cmpx (const TRACT *t, const int *items, int n)
{                               /* --- compare transactions */
  const    int *i;              /* to traverse the items */
  register int k, d;            /* loop variable, difference */

  assert(t && items);           /* check the function arguments */
  i = t->items;                 /* traverse the item array */
  for (k = (n < t->size) ? n : t->size; --k >= 0; i++, items++) {
    d = *i -*items;             /* compare corresponding items */
    if (d != 0) return d;       /* and abort the comparison */
  }                             /* if one of them is greater */
  return t->size -n;            /* return the size difference */
}  /* t_cmpx() */

/*----------------------------------------------------------------------
  Transaction Bag/Multiset Functions
----------------------------------------------------------------------*/

TABAG* tb_create (ITEMBASE *base)
{                               /* --- create a transaction bag */
  TABAG *bag;                   /* created transaction bag */

  assert(base);                 /* check the function argument */
  bag = malloc(sizeof(TABAG));  /* create a transaction bag/multiset */
  if (!bag) return NULL;        /* (just the base structure) */
  bag->base   = base;           /* store the underlying item base */
  bag->max    = bag->wgt = bag->size = bag->cnt = 0;
  bag->tracts = NULL;           /* initialize the other fields */
  return bag;                   /* return the created t.a. bag */
}  /* tb_create() */

/*--------------------------------------------------------------------*/

void tb_delete (TABAG *bag, int delis)
{                               /* --- delete a transaction bag */
  assert(bag);                  /* check the function argument */
  if (bag->tracts) {            /* if there are loaded transactions */
    while (--bag->cnt >= 0)     /* traverse the transaction array */
      free(bag->tracts[bag->cnt]);
    free(bag->tracts);          /* delete all transactions */
  }                             /* and the transaction array */
  if (bag->base && delis) ib_delete(bag->base);
  free(bag);                    /* delete the item base and */
}  /* tb_delete() */            /* the transaction bag body */

/*--------------------------------------------------------------------*/

int tb_add (TABAG *bag, TRACT *t)
{                               /* --- add a transaction */
  TRACT **p;                    /* new transaction array */
  int   n;                      /* new transaction array size */

  assert(bag);                  /* check the function arguments */
  n = bag->size;                /* get the transaction array size */
  if (bag->cnt >= n) {          /* if the transaction array is full */
    n += (n > BLKSIZE) ? (n >> 1) : BLKSIZE;
    p  = (TRACT**)realloc(bag->tracts, n *sizeof(TRACT*));
    if (!p) return -1;          /* enlarge the transaction array */
    bag->tracts = p; bag->size = n;
  }                             /* set the new array and its size */
  if (!t) { t = t_clone(ib_tract(bag->base));
            if (!t) return -1; }/* get trans. from item base if nec. */
  bag->tracts[bag->cnt++] = t;  /* store the transaction and */
  bag->wgt += t->wgt;           /* sum the transaction weight */
  if (t->size > bag->max) bag->max = t->size;
  return 0;                     /* update maximal transaction size */
}  /* tb_add() */               /* and return 'ok' */

/*--------------------------------------------------------------------*/

void tb_recode (TABAG *bag, int *map)
{                               /* --- recode items in transactions */
  int   k, n, x;                /* loop variable, item buffer */
  TRACT *t;                     /* to traverse the transactions */
  int   *s, *d;                 /* to traverse the items */

  assert(bag && map);           /* check the function arguments */
  bag->max = 0;                 /* clear maximal transaction size */
  for (n = bag->cnt; --n >= 0; ) {
    t = bag->tracts[n];         /* traverse the transactions */
    for (s = d = t->items; *s >= 0; s++) {
      x = map[*s];              /* traverse and recode the items */
      if (x >= 0) *d++ = x;     /* remove all items that are */
    }                           /* not mapped (mapped to id < 0) */
    k = (int)(d -t->items);     /* compute the new number of items */
    t->items[t->size = k] = -1; /* store a sentinel after the items */
    if (k > bag->max) bag->max = k;
  }                             /* update the maximal trans. size */
}  /* tb_recode() */

/*--------------------------------------------------------------------*/

void tb_filter (TABAG *bag, int min, const int *marks)
{                               /* --- filter (items in) transactions */
  int   k, n;                   /* loop variables */
  TRACT *t;                     /* to traverse the transactions */
  int   *s, *d;                 /* to traverse the items */

  assert(bag);                  /* check the function arguments */
  bag->max = 0;                 /* clear maximal transaction size */
  for (n = bag->cnt; --n >= 0; ) {
    t = bag->tracts[n];         /* traverse the transactions */
    if (marks) {                /* if item markers are given */
      for (s = d = t->items; *s >= 0; s++)
	if (marks[*s]) *d++ = *s;   /* remove unmarked items */
      t->size = k = (int)(d -t->items);
      t->items[k] = -1;         /* store the new number of items */
    }                           /* and a sentinel after the items */
    if (t->size < min) {        /* delete all items from */
      t->size     =  0;         /* those transactions that */
      t->items[0] = -1;         /* do not have the minimum size */
    }                           /* (cannot contribute to support) */
    if (t->size > bag->max)     /* update the maximal trans. size */
      bag->max = t->size;       /* (may differ from old size, because */
  }                             /* items may have been removed */
}  /* tb_filter() */

/*--------------------------------------------------------------------*/

void tb_itsort (TABAG *bag, int dir, int heap)
{                               /* --- sort items in transactions */
  int   i;                      /* loop variable */
  TRACT *t;                     /* to traverse the transactions */
  void  (*sortfn)(int*, int);   /* transaction sort function */

  assert(bag);                  /* check the function arguments */
  sortfn = (heap) ? int_heapsort : int_qsort;
  for (i = bag->cnt; --i >= 0; ) {
    t = bag->tracts[i];         /* traverse the transactions */
    sortfn(t->items, t->size);  /* and sort the items in them */
    if (dir < 0) int_reverse(t->items, t->size);
  }                             /* reverse the order if requested */
}  /* tb_itsort() */

/*--------------------------------------------------------------------*/

void tb_sort (TABAG *bag, int dir, int heap)
{                               /* --- sort a transaction bag */
  assert(bag);                  /* check the function arguments */
  if (heap) ptr_heapsort(bag->tracts, bag->cnt, t_cmp, NULL);
  else      ptr_qsort   (bag->tracts, bag->cnt, t_cmp, NULL);
  if (dir < 0) ptr_reverse(bag->tracts, bag->cnt);
}  /* tb_sort() */              /* sort the transactions */

/*--------------------------------------------------------------------*/

int tb_reduce (TABAG *bag)
{                               /* --- reduce a transaction bag */
  int   i;                      /* loop variable */
  TRACT **s, **d;               /* to traverse the transactions */

  assert(bag);                  /* check the function argument */
  if (bag->cnt <= 1) return 1;  /* deal only with two or more trans. */
  s = d = bag->tracts;          /* traverse the sorted transactions */
  for (i = bag->cnt; --i > 0; ) {
    if (t_cmp(*++s, *d, NULL) == 0) {
      (*d)->wgt += (*s)->wgt;   /* combine equal transactions */
      free(*s); }               /* by summing their weights */
    else {                      /* if transactions are not equal, */
      if ((*d)->wgt > 0) d++;   /* check weight of old transaction */
      else free(*d);            /* and delete it if it has no weight */
      *d = *s;                  /* copy the new transaction */
    }                           /* to close a possible gap */
  }                             /* (collect unique transactions) */
  if ((*d)->wgt > 0) d++;       /* check weight of last transaction */
  else free(*d);                /* and delete it if it has no weight */
  return bag->cnt = (int)(d -bag->tracts);
}  /* tb_reduce() */            /* return new number of transactions */

/*--------------------------------------------------------------------*/

int tb_occur (TABAG *bag, const int *items, int n)
{                               /* --- count transaction occurrences */
  int l, r, m, k;               /* index and loop variables */

  assert(bag && items);         /* check the function arguments */
  k = bag->cnt;                 /* get the number of transactions */
  for (r = m = 0; r < k; ) {    /* find right boundary */
    m = (r+k) >> 1;             /* by a binary search */
    if (t_cmpx(bag->tracts[m], items, n) > 0) k = m;
    else                                      r = m+1;
  }
  for (l = m = 0; l < k; ) {    /* find left boundary */
    m = (l+k) >> 1;             /* by a binary search */
    if (t_cmpx(bag->tracts[m], items, n) < 0) l = m+1;
    else                                      k = m;
  }
  for (k = 0; l < r; l++)       /* traverse the found section and */
    k += bag->tracts[l]->wgt;   /* sum the transaction weights */
  return k;                     /* return the number of occurrences */
}  /* tb_occur() */

/*--------------------------------------------------------------------*/
#ifndef NDEBUG

void tb_show (TABAG *bag, int wgt)
{                               /* --- show a transaction bag */
  int   i, k;                   /* loop variables */
  TRACT *t;                     /* to traverse the transactions */

  assert(bag);                  /* check the function argument */
  for (i = 0; i < bag->cnt; i++) {
    t = bag->tracts[i];         /* traverse the transactions */
    for (k = 0; k < t->size; k++) {    /* traverse the items */
      if (k > 0) fputc(bag->base->chars[1], stdout);
      printf(ib_name(bag->base, t->items[k]));
    }                           /* print the next item */
    if (wgt)                    /* print the transaction weight */
      printf("%c[%d]", bag->base->chars[1], t->wgt);
    fputc(bag->base->chars[2], stdout);
  }                             /* terminate the transaction */
  printf("%d/%d transaction(s)\n", bag->cnt, bag->wgt);
}  /* tb_show() */              /* finally print the number of t.a. */

#endif
/*----------------------------------------------------------------------
  Transaction Tree Functions
----------------------------------------------------------------------*/

void _delete (TTNODE *root)
{                               /* --- delete a transaction (sub)tree */
  int    i;                     /* loop variable */
  TTNODE **cnds;                /* array of child nodes */

  assert(root);                 /* check the function argument */
  cnds = (TTNODE**)(root->items +CHOFF(root->size));
  for (i = root->size; --i >= 0; )
    _delete(cnds[i]);           /* recursively delete the subtrees */
  free(root);                   /* and the tree node itself */
}  /* _delete() */

/*--------------------------------------------------------------------*/

TTNODE* _create (TRACT **tracts, int cnt, int index)
{                               /* --- recursive part of tt_create() */
  int    i, k, t, w;            /* loop variables, buffers */
  int    item, n;               /* item and item counter */
  TTNODE *root;                 /* root of created transaction tree */
  TTNODE **cnds;                /* array of child nodes */
  int    *s, *d;                /* to traverse the items */

  assert(tracts                 /* check the function arguments */
     && (cnt >= 0) && (index >= 0));
  for (w = 0, k = cnt; --k >= 0; )
    w += tracts[k]->wgt;        /* determine the total trans. weight */
  if (w <= 0) cnt = 0;          /* check for weightless transactions */

  if (cnt <= 1) {               /* if only one transaction left */
    n    = (cnt > 0) ? (*tracts)->size -index : 0;
    root = (TTNODE*)malloc(sizeof(TTNODE) +n *sizeof(int));
    if (!root) return NULL;     /* create a transaction tree node */
    root->wgt  =  w;            /* and initialize the fields */
    root->max  =  n;
    root->size = -n;
    d = root->items; d[n] = -1; /* place a sentinel at the end */
    s = (*tracts)->items +index;
    while (--n >= 0) d[n] = s[n];
    return root;                /* copy the remaining items and */
  }                             /* return the created leaf node */

  for (k = cnt; (--k >= 0) && ((*tracts)->size <= index); )
    tracts++;                   /* skip t.a. that are too short */
  for (n = 0, item = -1, i = ++k; --i >= 0; ) {
    t = tracts[i]->items[index];/* traverse the transactions */
    if (t != item) { item = t; n++; }
  }                             /* count the different items */
  i    = CHOFF(n);              /* get offset to the child pointers */
  root = (TTNODE*)malloc(sizeof(TTNODE) +(i-1) *sizeof(int)
                                        + n    *sizeof(TTNODE*));
  if (!root) return NULL;       /* create a transaction tree node */
  root->wgt  = w;               /* and initialize its fields */
  root->max  = 0;
  root->size = n;               /* if transactions are captured, */
  if (n <= 0) return root;      /* return the created tree */

  cnds = (TTNODE**)(root->items +i);
  for (--k; --n >= 0; k = i) {  /* traverse the different items */
    root->items[n] = item = tracts[k]->items[index];
    for (i = k; --i >= 0; )     /* find trans. with the current item */
      if (tracts[i]->items[index] != item) break;
    cnds[n] = _create(tracts +i+1, k-i, index+1);
    if (!cnds[n]) break;        /* recursively create a subtree */
    t = cnds[n]->max +1;        /* adapt the maximal remaining size */
    if (t > root->max) root->max = t;
  }                             /* (create all subtrees) */
  if (n < 0) return root;       /* if successful, return created tree */

  for (i = root->size; --i > n; )
    _delete(cnds[i]);           /* on error delete created subtrees */
  free(root);                   /* and the transaction tree node */
  return NULL;                  /* return 'failure' */
}  /* _create() */

/*--------------------------------------------------------------------*/

TATREE* tt_create (TABAG *bag)
{                               /* --- create a transactions tree */
  TATREE *tree;                 /* created transaction tree */

  assert(bag);                  /* check the function argument */
  tree = (TATREE*)malloc(sizeof(TATREE));
  if (!tree) return NULL;       /* create the transaction tree body */
  tree->base = bag->base;       /* note the underlying item set */
  tree->root = _create(bag->tracts, bag->cnt, 0);
  if (!tree->root) { free(tree); return NULL; }
  return tree;                  /* recursively build the tree */
}  /* tt_create() */

/*--------------------------------------------------------------------*/

void tt_delete (TATREE *tree, int delis)
{                               /* --- delete a transaction tree */
  assert(tree);                 /* check the function argument */
  _delete(tree->root);          /* delete the nodes of the tree */
  if (tree->base && delis) ib_delete(tree->base);
  free(tree);                   /* delete the item base and */
}  /* tt_delete() */            /* the transaction tree body */

/*--------------------------------------------------------------------*/

static int _nodecnt (TTNODE *root)
{                               /* --- count the nodes */
  int    i, n;                  /* loop variable, number of nodes */
  TTNODE **cnds;                /* array of child nodes */

  assert(root);                 /* check the function argument */
  if (root->size <= 0)          /* if this is a leaf node, */
    return 1;                   /* there is only one node */
  cnds = (TTNODE**)(root->items +CHOFF(root->size));
  for (n = 0, i = root->size; --i >= 0; )
    n += _nodecnt(cnds[i]);     /* recursively count the nodes */
  return n+1;                   /* return number of nodes in tree */
}  /* _nodecnt() */

/*--------------------------------------------------------------------*/

int tt_nodecnt (TATREE *tree)
{ return _nodecnt(tree->root); }

/*--------------------------------------------------------------------*/

static int _extcnt (TTNODE *root)
{                               /* --- extended node counting */
  int    i, n;                  /* loop variable, number of nodes */
  TTNODE **cnds;                /* array of child nodes */

  assert(root);                 /* check the function argument */
  if (root->size <= 0)          /* if this is a leaf node, */
    return -root->size;         /* return the number of items */
  cnds = (TTNODE**)(root->items +CHOFF(root->size));
  for (n = 0, i = root->size; --i >= 0; )
    n += _extcnt(cnds[i]);      /* recursively count the nodes */
  return n+1;                   /* return number of nodes in tree */
}  /* _extcnt() */

/*--------------------------------------------------------------------*/

int tt_extcnt (TATREE *tree)
{ return _extcnt(tree->root); }

/*--------------------------------------------------------------------*/
#ifdef ARCH64

TTNODE* ttn_child (TTNODE *node, int index)
{                               /* --- go to a child node */
  assert(node                   /* check the function arguments */
     && (index >= 0) && (index < node->size));
  return ((TTNODE**)(node->items +CHOFF(node->size)))[index];
}  /* ttn_child */              /* return the child node/subtree */

#endif
/*--------------------------------------------------------------------*/
#ifndef NDEBUG

void _show (TTNODE *node, ITEMBASE *base, int ind)
{                               /* --- rekursive part of tt_show() */
  int    i, k;                  /* loop variables */
  TTNODE **cnds;                /* array of child nodes */

  assert(node && (ind >= 0));   /* check the function arguments */
  if (node->size <= 0) {        /* if this is a leaf node */
    for (i = 0; i < node->max; i++)
      printf("%s ", ib_name(base, node->items[i]));
    printf("[%d]\n", node->wgt);
    return;                     /* print the items in the */
  }                             /* (rest of) the transaction */
  cnds = (TTNODE**)(node->items +node->size);
  for (i = 0; i < node->size; i++) {
    if (i > 0) for (k = ind; --k >= 0; ) printf("  ");
    printf("%s ", ib_name(base, node->items[i]));
    _show(cnds[i], base, ind+1);/* traverse the items, print them, */
  }                             /* and show the children recursively */
}  /* _show() */

/*--------------------------------------------------------------------*/

void tt_show (TATREE *tree)
{                               /* --- show a transaction tree */
  assert(tree);                 /* check the function argument */
  _show(tree->root, tree->base, 0);
}  /* tt_show() */              /* call the recursive function */

#endif
