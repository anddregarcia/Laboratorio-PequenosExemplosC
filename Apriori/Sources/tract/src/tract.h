
#ifndef __TRACT__
#define __TRACT__
#ifndef NIMAPFN
#define NIMAPFN
#endif
#include "arrays.h"
#include "symtab.h"
#include "tabscan.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- item appearance flags --- */
#define APP_NONE    0x00        /* item should be ignored */
#define APP_BODY    0x01        /* item may appear in rule body */
#define APP_HEAD    0x02        /* item may appear in rule head */
#define APP_BOTH    (APP_HEAD|APP_BODY)

/* --- error codes --- */
#define E_NONE         0        /* no error */
#define E_NOMEM      (-1)       /* not enough memory */
#define E_FOPEN      (-2)       /* cannot open file */
#define E_FREAD      (-3)       /* read error on file */
#define E_FWRITE     (-4)       /* write error on file */

#define E_ITEMEXP   (-16)       /* item expected */
#define E_DUPITEM   (-17)       /* duplicate item */
#define E_FLDCNT    (-18)       /* too many fields */
#define E_APPEXP    (-19)       /* appearance indicator expected */
#define E_UNKAPP    (-20)       /* unknown appearance indicator */
#define E_PENEXP    (-21)       /* insertion penalty expected */
#define E_PENALTY   (-22)       /* invalid insertion penalty */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- an item --- */
  int      id;                  /* item identifier */
  int      frq;                 /* frequency in transactions */
  int      xfq;                 /* extended frequency (trans. sizes) */
  int      app;                 /* appearance indicator */
  double   pen;                 /* insertion penalty */
} ITEM;                         /* (item) */

typedef struct {                /* --- a transaction --- */
  int      size;                /* size   (number of items) */
  int      wgt;                 /* weight (number of occurrences) */
  int      items[1];            /* items in the transaction */
} TRACT;                        /* (transaction) */

typedef struct {                /* --- an item base --- */
  NIMAP    *nimap;              /* name/identifier map */
  TABSCAN  *tscan;              /* table scanner */
  char     chars[4];            /* special characters */
  int      wgt;                 /* total weight of transactions */
  int      app;                 /* default appearance indicator */
  double   pen;                 /* default insertion penalty */
  int      size;                /* size of the transaction buffer */
  TRACT    *tract;              /* buffer for a transaction */
} ITEMBASE;                     /* (item base) */

typedef struct {                /* --- a transaction bag/multiset --- */
  ITEMBASE *base;               /* underlying item base */
  int      max;                 /* number of items in largest trans. */
  int      wgt;                 /* total weight of transactions */
  int      size;                /* size of the transaction array */
  int      cnt;                 /* number of transactions */
  TRACT    **tracts;            /* transaction array */
} TABAG;                        /* (transaction bag/multiset) */

typedef struct {                /* --- a transaction tree node --- */
  int      wgt;                 /* weight (number of transactions) */
  int      max;                 /* number of items in largest trans. */
  int      size;                /* node size (number of children) */
  int      items[1];            /* next items in rep. transactions */
} TTNODE;                       /* (transaction tree node) */

typedef struct {                /* --- a transaction tree --- */
  ITEMBASE *base;               /* underlying item base */
  TTNODE   *root;               /* root of the transaction tree */
} TATREE;                       /* (transaction tree) */

/*----------------------------------------------------------------------
  Item Base Functions
----------------------------------------------------------------------*/
extern ITEMBASE*   ib_create  (int size);
extern void        ib_delete  (ITEMBASE *base);
extern TABSCAN*    ib_tabscan (ITEMBASE *base);
extern void        ib_chars   (ITEMBASE *base, const char *blanks,
                                               const char *fldseps,
                                               const char *recseps,
                                               const char *comment);

extern int         ib_cnt     (ITEMBASE *base);
extern int         ib_add     (ITEMBASE *base, const char *name);
extern int         ib_item    (ITEMBASE *base, const char *name);
extern const char* ib_name    (ITEMBASE *base, int item);

extern int         ib_getwgt  (ITEMBASE *base);
extern int         ib_setwgt  (ITEMBASE *base, int cnt);
extern int         ib_incwgt  (ITEMBASE *base, int cnt);

extern int         ib_getfrq  (ITEMBASE *base, int item);
extern int         ib_setfrq  (ITEMBASE *base, int item, int frq);
extern int         ib_incfrq  (ITEMBASE *base, int item, int frq);
extern int         ib_getxfq  (ITEMBASE *base, int item);
extern int         ib_setxfq  (ITEMBASE *base, int item, int frq);
extern int         ib_incxfq  (ITEMBASE *base, int item, int frq);
extern int         ib_getapp  (ITEMBASE *base, int item);
extern int         ib_setapp  (ITEMBASE *base, int item, int app);
extern double      ib_getpen  (ITEMBASE *base, int item);
extern double      ib_setpen  (ITEMBASE *base, int item, double pen);

extern int         ib_readapp (ITEMBASE *base, FILE *file);
extern int         ib_readpen (ITEMBASE *base, FILE *file);
extern int         ib_read    (ITEMBASE *base, FILE *file);

extern void        ib_penfrq  (ITEMBASE *base);
extern int         ib_recode  (ITEMBASE *base, int minfrq,
                               int dir, int *map);
extern void        ib_trunc   (ITEMBASE *base, int cnt);

extern TRACT*      ib_tract   (ITEMBASE *base);

/*----------------------------------------------------------------------
  Transaction Functions
----------------------------------------------------------------------*/
extern TRACT*      t_create   (const int *items, int n, int wgt);
extern void        t_delete   (TRACT *t);
extern TRACT*      t_clone    (const TRACT *t);

extern const int*  t_items    (const TRACT *t);
extern int         t_size     (const TRACT *t);
extern int         t_wgt      (const TRACT *t);

extern void        t_sort     (TRACT *t);
extern void        t_reverse  (TRACT *t);

extern int         t_cmp      (const void *p1,
                               const void *p2, void *data);
extern int         t_cmpx     (const TRACT *t, const int *items, int n);

/*----------------------------------------------------------------------
  Transaction Bag/Multiset Functions
----------------------------------------------------------------------*/
extern TABAG*      tb_create  (ITEMBASE *base);
extern void        tb_delete  (TABAG *bag, int delis);
extern ITEMBASE*   tb_base    (TABAG *bag);

extern int         tb_cnt     (TABAG *bag);
extern int         tb_wgt     (TABAG *bag);
extern int         tb_max     (TABAG *bag);

extern int         tb_add     (TABAG *bag, TRACT *t);
extern int         tb_addx    (TABAG *bag,
                               const int *items, int n, int wgt);
extern TRACT*      tb_tract   (TABAG *bag, int index);

extern void        tb_recode  (TABAG *bag, int *map);
extern void        tb_filter  (TABAG *bag, int min, const int *marks);
extern void        tb_itsort  (TABAG *bag, int dir, int heap);
extern void        tb_sort    (TABAG *bag, int dir, int heap);
extern int         tb_reduce  (TABAG *bag);
extern void        tb_shuffle (TABAG *bag, double randfn(void));
extern int         tb_occur   (TABAG *bag, const int *items, int n);

#ifndef NDEBUG
extern void        tb_show    (TABAG *bag, int wgt);
#endif

/*----------------------------------------------------------------------
  Transaction Tree Functions
----------------------------------------------------------------------*/
extern TATREE*     tt_create  (TABAG *bag);
extern void        tt_delete  (TATREE *tree, int delis);
extern ITEMBASE*   tt_base    (TATREE *tree);

extern TTNODE*     tt_root    (TATREE *tree);
extern int         tt_nodecnt (TATREE *tree);
extern int         tt_extcnt  (TATREE *tree);
extern int         tt_wgt     (TATREE *tree);
extern int         tt_max     (TATREE *tree);

#ifndef NDEBUG
extern void        tt_show    (TATREE *tree);
#endif

/*----------------------------------------------------------------------
  Transaction Tree Node Functions
----------------------------------------------------------------------*/
extern int         ttn_wgt    (TATREE *tree);
extern int         ttn_max    (TATREE *tree);
extern int         ttn_size   (TATREE *tree);
extern int*        ttn_items  (TATREE *tree);
extern int         ttn_item   (TATREE *tree, int index);
extern TTNODE*     ttn_child  (TATREE *tree, int index);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define ib_tabscan(s)     ((s)->tscan)

#define ib_cnt(s)         nim_cnt((s)->nimap)
#define ib_item(s,n)      nim_getid((s)->nimap, n)
#define ib_name(s,i)      nim_name(nim_byid((s)->nimap, i))

#define ib_getwgt(s)      ((s)->wgt)
#define ib_setwgt(s,n)    ((s)->wgt  = (n))
#define ib_incwgt(s,n)    ((s)->wgt += (n))

#define ib_getfrq(s,i)    (((ITEM*)nim_byid((s)->nimap, i))->frq)
#define ib_setfrq(s,i,n)  (((ITEM*)nim_byid((s)->nimap, i))->frq  = (n))
#define ib_incfrq(s,i,n)  (((ITEM*)nim_byid((s)->nimap, i))->frq += (n))
#define ib_getxfq(s,i)    (((ITEM*)nim_byid((s)->nimap, i))->xfq)
#define ib_setxfq(s,i,n)  (((ITEM*)nim_byid((s)->nimap, i))->xfq  = (n))
#define ib_incxfq(s,i,n)  (((ITEM*)nim_byid((s)->nimap, i))->xfq += (n))
#define ib_getapp(s,i)    (((ITEM*)nim_byid((s)->nimap, i))->app)
#define ib_setapp(s,i,a)  (((ITEM*)nim_byid((s)->nimap, i))->app  = (a))
#define ib_getpen(s,i)    (((ITEM*)nim_byid((s)->nimap, i))->pen)
#define ib_setpen(s,i,p)  (((ITEM*)nim_byid((s)->nimap, i))->pen  = (p))

#define ib_tract(s)       ((s)->tract)

/*--------------------------------------------------------------------*/
#define t_delete(t)       free(t)
#define t_sort(t)         int_qsort((t)->items, (t)->size)
#define t_reverse(t)      int_reverse((t)->items, (t)->size)
#define t_items(t)        ((t)->items)
#define t_size(t)         ((t)->size)
#define t_wgt(t)          ((t)->wgt)

/*--------------------------------------------------------------------*/
#define tb_base(b)        ((b)->base)

#define tb_cnt(b)         ((b)->cnt)
#define tb_wgt(b)         ((b)->wgt)
#define tb_max(b)         ((b)->max)
#define tb_addx(b,i,n,w)  tb_add(b, t_create(i,n,w))

#define tb_tract(b,i)     ((b)->tracts[i])

#define tb_shuffle(b,f)   ptr_shuffle((b)->tracts, (b)->cnt, f)

/*--------------------------------------------------------------------*/
#define tt_base(t)        ((t)->base)
#define tt_root(t)        ((t)->root)
#define tt_wgt(t)         ((t)->root->wgt)
#define tt_max(t)         ((t)->root->max)

/*--------------------------------------------------------------------*/
#define ttn_wgt(n)        ((n)->wgt)
#define ttn_max(n)        ((n)->max)
#define ttn_size(n)       ((n)->size)
#define ttn_item(n,i)     ((n)->items[i])
#define ttn_items(n)      ((n)->items)
#ifndef ARCH64
#define ttn_child(n,i)    (((TTNODE**)((n)->items +(n)->size))[i])
#endif

#endif
