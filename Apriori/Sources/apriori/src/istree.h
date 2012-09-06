/*
Projeto de Organização e recuperação da Informação (ORI)
Escola de Engenharia de Piracicaba
Authors: Adilson Perecin(a.perecin(at)hotmail.com) & Cristiano Benato(benato(at)hst.com.br)

Tema: Organização e verificação de itemsets e ordenação de coincidencias através
do algoritmo de ordenação Apriori

Instrutor: Luiz Camolesi
 


*/

#ifndef __ISTREE__
#define __ISTREE__
#include <limits.h>
#include "report.h"


/* --- operation modes --- */
#define IST_PERFECT  INT_MIN    /* prune with perfect extensions */

/* --- additional evaluation measures --- */
#define IST_NONE       0        /* no measure */
#define IST_CONF       1        /* rule confidence */
#define IST_DIFF       2        /* abs. confidence diff. to prior */
#define IST_LIFT       3        /* lift value (confidence/prior) */
#define IST_LD21       4        /* abs. diff. of lift value to 1 */
#define IST_QUOT       5        /* difference of lift quotient to 1 */
#define IST_CHI2       6        /* normalized chi^2 measure */
#define IST_PVAL       7        /* p-value from chi^2 measure */
#define IST_INFO       8        /* information difference to prior */
#define IST_PGST       9        /* p-value from G statistic */
#define IST_LOGQ      10        /* binary log. of support quotient */

/* --- evaluation measure aggregation modes --- */
/*      IST_NONE       0 */     /* no aggregation (use first value) */
#define IST_MIN        1        /* minimum of measure values */
#define IST_MAX        2        /* maximum of measure values */
#define IST_AVG        3        /* average of measure values */

/* --- item set mark modes --- */
#define IST_CLEAR      0        /* clear markers */
#define IST_CLOSED     1        /* closed  item sets (mark complement)*/
#define IST_MAXIMAL    2        /* maximal item sets (mark complement)*/
#define IST_EVAL       4        /* filter also with evaluation */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _isnode {        /* --- item set node --- */
  struct _isnode *parent;       /* parent node */
  struct _isnode *succ;         /* successor node on same level */
  int            id;            /* item id used in parent node */
  int            size;          /* size   of counter array */
  int            offset;        /* offset of counter array */
  int            chcnt;         /* number of child nodes */
  int            cnts[1];       /* counter array (weights) */
} ISNODE;                       /* (item set node) */

typedef struct {                /* --- item set tree --- */
  ITEMBASE *base;               /* underlying item base */
  int      mode;                /* search mode (e.g. support def.) */
  int      wgt;                 /* total weight of transactions */
  int      height;              /* tree height (number of levels) */
  int      maxht;               /* max. height (size of level array) */
  ISNODE   **lvls;              /* first node of each level */
  int      rule;                /* minimal support of an assoc. rule */
  int      supp;                /* minimal support of an item set */
  int      smax;                /* maximal support of an item set */
  double   conf;                /* minimal confidence of a rule */
  int      eval;                /* additional evaluation measure */
  int      agg;                 /* aggregation mode of measure values */
  double   minval;              /* minimal evaluation measure value */
  ISNODE   *curr;               /* current node for traversal */
  int      size;                /* current size of an item set */
  int      minsz;               /* minimal size of an item set */
  int      maxsz;               /* maximal size of an item set */
  int      dir;                 /* direction for output acc. to size */
  ISNODE   *node;               /* item set node for extraction */
  int      index;               /* index in item set node */
  ISNODE   *head;               /* head item node for extraction */
  int      prune;               /* start level for evaluation pruning */
  int      item;                /* head item of previous rule */
  int      *buf;                /* buffer for paths (support check) */
  int      *path;               /* current path / (partial) item set */
  int      hdonly;              /* head only item in current set */
  int      *map;                /* to create identifier maps */
#ifdef BENCH                    /* if benchmark version */
  int      ndcnt;               /* number of item set tree nodes */
  int      ndprn;               /* number of pruned tree nodes */
  int      mapsz;               /* number of elements in id. maps */
  int      sccnt;               /* number of created support counters */
  int      scnec;               /* number of necessary supp. counters */
  int      scprn;               /* number of pruned support counters */
  int      cpcnt;               /* number of created child pointers */
  int      cpnec;               /* number of necessary child pointers */
  int      cpprn;               /* number of pruned child pointers */
#endif
} ISTREE;                       /* (item set tree) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern ISTREE* ist_create  (ITEMBASE *base, int mode,
                            int supp, int smax, double conf);
extern void    ist_delete  (ISTREE *ist);
extern int     ist_itemcnt (ISTREE *ist);

extern void    ist_count   (ISTREE *ist,
                            const int *items, int n, int wgt);
extern void    ist_countt  (ISTREE *ist, const TRACT  *tract);
extern void    ist_countb  (ISTREE *ist, const TABAG  *bag);
extern void    ist_countx  (ISTREE *ist, const TATREE *tree);

extern void    ist_prune   (ISTREE *ist);
extern int     ist_check   (ISTREE *ist, int *marks);
extern int     ist_addlvl  (ISTREE *ist);

extern int     ist_height  (ISTREE *ist);
extern int     ist_getwgt  (ISTREE *ist);
extern int     ist_setwgt  (ISTREE *ist, int wgt);
extern int     ist_incwgt  (ISTREE *ist, int wgt);

extern void    ist_up      (ISTREE *ist, int root);
extern int     ist_down    (ISTREE *ist, int item);
extern int     ist_next    (ISTREE *ist, int item);
extern int     ist_supp    (ISTREE *ist, int item);
extern int     ist_suppx   (ISTREE *ist, int *items, int cnt);

extern void    ist_mark    (ISTREE *ist, int mode);
extern void    ist_setsize (ISTREE *ist, int min, int max, int dir);
extern void    ist_seteval (ISTREE *ist, int eval, int agg, double min,
                            int prune);

extern void    ist_init    (ISTREE *ist);
extern int     ist_set     (ISTREE *ist, int *items, int *supp,
                            double *eval);
extern int     ist_rule    (ISTREE *ist, int *rule,  int *supp,
                            int *body, int *head, double *eval);

extern int     ist_report  (ISTREE *ist, ISREPORT *rep);
extern double  ist_eval    (ISTREE *ist);
extern double  ist_evalx   (ISREPORT *rep, void *data);

#ifndef NDEBUG
extern void    ist_show    (ISTREE *ist);
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define ist_itemcnt(t)     ((t)->levels[0]->size)
#define ist_height(t)      ((t)->height)
#define ist_getwgt(t)      ((t)->wgt & ~INT_MIN)
#define ist_setwgt(t,n)    ((t)->wgt = (n))
#define ist_incwgt(t,n)    ((t)->wgt = ((t)->wgt & ~INT_MIN) +(n))

#endif
