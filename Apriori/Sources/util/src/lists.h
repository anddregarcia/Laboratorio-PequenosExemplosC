
#ifndef __LISTS__
#define __LISTS__
#include "fntypes.h"

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _dlle {          /* --- doubly linked list element --- */
  struct _dlle *succ;           /* pointer to successor */
  struct _dlle *pred;           /* pointer to predecessor */
} DLLE;                         /* (doubly linked list element) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern void  l_delete  (void *list, OBJFN delfn);
extern void* l_reverse (void *list);
extern void* l_append  (void *dst, void *src);
extern void* l_merge   (void *in1, void *in2, CMPFN cmpfn, void *data);
extern void* l_sort    (void *list, CMPFN cmpfn, void *data);

#endif
