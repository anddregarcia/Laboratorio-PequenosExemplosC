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
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "scan.h"
#include "istree.h"
#ifdef STORAGE
#include "storage.h"
#endif


#define TT_SET      IST_CLEAR
#define TT_CLOSED   IST_CLOSED
#define TT_MAXIMAL  IST_MAXIMAL
#define TT_RULE     IST_EVAL

/* --- error codes --- */
#define E_OPTION     (-5)
#define E_OPTARG     (-6)
#define E_ARGCNT     (-7)
#define E_STDIN      (-8)
#define E_TARGET     (-9)
#define E_SIZE      (-10)
#define E_SUPP      (-11)
#define E_CONF      (-12)
#define E_MEASURE   (-13)
#define E_NOTRANS   (-14)
#define E_NOFREQ    (-15)
#define E_UNKNOWN   (-23)
#define PRGNAME     "\n\nApriori"
#define DESCRIPTION "\n*_________________________________________*\nDeveloped by Cristiano Benato & Adilson Perecin\n\n"
#define VERSION     "Computer Science  \n" \
                    "Escola de Eng. de Piracicaba (EEP)\n_________________________________________\n"

#ifndef QUIET
#define MSG         fprintf
#else
#define MSG(...)
#endif

#define SEC_SINCE(t)  ((clock()-(t)) /(double)CLOCKS_PER_SEC)
#define RECCNT(s)     (ts_reccnt(ib_tabscan(s)) \
                      - ((ts_delim(ib_tabscan(s)) == TS_REC) ? 1 : 0))
#define BUFFER(s)     ts_buf(ib_tabscan(s))

#ifndef QUIET
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error\n",
  /* E_NOMEM    -1 */  "not enough memory\n",
  /* E_FOPEN    -2 */  "cannot open file %s\n",
  /* E_FREAD    -3 */  "read error on file %s\n",
  /* E_FWRITE   -4 */  "write error on file %s\n",
  /* E_OPTION   -5 */  "unknown option -%c\n",
  /* E_OPTARG   -6 */  "missing option argument\n",
  /* E_ARGCNT   -7 */  "wrong number of arguments\n",
  /* E_STDIN    -8 */  "double assignment of standard input\n",
  /* E_TARGET   -9 */  "invalid target type '%c'\n",
  /* E_SIZE    -10 */  "invalid item set or rule size %d\n",
  /* E_SUPP    -11 */  "invalid minimum support %g%%\n",
  /* E_CONF    -12 */  "invalid minimum confidence %g%%\n",
  /* E_MEASURE -13 */  "invalid evaluation measure "
                         "or aggregation mode %c\n",
  /* E_NOTRANS -14 */  "no items or transactions to work on\n",
  /* E_NOFREQ  -15 */  "no frequent items found\n",
  /* E_ITEMEXP -16 */  "file %s, record %d: item expected\n",
  /* E_DUPITEM -17 */  "file %s, record %d: duplicate item %s\n",
  /* E_FLDCNT  -18 */  "file %s, record %d: too many fields\n",
  /* E_APPEXP  -19 */  "file %s, record %d: "
                         "appearance indicator expected\n",
  /* E_UNKAPP  -20 */  "file %s, record %d: "
                         "unknown appearance indicator %s\n",
  /*    -21 to -22 */  NULL, NULL,
  /* E_UNKNOWN -23 */  "unknown error\n"
};
#endif

#ifndef QUIET
static char     *prgname;
#endif
static ITEMBASE *ibase  = NULL;
static TABAG    *tabag  = NULL;
static TATREE   *tatree = NULL;
static ISTREE   *istree = NULL;
static ISREPORT *isrep  = NULL;
static int      *map    = NULL;
static FILE     *in     = NULL;
static FILE     *out    = NULL;

static void help (void)
{
  #ifndef QUIET
  fprintf(stderr, "\n");
  printf("additional evaluation measures (option -e#)\n");
  printf("association rules:\n");
  printf("  x   no measure\n");
  printf("  c   rule confidence\n");
  printf("  d   absolute confidence difference to prior\n");
  printf("  l   lift value (confidence divided by prior)\n");
  printf("  a   absolute difference of lift value to 1\n");
  printf("  r   difference of lift quotient to 1\n");
  printf("  n   normalized chi^2 measure\n");
  printf("  p   p-value computed from chi^2 measure\n");
  printf("  i   information difference to prior\n");
  printf("  g   p-value computed from G statistic\n");
  printf("frequent item sets:\n");
  printf("  x   no measure\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("All measures for association rules are also applicable,\n");
  printf("which are then aggregated over all possible rules.\n");
  printf("The aggregation mode can be set with the option -x#.\n");
  printf("\n");
  printf("evaluation measure aggregation modes (option -x#)\n");
  printf("  x   no aggregation (use first value)\n");
  printf("  m   minimum of individual measure values\n");
  printf("  n   maximum of individual measure values\n");
  printf("  a   average of individual measure values\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%a  absolute item set  support\n");
  printf("  %%s  relative item set  support as a fraction\n");
  printf("  %%S  relative item set  support as a percentage\n");
  printf("  %%b  absolute body set  support\n");
  printf("  %%x  relative body set  support as a fraction\n");
  printf("  %%X  relative body set  support as a percentage\n");
  printf("  %%h  absolute head item support\n");
  printf("  %%y  relative head item support as a fraction\n");
  printf("  %%Y  relative head item support as a percentage\n");
  printf("  %%c  rule confidence as a fraction\n");
  printf("  %%C  rule confidence as a percentage\n");
  printf("  %%l  lift value of a rule (confidence/prior)\n");
  printf("  %%L  lift value of a rule as a percentage\n");
  printf("  %%e  additional evaluation measure\n");
  printf("  %%E  additional evaluation measure as a percentage\n");
  printf("s,S,x,X,y,Y,c,C,l,L,e,E can be preceded by the number\n");
  printf("of decimal places to be printed (at most 32 places).\n");
  #endif
  exit(0);
}


static void error (int code, ...)
{                               
  #ifndef QUIET                
  va_list    args;           
  const char *msg;            

  assert(prgname);            
  if (code < E_UNKNOWN) code = E_UNKNOWN;
  if (code < 0) {              
    msg = errmsgs[-code];      
    if (!msg) msg = errmsgs[-E_UNKNOWN];
    fprintf(stderr, "\n%s: ", prgname);
    va_start(args, code);      
    vfprintf(stderr, msg, args);
    va_end(args);               
  }
  #endif
  #ifndef NDEBUG                
  if (map)    free(map);        
  if (isrep)  isr_delete(isrep, 0); 
  if (istree) ist_delete(istree);
  if (tatree) tt_delete(tatree, 0);
  if (tabag)  tb_delete(tabag, 0);
  if (ibase)  ib_delete(ibase);
  if (in  && (in  != stdin))  fclose(in);
  if (out && (out != stdout)) fclose(out);
  #endif
  #ifdef STORAGE       
  showmem("at end of program"); 
  #endif
  exit(code);                  
}  /* error() */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                              
  int     i, k = 0, n;          
  char    *s;                   
  char    **optarg = NULL;     
  char    *fn_in   = NULL;      
  char    *fn_out  = NULL;    
  char    *fn_app  = NULL;     
  char    *blanks  = NULL;   
  char    *fldseps = NULL;     
  char    *recseps = NULL;      
  char    *comment = NULL;      
  char    *isep    = " ";      
  char    *impl    = " <- ";   
  char    *dflt    = "  (%1S)"; 
  char    *format  = dflt;      
  int     target   = 's';     
  int     min      = 1;        
  int     max      = INT_MAX;   
  double  supp     = 0.1;       
  double  smax     = 1.0;     
  double  conf     = 0.8;       
  int     dir      = 0;        
  int     eval     = 0;       
  int     aggm     = 0;        
  double  minval   = 0.1;       
  int     prune    = 0;         
  double  filter   = 0.1;       
  int     sort     = 2;         
  int     tree     = 1;         
  int     heap     = 1;         
  int     post     = 0;        
  int     report   = 0;       
  int     mode     = APP_BODY|IST_PERFECT;  
  int     size;               
  int     wgt;                  
  int     frq, body, head;     
  int     *items;             
  clock_t t, tt, tc, x;         

  #ifndef QUIET               
  prgname = argv[0];           

  if (argc > 1) {          
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } 
  else {                      
    printf("usage: %s [options] infile outfile\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent item sets, c: closed item sets,\n"
           "          m: maximal item sets,  r: association rules)\n");
    printf("-m#      minimum number of items per set/rule     "
                    "(default: %d)\n", min);
    printf("-n#      maximum number of items per set/rule     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of a set/rule     "
                    "(default: %g%%)\n", supp *100);
    printf("-S#      maximum support of a set/rule     "
                    "(default: %g%%)\n", smax *100);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-c#      minimum confidence of a     rule         "
                    "(default: %g%%)\n", conf *100);
    printf("infile   file to read transactions from\n");
    printf("outfile  file to write item sets to\n");
    return 0;                
  }                            
  #endif  
  for (i = 1; i < argc; i++) {  
    s = argv[i];                
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  
      while (*s) {             
        switch (*s++) {        
          case '!': help();                         break;
          case 't': target = (*s) ? *s++ : 's';     break;
          case 'm': min    = (int)strtol(s, &s, 0); break;
          case 'n': max    = (int)strtol(s, &s, 0); break;
          case 's': supp   = 0.01*strtod(s, &s);    break;
          case 'S': smax   = 0.01*strtod(s, &s);    break;
          case 'c': conf   = 0.01*strtod(s, &s);    break;
          case 'o': mode  |= APP_BOTH;              break;
          case 'e': eval   = (*s) ? *s++ : 0;       break;
          case 'a': aggm   = (*s) ? *s++ : 0;       break;
          case 'd': minval = 0.01*strtod(s, &s);    break;
          case 'p': prune  = (int)strtol(s, &s, 0); break;
          case 'g': report = ISR_SCAN;              break;
          case 'k': optarg = &isep;                 break;
          case 'i': optarg = &impl;                 break;
          case 'v': optarg = &format;               break;
          case 'l': dir    = (int)strtol(s, &s, 0); break;
          case 'q': sort   = (int)strtol(s, &s, 0); break;
          case 'u': filter =      strtod(s, &s);    break;
          case 'h': tree   = 0;                     break;
          case 'j': heap   = 0;                     break;
          case 'x': mode  &= ~IST_PERFECT;          break;
          case 'y': post   = 1;                     break;
          case 'b': optarg = &blanks;               break;
          case 'f': optarg = &fldseps;              break;
          case 'r': optarg = &recseps;              break;
          case 'C': optarg = &comment;              break;
          default : error(E_OPTION, *--s);          break;
        }                       
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       
    else {                     
      switch (k++) {            
        case  0: fn_in  = s;      break;
        case  1: fn_out = s;      break;
        case  2: fn_app = s;      break;
        default: error(E_ARGCNT); break;
      }                         
    }
  }
  if (optarg) error(E_OPTARG);  
  if ((k < 2) || (k > 3))       
    error(E_ARGCNT);           
  if ((!fn_in || !*fn_in) && (fn_app && !*fn_app))
    error(E_STDIN);             
  switch (target) {             
    case 's': target = TT_SET;               break;
    case 'c': target = TT_CLOSED;            break;
    case 'm': target = TT_MAXIMAL;           break;
    case 'r': target = TT_RULE;              break;
    default : error(E_TARGET, (char)target); break;
  }
  if (min < 0) error(E_SIZE, min); 
  if (max < 0) error(E_SIZE, max); 
  if (supp  > 1)                
    error(E_SUPP, supp);        
  if ((conf  < 0) || (conf > 1))
    error(E_CONF, conf);       
  switch (eval) {              
    case 'x': case 0: eval = IST_NONE;      break;
    case 'c': eval = IST_CONF;              break;
    case 'd': eval = IST_DIFF;              break;
    case 'l': eval = IST_LIFT;              break;
    case 'a': eval = IST_LD21;              break;
    case 'q': eval = IST_QUOT;              break;
    case 'n': eval = IST_CHI2;              break;
    case 'p': eval = IST_PVAL;              break;
    case 'i': eval = IST_INFO;              break;
    case 'g': eval = IST_PGST;              break;
    case 'b': eval = IST_LOGQ;              break;
    default : error(E_MEASURE, (char)eval); break;
  }
  switch (aggm) {            
    case 'x': case 0: aggm = IST_NONE;      break;
    case 'm': aggm = IST_MIN;               break;
    case 'n': aggm = IST_MAX;               break;
    case 'a': aggm = IST_AVG;               break;
    default : error(E_MEASURE, (char)aggm); break;
  }
  if ((target > TT_SET)         
  || ((eval > IST_NONE) && (eval < IST_LOGQ)))
    mode &= ~IST_PERFECT;      
  if (target <= TT_MAXIMAL) {  
    mode |= APP_BOTH; conf = 1;}
  if ((filter <= -1) || (filter >= 1))
    filter = 0;                 

  ibase = ib_create(-1);       
  if (!ibase) error(E_NOMEM);  
  ib_chars(ibase, blanks, fldseps, recseps, comment);
  MSG(stderr, "\n");          

  if (fn_app) {                
    t = clock();                
    if (*fn_app)            
      in = fopen(fn_app, "r");  
    else {                      
      in = stdin; fn_app = "<stdin>"; }   
    MSG(stderr, "reading %s ... ", fn_app);
    if (!in) error(E_FOPEN, fn_app);
    k = ib_readapp(ibase, in); 
    if (k  != 0) error(k, fn_app, RECCNT(ibase), BUFFER(ibase));
    if (in != stdin) fclose(in);
    in = NULL;                  
    MSG(stderr, "[%d item(s)]", ib_cnt(ibase));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                           

  t = clock();                 
  if (fn_in && *fn_in)        
    in = fopen(fn_in, "r");     
  else {                       
    in = stdin; fn_in = "<stdin>"; }   
  MSG(stderr, "reading %s ... ", fn_in);
  if (!in) error(E_FOPEN, fn_in);
  tabag = tb_create(ibase);     
  if (!tabag) error(E_NOMEM);   
  while (1) {                
    k = ib_read(ibase, in);     
    if (k) { if (k > 0) break; 
      error(k, fn_in, RECCNT(ibase), BUFFER(ibase)); }
    if (tb_add(tabag, NULL) != 0) error(E_NOMEM);
  }                            
  if (in != stdin) fclose(in);  
  in  = NULL;                  
  n   = ib_cnt(ibase);          
  k   = tb_cnt(tabag);         
  wgt = tb_wgt(tabag);          
  MSG(stderr, "[%d item(s), ", n);
  if (k == wgt) MSG(stderr,    "%d transaction(s)]", k);
  else          MSG(stderr, "%d/%d transaction(s)]", k, wgt);
  MSG(stderr, " done [%.2fs].", SEC_SINCE(t));
  if ((n <= 0) || (wgt <= 0))  
    error(E_NOTRANS);           
  MSG(stderr, "\n");            
  if (format == dflt) {       
    if (target != TT_RULE) format = (supp < 0) ? "  (%a)" : "  (%1S)";
    else format = (supp < 0) ? "  (%b, %1C)" : "  (%1X, %1C)";
  }                            
  supp = ceil (((supp < 0) ? -100 : wgt) *supp);
  smax = floor(((smax < 0) ? -100 : wgt) *smax);

  
  t = clock();                  
  MSG(stderr, "filtering, sorting and recoding items ... ");
  map = (int*)malloc(n *sizeof(int));
  if (!map) error(E_NOMEM);     
  k = (int)((mode & APP_HEAD) ? supp : ceil(supp *conf));
  n = ib_recode(ibase, k, sort, map);
  tb_recode(tabag, map);       
  tb_itsort(tabag, 1, heap);    
  free(map); map = NULL;        
  MSG(stderr, "[%d item(s)] done [%.2fs].", n, SEC_SINCE(t));
  if (n <= 0) error(E_NOFREQ); 
  MSG(stderr, "\n");            
  k   = tb_max(tabag);         
  if (max > k) max = k;         


  t = clock();                  
  MSG(stderr, "reducing transactions ... ");
  tb_filter(tabag, min, NULL);  
  tb_sort(tabag, 1, heap);      
  k = tb_reduce(tabag);         
  if (k == wgt) MSG(stderr,    "[%d transaction(s)]", k);
  else          MSG(stderr, "[%d/%d transaction(s)]", k, wgt);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  
  tt = 0;                      
  if (tree) {                   
    t = clock();               
    MSG(stderr, "building transaction tree ... ");
    tatree = tt_create(tabag);  
    if (!tatree) error(E_NOMEM);
    if (filter == 0) {          
      tb_delete(tabag, 0);      
      tabag = NULL;             
    }
    MSG(stderr, "[%d node(s)]", tt_nodecnt(tatree));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
    tt = clock() -t;            
  }                             


  t = clock(); tc = 0;         
  istree = ist_create(ibase, mode, (int)supp, (int)smax, conf);
  if (!istree) error(E_NOMEM);  
  ist_seteval(istree, eval, aggm, minval, prune);

  /* --- check item subsets --- */
  MSG(stderr, "checking subsets of size 1");
  map = (int*)malloc(n *sizeof(int));
  if (!map) error(E_NOMEM);     
  while (1) {                   
    size = ist_height(istree);  
    if (size >= max) break;     
    if ((filter != 0)        
    &&  (ist_check(istree, map) <= size))
      break;                  
    if (post)                  
      ist_prune(istree);       
    k = ist_addlvl(istree);     
    if (k) { if (k > 0) break;
             error(E_NOMEM);  } 
    if (((filter < 0)           
    &&   (i < -filter *n))      
    ||  ((filter > 0)          
    &&   (i < n) && (i *(double)tt < filter *n *tc))) {
      n = i;                   
      x = clock();             
      tb_filter(tabag, size+1, map);
      tb_sort(tabag, 0, heap);  
      tb_reduce(tabag);         
      if (tatree) {             
        tt_delete(tatree, 0);   
        tatree = tt_create(tabag);
        if (!tatree) error(E_NOMEM);
      }                         
      tt = clock() -x;          
    }
    MSG(stderr, " %d", ++size); 
    x = clock();             
    if (tatree) ist_countx(istree, tatree);
    else        ist_countb(istree, tabag);
    tc = clock() -x;           
  }                             
  free(map); map = NULL;        
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  if ((target == TT_CLOSED) || (target == TT_MAXIMAL)) {
    t = clock();               
    MSG(stderr, "filtering for %s item sets ... ",
        (target == TT_MAXIMAL) ? "maximal" : "closed");
    k = target | ((prune < 0) ? IST_EVAL : 0);
    ist_mark(istree, k);        
    MSG(stderr, "done [%.2fs].\n", SEC_SINCE(t));
  }      

  t = clock();                  
  if (fn_out && *fn_out)        
    out = fopen(fn_out, "w");  
  else {                        
    out = stdout; fn_out = "<stdout>"; }    
  MSG(stderr, "writing %s ... ", fn_out);
  if (!out) error(E_FOPEN, fn_out);
  if (eval == IST_LOGQ) report |= ISR_LOGS;
  if ((target == TT_CLOSED) || (target == TT_MAXIMAL))
    report |= ISR_CLOSED;      
  isrep = isr_create(ibase, out, report, isep, impl);
  if (!isrep) error(E_NOMEM);  
  isr_setfmt (isrep, format);   
  isr_setsize(isrep,  min, max);
  ist_setsize(istree, min, max, dir);
  ist_init   (istree);          
  items = t_items(ib_tract(ibase));
  if ((target <= TT_MAXIMAL)    
  &&  (dir == 0)) {            
    if      (eval == IST_LOGQ)  
      isr_seteval(isrep, isr_logq,  NULL,   minval);
    else if (eval >  IST_NONE)  
      isr_seteval(isrep, ist_evalx, istree, minval);
    n = ist_report(istree, isrep); } 
  else if (target <= TT_MAXIMAL) { 
    for (n = 0; 1; ) {         
      k = ist_set(istree, items, &frq, &minval);
      if (k < 0) break;         
      if (k > 0) fputs(isr_name(isrep, items[0]), out);
      for (i = 0; ++i < k; ) {  
        fputs(isep, out); fputs(isr_name(isrep, items[i]), out); }
      if (format)               
        isr_sinfo(isrep, frq, minval);
      fputc('\n', out); n++;    
    } }                       
  else if (target == TT_RULE) { 
    for (n = 0; 1; ) {          
      k = ist_rule(istree, items, &frq, &body, &head, &minval);
      if (k < 0) break;         
      fputs(isr_name(isrep, items[0]), out);
      fputs(impl, out);         
      if (k > 1) fputs(isr_name(isrep, items[1]), out);
      for (i = 1; ++i < k; ) {  
        fputs(isep, out); fputs(isr_name(isrep, items[i]), out); }
      if (format)               
        isr_rinfo(isrep, frq, body, head, minval);
      fputc('\n', out); n++;    
    }                           
  }  /
  if (fflush(out) != 0) error(E_FWRITE, fn_out);
  if (out != stdout) fclose(out);
  out = NULL;                   
  MSG(stderr, "[%d %s(s)] done ", n,
              (target == TT_RULE) ? "rule" : "set");
  MSG(stderr, "[%.2fs].\n", SEC_SINCE(t));
  #ifdef BENCH
  printf("number of created nodes    : %d\n", istree->ndcnt);
  printf("number of pruned  nodes    : %d\n", istree->ndprn);
  printf("number of item map elements: %d\n", istree->mapsz);
  printf("number of support counters : %d\n", istree->sccnt);
  printf("necessary support counters : %d\n", istree->scnec);
  printf("pruned    support counters : %d\n", istree->scprn);
  printf("number of child pointers   : %d\n", istree->cpcnt);
  printf("necessary child pointers   : %d\n", istree->cpnec);
  printf("pruned    child pointers   : %d\n", istree->cpprn);
  #endif

 
  #ifndef NDEBUG                
  isr_delete(isrep, 0);            
  ist_delete(istree);              
  if (tatree) tt_delete(tatree, 0);
  if (tabag)  tb_delete(tabag, 0); 
  ib_delete(ibase);              
  #endif
  #ifdef STORAGE                
  showmem("at end of program");
  #endif
  return 0;                     
}  
