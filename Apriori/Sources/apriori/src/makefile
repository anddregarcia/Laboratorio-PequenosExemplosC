#-----------------------------------------------------------------------
# File    : makefile
#
#Projeto de Organização e recuperação da Informação (ORI)
#Escola de Engenharia de Piracicaba
#Authors: Adilson Perecin(a.perecin(at)hotmail.com) & #Cristiano Benato(benato(at)hst.com.br)
#
#Tema: Organização e verificação de itemsets e ordenação de #coincidencias através
#do algoritmo de ordenação Apriori
#
#Instructor: Luiz Camolesi
# 
#
#
#-----------------------------------------------------------------------
CC       = gcc
CFBASE   = -ansi -Wall -pedantic $(ADDFLAGS) \
           -I$(UTILDIR) -I$(MATHDIR) -I$(TRACTDIR)
CFLAGS   = $(CFBASE) -DNDEBUG -O3
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -DBENCH
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -DARCH64
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DARCH64
# CFLAGS   = $(CFBASE) -g -DSTORAGE $(ADDINC)
LDFLAGS  =
LIBS     = -lm
# ADDINC   = -I../../misc/src
# ADDOBJ   = storage.o

UTILDIR  = ../../util/src
MATHDIR  = ../../math/src
TRACTDIR = ../../tract/src

HDRS     = $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
           $(UTILDIR)/tabscan.h $(UTILDIR)/scan.h \
           $(MATHDIR)/gamma.h   $(MATHDIR)/chi2.h \
           $(TRACTDIR)/tract.h  $(TRACTDIR)/report.h \
           istree.h
OBJS     = $(UTILDIR)/arrays.o  $(UTILDIR)/nimap.o \
           $(UTILDIR)/tabscan.o $(UTILDIR)/scform.o \
           $(MATHDIR)/gamma.o   $(MATHDIR)/chi2.o \
           $(TRACTDIR)/tract.o  $(TRACTDIR)/report.o \
           istree.o apriori.o $(ADDOBJ)
PRGS     = apriori

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:       $(PRGS)

apriori:   $(OBJS) makefile
	$(CC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
apriori.o: $(HDRS)
apriori.o: apriori.c makefile
	$(CC) $(CFLAGS) -c apriori.c -o $@

#-----------------------------------------------------------------------
# Frequent Item Set Tree Management
#-----------------------------------------------------------------------
istree.o:  $(HDRS)
istree.o:  istree.c makefile
	$(CC) $(CFLAGS) -c istree.c -o $@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o  ADDFLAGS=$(ADDFLAGS)
$(UTILDIR)/nimap.o:
	cd $(UTILDIR);  $(MAKE) nimap.o   ADDFLAGS=$(ADDFLAGS)
$(UTILDIR)/tabscan.o:
	cd $(UTILDIR);  $(MAKE) tabscan.o ADDFLAGS=$(ADDFLAGS)
$(UTILDIR)/scform.o:
	cd $(UTILDIR);  $(MAKE) scform.o  ADDFLAGS=$(ADDFLAGS)
$(MATHDIR)/gamma.o:
	cd $(MATHDIR);  $(MAKE) gamma.o   ADDFLAGS=$(ADDFLAGS)
$(MATHDIR)/chi2.o:
	cd $(MATHDIR);  $(MAKE) chi2.o    ADDFLAGS=$(ADDFLAGS)
$(TRACTDIR)/tract.o:
	cd $(TRACTDIR); $(MAKE) tract.o   ADDFLAGS=$(ADDFLAGS)
$(TRACTDIR)/report.o:
	cd $(TRACTDIR); $(MAKE) report.o  ADDFLAGS=$(ADDFLAGS)

#-----------------------------------------------------------------------
# Storage Debugging
#-----------------------------------------------------------------------
storage.o:  ../../misc/src/storage.h
storage.o:  ../../misc/src/storage.c
	$(CC) $(CFLAGS) -c ../../misc/src/storage.c -o $@

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	cd ../..; rm -f apriori.zip apriori.tar.gz; \
        zip -rq apriori.zip    apriori/{src} \
          tract/{src} util/{src} math/{src}; \
        tar cfz apriori.tar.gz apriori/{src} \
          tract/{src} util/{src} math/{src}

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) $(HOME)/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.o *~ *.flc core $(PRGS)

clean:
	$(MAKE) localclean
	cd $(UTILDIR);  $(MAKE) clean
	cd $(MATHDIR);  $(MAKE) clean
	cd $(TRACTDIR); $(MAKE) clean
