#-----------------------------------------------------------------------
# File    : makefile   (directory: tract)
#-----------------------------------------------------------------------
CC      = gcc
CFBASE  = -ansi -Wall -pedantic -I$(UTILDIR) $(ADDFLAGS)
CFLAGS  = $(CFBASE) -DNDEBUG -O3
# CFLAGS  = $(CFBASE) -DNDEBUG -O3 -DARCH64
# CFLAGS  = $(CFBASE) -g
# CFLAGS  = $(CFBASE) -g -DARCH64
# CFLAGS  = $(CFBASE) -g -DSTORAGE $(ADDINC)
LDFLAGS =
LIBS    = -lm
# ADDINC  = -I../../misc/src
# ADDOBJ  = storage.o

UTILDIR = ../../util/src

#-----------------------------------------------------------------------
# Build Module
#-----------------------------------------------------------------------
all:       tract.o report.o

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tract.o:   tract.h $(UTILDIR)/symtab.h
tract.o:   tract.c makefile
	$(CC) $(CFLAGS) -c tract.c -o $@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
report.o:  report.h tract.h $(UTILDIR)/symtab.h
report.o:  report.c makefile
	$(CC) $(CFLAGS) -c report.c -o $@

#-----------------------------------------------------------------------
# Storage Debugging
#-----------------------------------------------------------------------
storage.o:  ../../misc/src/storage.h
storage.o:  ../../misc/src/storage.c
	$(CC) $(CFLAGS) -c ../../misc/src/storage.c -o $@

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	rm -f *.o *~ *.flc core
