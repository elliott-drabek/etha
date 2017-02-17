# -*- makefile -*-
# just make headings for each thing and put the headings
# on the line for ALL
PROJECTDIR = .
BINDIR = .

CC = gcc
CPPC = g++
CFLAGS = -g -Wall
LDFLAGS = -lm

DEPEND_FILES = *.cc *.c *.h
CLEANABLE_FILES = *.o *~

ALL = kmer-repair multi-trace multi-walk primer-pair-matches unitig

.SUFFIXES: .cc .c

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o 

.cc.o:
	$(CPPC) $(CFLAGS) -c $*.cc -o $*.o 

all:	$(ALL)

kmer-repair:	kmer-repair.o delcher.o fasta.o kmer-hash.o
	$(CPPC) -o $(BINDIR)/$@ kmer-repair.o delcher.o fasta.o kmer-hash.o $(LDFLAGS)

multi-trace:	multi-trace.o delcher.o fasta.o kmer-hash.o
	$(CPPC) -o $(BINDIR)/$@ multi-trace.o delcher.o fasta.o kmer-hash.o $(LDFLAGS)

multi-walk:	multi-walk.o delcher.o fasta.o kmer-hash.o
	$(CPPC) -o $(BINDIR)/$@ multi-walk.o delcher.o fasta.o kmer-hash.o $(LDFLAGS)

primer-pair-matches:	primer-pair-matches.o delcher.o fasta.o
	$(CPPC) -o $(BINDIR)/$@ primer-pair-matches.o delcher.o fasta.o kmer-hash.o $(LDFLAGS)

unitig:	unitig.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ unitig.o delcher.o $(LDFLAGS)

depend:
	makedepend $(DEPEND_FILES)

clean:
	/bin/rm -f $(CLEANABLE_FILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.



