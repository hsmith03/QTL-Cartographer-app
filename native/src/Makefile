#
# Makefile for the QTL Cartographer suite
# Written by Chris Basten  2 March 1995
# This version last edited 2 October 2002
#
#    Look for the word MODIFY to indicate what you might need to change.
#
#  SOURCES are C source files
#  HEADERS are header files
#  NOCOPY are files that don't get a copyright notice
#  OTHERS are the Makefile and conversion script

PUBSRC=Blas.c Linpak.c NumRec.c LAPACK.c f.c ranlib.c 

BWZSRC= Otraits.c  RCmain.c   Zfunc.c  EQfunc.c   Preplot.c  RMfunc.c   Zmain.c   \
	EQmain.c   MLnpkws.c  Prune.c    RMmain.c  Genome.c   MZfunc.c   QSfunc.c   RQfunc.c  \
	Idatain.c  MZmain.c   QSmain.c   RQmain.c  LRfunc.c   Mdatain.c  QTLcart.c  SRfunc.c  \
	LRmain.c   MissMark.c Qdatain.c  SRmain.c  RCfunc.c   Utilities.c  \
	MImain.c   MIfunc.c   JZmain.c   params.c  MRmain.c MRfunc.c   Emain.c    Efunc.c    

MWBSRC= BTL_LRfunc.c BTL_LRmain.c

GYBSRC= Bmain.c accept_birth.c birth.c death.c INITVALS.c mcmc.c QTL_MOVE.c read_data.c\
	update_parms.c update_qtl.c

SOURCES= $(BWZSRC) $(MWBSRC) $(GYBSRC) $(PUBSRC)

PUBHEAD=NumRec.h Blas.h Linpak.h f2c.h ranlib.h  

BWZHEAD=Mdatain.h  Qdatain.h  Utilities.h \
	EQfunc.h   LocalD.h   MissMark.h RCfunc.h   Zfunc.h     \
	Genome.h   RMfunc.h   params.h    \
	Idatain.h  MZfunc.h   Otraits.h  RQfunc.h   MIfunc.h    \
	LRfunc.h   Main.h     QSfunc.h   SRfunc.h   MRfunc.h Efunc.h

MWBHEAD=BTL_LRfunc.h

GYBHEAD=chrom.h qtl.h  revjump.h mygenome.h

HEADERS=$(BWZHEAD) $(MWBHEAD) $(GYBHEAD) $(PUBHEAD)

NOCOPY=Blas.h Linpak.h Blas.c Linpak.c NumRec.c NumRec.h BTL_LRfunc.c \
	BTL_LRmain.c BTL_LRfunc.h

OTHERS=LocalVert.pl Makefile MemMatch.pl Test.csh

#
#####################################################################
#                                    MODIFY
#####################################################################
#   LDFLAGS   for static binaries, add -Bstatic (probably obsolete)
#   CPPFLAGS  use -DAQUA  for MacOS X only
#   CFLAGS    Sometimes -O  will create buggy binaries.  You may 
#               need to eliminate it.
#   CC        Your compiler may be gcc, or you may need to 
#               specify a complete path.   
#   debug     Only works with Sun make, so it is commented out here.  
#debug := CFLAGS = -g
#  These work for Mac OSX 10.2 and higher.
LDFLAGS= -lm 
CPPFLAGS=-DUNIX -DITOA -DAQUA
CFLAGS=-O   
CC=cc
#  for Solaris 7 and higher, if you want 64 bit binaries using the SUNWspro
# compiler, then use
#CPPFLAGS=-DUNIX -DITOA 
#CC=cc
#LDFLAGS= -xarch=v9 -lm 
#CFLAGS=-xarch=v9 -O   
#  you can compile 32 bit binaries by eliminating the -xarch=v9 in the above
#  Someone used these options on a Unix Tru64 v5.1 machine
#CFLAGS=-O3 -ieee -arch ev7 -tune ev7
#
# I use cp for installation because install doesn't always work the same 
# way on different systems.   You may also want to change the OWNER and GROUP
# Set the BINDIR and LBINDIR.  BINDIR is where binaries are installed.  LBINDIR
# is a directory for binaries using the -g compile flag, and you don't need to 
# worry about it.    You may want to use /usr/local/bin for BINDIR so that
# everyone has access to the programs.
#
INSTALL = /bin/cp
REMOVE = /bin/rm -f
OWNER = bin
GROUP = bin
BINDIR = ~/bin
LBINDIR = ../bin
###################################################################
#                                 END  MODIFY
#################################################################### 
#
#
LINTFILES= $(SOURCES:.c=.ln)
OBJECTS= $(SOURCES:.c=.o)
TARGETS= Rmap Rqtl Eqtl Rcross Emap Qstats LRmapqtl BTmapqtl SRmapqtl Zmapqtl \
	Preplot Prune JZmapqtl MImapqtl MultiRegress Bmapqtl
LINTFLAGS= -chu  
all debug : $(TARGETS)
$(OBJECTS) : $(HEADERS)         # objects depend on all headers.
#
.KEEP_STATE:
#
# A rule for each target:
#
Preplot: Preplot.o EQfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o
	$(CC) -o Preplot Preplot.o EQfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

QTLcart : $(OBJECTS)
	$(CC) -o QTLcart QTLcart.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o $(LDFLAGS)

Rmap: RMmain.o RMfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o
	$(CC) -o Rmap RMmain.o RMfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

Rqtl: RQmain.o RQfunc.o   params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o
	$(CC) -o Rqtl RQmain.o RQfunc.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o $(LDFLAGS)

Eqtl: EQmain.o EQfunc.o MIfunc.o MissMark.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o 
	$(CC) -o Eqtl EQmain.o EQfunc.o MIfunc.o MissMark.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

Rcross: RCmain.o RCfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o 
	$(CC) -o Rcross RCmain.o RCfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o $(LDFLAGS)

Emap:  Emain.o Efunc.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MissMark.o Blas.o
	$(CC) -o Emap Emain.o Efunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MissMark.o Blas.o $(LDFLAGS) 

Prune: Prune.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o
	$(CC) -o Prune Prune.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

Qstats: QSmain.o QSfunc.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o 
	$(CC) -o Qstats QSmain.o QSfunc.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o $(LDFLAGS)

LRmapqtl: LRmain.o LRfunc.o Blas.o Linpak.o MLnpkws.o Otraits.o MissMark.o QSfunc.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o
	$(CC) -o LRmapqtl LRmain.o LRfunc.o Blas.o Linpak.o MLnpkws.o Otraits.o  MissMark.o QSfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

BTmapqtl: BTL_LRmain.o BTL_LRfunc.o LRfunc.o  Blas.o Linpak.o MLnpkws.o Otraits.o MissMark.o QSfunc.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o 
	$(CC) -o BTmapqtl BTL_LRmain.o BTL_LRfunc.o LRfunc.o Blas.o Linpak.o MLnpkws.o Otraits.o MissMark.o QSfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o  $(LDFLAGS)

Zmapqtl: Zmain.o Zfunc.o QSfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o   
	$(CC) -o Zmapqtl Zmain.o Zfunc.o QSfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o    $(LDFLAGS)

SRmapqtl: SRmain.o SRfunc.o Otraits.o  MissMark.o  params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o 
	$(CC) -o SRmapqtl SRmain.o SRfunc.o Otraits.o  MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o  $(LDFLAGS)

JZmapqtl: MZmain.o MZfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o  
	$(CC) -o JZmapqtl MZmain.o MZfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o $(LDFLAGS)

JZtest: JZmain.o MZfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o  
	$(CC) -o JZtest JZmain.o MZfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o Linpak.o MLnpkws.o Blas.o $(LDFLAGS)

MImapqtl: MImain.o MIfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MLnpkws.o 
	$(CC) -o MImapqtl MImain.o MIfunc.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MLnpkws.o  $(LDFLAGS)

MultiRegress: MRmain.o MRfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MLnpkws.o Linpak.o Blas.o
	$(CC) -o MultiRegress MRmain.o MRfunc.o params.o Utilities.o Mdatain.o Genome.o Qdatain.o Idatain.o NumRec.o MLnpkws.o Linpak.o Blas.o $(LDFLAGS)

Bmapqtl: Bmain.o QTL_MOVE.o birth.o mcmc.o update_parms.o INITVALS.o accept_birth.o death.o read_data.o  update_qtl.o \
	LAPACK.o f.o ranlib.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Idatain.o NumRec.o RCfunc.o   
	$(CC) -o Bmapqtl Bmain.o QTL_MOVE.o birth.o mcmc.o update_parms.o INITVALS.o accept_birth.o death.o read_data.o  update_qtl.o \
	LAPACK.o f.o ranlib.o Otraits.o MissMark.o params.o Utilities.o Mdatain.o Genome.o Idatain.o NumRec.o RCfunc.o   $(LDFLAGS)


lint: $(LINTFILES)
	$(LINT.c) $(LINTFILES)

QTLCart.tar : $(SOURCES) $(OTHERS) $(HEADERS)
	tar cvf QTLCart.tar $(SOURCES) $(OTHERS) $(HEADERS) 

QTLCart.tar.Z : QTLCart.tar
	compress QTLCart.tar

QTLCart.tar.gz : QTLCart.tar
	gzip QTLCart.tar

#  Note that you may need root access to chown and chgrp
install: all
	for i in ${TARGETS}; do ( $(INSTALL) $$i $(BINDIR) ); done
	for i in ${TARGETS}; do ( chown $(OWNER) $(BINDIR)/$$i ); done
	for i in ${TARGETS}; do ( chgrp $(GROUP) $(BINDIR)/$$i ); done

linstall: all
	for i in ${TARGETS}; do ( $(INSTALL) $$i $(LBINDIR) ); done

install2: debug  
	for i in ${TARGETS}; do ( $(INSTALL) $$i $(LBINDIR) ); done

uninstall:
	cd $(BINDIR) ; for i in ${TARGETS};  do ( rm -f $$i ); done

uninstall2: 
	cd $(LBINDIR) ; for i in ${TARGETS};  do ( rm -f $$i ); done

# pretty-print format the source files, but make an archive for safety first.
pretty: tar
	for i in ${SOURCES}; do ( indent $$i $$i.t -l100 -lc60 -nbc -nce -nfc1 -nps1 -i2 -br -v ; mv $$i.t $$i ); done

clean:
	rm -f $(LINTFILES) $(OBJECTS) *.BAK core .core ..core

veryclean: clean
	rm -f *.err *.out *% *.i.c *.a.c lint.out *.com $(TARGETS)

cleantargets :
	rm -f $(TARGETS)

#
#  The following targets only make sense from the original distribution subdirectory on 
#  the Mac OSX machine where development takes place.   If you know Perl and how to modify
#  makefiles, then you can read the LocalVert.pl script and determine how to convert all
#  the sources and headers to Macintosh or Dos line endings.   
#
UNIXDIR = ../QTLCartUnix/src
DOSDIR = ./winsrc

unix-src: $(SOURCES) $(HEADERS) $(OTHERS)
	echo " " ; echo "       Removing old UNIX Sources " ; echo " "
	for i in ${SOURCES}; do ( $(REMOVE) $(UNIXDIR)/$$i ); done
	for i in ${HEADERS}; do ( $(REMOVE) $(UNIXDIR)/$$i ); done
	for i in ${OTHERS}; do ( $(REMOVE) $(UNIXDIR)/$$i ); done
	echo " " ; echo "       Copying new  UNIX Sources " ; echo " " 
	for i in ${PUBSRC};  do ( ./LocalVert.pl -s 0 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${PUBHEAD}; do ( ./LocalVert.pl -s 0 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${BWZSRC};  do ( ./LocalVert.pl -s 1 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${BWZHEAD}; do ( ./LocalVert.pl -s 1 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${MWBSRC};  do ( ./LocalVert.pl -s 2 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${MWBHEAD}; do ( ./LocalVert.pl -s 2 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${GYBSRC};  do ( ./LocalVert.pl -s 3 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${GYBHEAD}; do ( ./LocalVert.pl -s 3 -c -i m -o u -f $$i <  $$i > $(UNIXDIR)/$$i ); done
	for i in ${OTHERS}; do ( $(INSTALL) $$i   $(UNIXDIR)/$$i ); done	
	echo " " ; echo "       Finished updating UNIX Sources " ; echo " "

dos-src:
	echo " " ; echo "       Removing old DOS Sources " ; echo " "
	for i in ${SOURCES}; do ( $(REMOVE) $(DOSDIR)/$$i ); done
	for i in ${HEADERS}; do ( $(REMOVE) $(DOSDIR)/$$i ); done
	for i in ${OTHERS}; do ( $(REMOVE) $(DOSDIR)/$$i ); done
	echo " " ; echo "       Copying new DOS Sources " ; echo " "
	for i in ${PUBSRC};  do ( ./LocalVert.pl -s 0 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${PUBHEAD}; do ( ./LocalVert.pl -s 0 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${BWZSRC};  do ( ./LocalVert.pl -s 1 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${BWZHEAD}; do ( ./LocalVert.pl -s 1 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${MWBSRC};  do ( ./LocalVert.pl -s 2 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${MWBHEAD}; do ( ./LocalVert.pl -s 2 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${GYBSRC};  do ( ./LocalVert.pl -s 3 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${GYBHEAD}; do ( ./LocalVert.pl -s 3 -c -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	for i in ${OTHERS}; do ( ./LocalVert.pl -i m -o d -f $$i <  $$i > $(DOSDIR)/$$i ); done
	echo " " ; echo "       Finished updating DOS Sources " ; echo " "

unix:  unix-src
	echo " " ; echo "       Changing to UNIX directory to make binaries " ; echo " "
	cd $(UNIXDIR) ; make install

