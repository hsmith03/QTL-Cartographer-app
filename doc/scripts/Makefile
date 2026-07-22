#
#  Makefile for the  scripts
#  Chris Basten, 17 May 2001
#
.SUFFIXES: .txt .pl .v .ps .1 .print .html .tex .php

MAKEFILE        = Makefile

all: install

# ======================>  Check these for your system
REMOVE = /bin/rm -f
CAT = /bin/cat
INSTALL = /bin/cp
MOVE = /bin/mv -f
PERL = /usr/bin/perl
CSH = /bin/csh
TOUCH = /usr/bin/touch
PAPERSIZE = letter
RCP = /bin/rcp
WEBDIR = statgen.ncsu.edu:/local/www/qtlcart
#  I'm guessing that you could use a4 instead of letter.
#  These are defined for making the distribution: do not change
UNIXBASE = ../../QTLCartUnix/doc
MACBASE = ../../QTLCartMac/doc
WINBASE = ../../QTLCartWin/doc
TEXDIR = ../tex
# ======================>   Ultimate location of the scripts
#   Change this to suit your needs...it should be the same as 
#   where the QTL Cart binaries go.
BINDIR =  /Users/basten/bin

# ============> Probably no need to change anything below.
#  
#  Here are the rules for converting the perl scripts (containing 
#  pod documentation) into various formats.
#
#  wrap php around html file 
.html.php :
	./phpwrap.pl < $< > $*.php
	
#  into plain text
.pl.txt :
	pod2text $< > $*.txt

# into nroff format
.pl.1 :
	pod2man --center="QTL Cartographer 1.17 Perl Script" $< > $*.1

#  into latex2e  pod2latex appends .tex to the .pl file, so move it back
.pl.tex :
	pod2latex $<  
#	$(MOVE) $*.pl.tex $*.tex

# into html
.pl.html :
	pod2html --noindex --title=$*  $< > $*.html

# nroff the man page and pip it into less (change to more if you like) 
.1.v :
	nroff -man $< | less

# nroff the man page and pipe it through a2ps to lpr
.1.print :
	nroff -man $< | a2ps -m --medium=$(PAPERSIZE) --portrait --columns 1 

# nroff the man page and pipe it through a2ps to a file
.1.ps :
	nroff -man $< | a2ps -m --medium=$(PAPERSIZE) --portrait --columns 1 --output=- > $*.ps

# Define the source files
PLS =   SSupdate.pl Prepraw.pl EWThreshold.pl GetMaxLR.pl Permute.pl \
	Vert.pl CWTupdate.pl Ztrim.pl SRcompare.pl Model8.pl Bootstrap.pl \
	Dobasics.pl Ttransform.pl TestExamples.pl RCpermute.pl SumLR.pl CSSupdate.pl

OTHERS = COPYING Makefile INSTALL scripts.cfg

SCRIPTS = Bootstrap.csh Permute.csh PermuteMI.csh PermuteJZ.csh  Model8.csh Sample.csh Test.csh

VPAGES = $(MAN:.1=.v)

HPAGES = $(MAN:.1=.print)

PSPAGES = $(MAN:.1=.ps)

MAN = $(PLS:.pl=.1)

HTMLS = $(PLS:.pl=.html)

TEXS = $(PLS:.pl=.tex)

TEXTS = $(PLS:.pl=.txt)

PHPS = $(HTMLS:.html=.php)

hardcopy: ${HPAGES}

# Create a local, usable copy of Vert.pl so that we may convert
# all unix formatted files into macintosh and dos formats.
Vrt: Vert.pl
	$(REMOVE) Vrt
	$(CAT) scripts.cfg Vert.pl > Vrt
	chmod +x Vrt
	cp Vrt ..	

textfiles: $(TEXTS) Vrt
	echo "==========>  Installing text documentation for scripts"
	$(INSTALL) $(TEXTS) $(UNIXBASE)/txt
	for i in ${TEXTS}; do ( $(REMOVE) $(MACBASE)/txt/$$i ); done
	for i in ${TEXTS}; do ( ../Vrt -i u -o m < $$i >  $(MACBASE)/txt/$$i ); done
	for i in ${TEXTS}; do ( $(REMOVE) $(WINBASE)/txt/$$i ); done           
	for i in ${TEXTS}; do ( ../Vrt -i u -o d < $$i >  $(WINBASE)/txt/$$i ); done

#postscript: ${PSPAGES}
#	cp ${PSPAGES} $(PSDIR)

htmls:  $(HTMLS)
	echo "==========>  Installing html documentation for scripts"
	$(INSTALL) $(HTMLS) $(UNIXBASE)/html
	for i in ${HTMLS}; do ( $(REMOVE) $(MACBASE)/html/$$i ); done
	for i in ${HTMLS}; do ( ../Vrt -i u -o m < $$i >  $(MACBASE)/html/$$i ); done
	for i in ${HTMLS}; do ( $(REMOVE) $(WINBASE)/html/$$i ); done
	for i in ${HTMLS}; do ( ../Vrt -i u -o d < $$i >  $(WINBASE)/html/$$i ); done
	
phps: $(PHPS)
	for i in ${PHPS}; do ( $(RCP) $$i $(WEBDIR)/$$i ); done

texs:   $(TEXS)
	$(INSTALL) $(TEXS) $(TEXDIR)

bins:  $(PLS) $(SCRIPTS) $(OTHERS)
	echo "==========>  Installing scripts into distributions"
	$(INSTALL) $(PLS) $(UNIXBASE)/scripts
	$(INSTALL) $(OTHERS) $(UNIXBASE)/scripts
	$(INSTALL) $(SCRIPTS) $(UNIXBASE)/scripts
	for i in ${PLS}; do ( $(REMOVE) $(MACBASE)/scripts/$$i ); done
	for i in ${OTHERS}; do ( $(REMOVE) $(MACBASE)/scripts/$$i ); done
	for i in ${SCRIPTS}; do ( $(REMOVE) $(MACBASE)/scripts/$$i ); done
	for i in ${PLS}; do ( ../Vrt -i u -o m < $$i >  $(MACBASE)/scripts/$$i ); done
	for i in ${OTHERS}; do ( ../Vrt -i u -o m < $$i >  $(MACBASE)/scripts/$$i ); done
	for i in ${SCRIPTS}; do ( ../Vrt -i u -o m < $$i >  $(MACBASE)/scripts/$$i); done
	for i in ${PLS}; do ( $(REMOVE) $(WINBASE)/scripts/$$i ); done
	for i in ${OTHERS}; do ( $(REMOVE) $(WINBASE)/scripts/$$i ); done
	for i in ${SCRIPTS}; do ( $(REMOVE) $(WINBASE)/scripts/$$i ); done
	for i in ${PLS}; do ( ../Vrt -i u -o d < $$i >  $(WINBASE)/scripts/$$i ); done
	for i in ${OTHERS}; do ( ../Vrt -i u -o d < $$i >  $(WINBASE)/scripts/$$i ); done
	for i in ${SCRIPTS}; do ( ../Vrt -i u -o d < $$i >  $(WINBASE)/scripts/$$i ); done

SSupdate : SSupdate.1 SSupdate.html SSupdate.tex SSupdate.txt SSupdate.pl SSupdate.php

Prepraw : Prepraw.1 Prepraw.html Prepraw.tex Prepraw.txt Prepraw.pl Prepraw.php

EWThreshold : EWThreshold.1 EWThreshold.html EWThreshold.tex EWThreshold.txt EWThreshold.pl EWThreshold.php

GetMaxLR : GetMaxLR.1 GetMaxLR.html GetMaxLR.tex GetMaxLR.txt GetMaxLR.pl GetMaxLR.php

Vert : Vert.1 Vert.html Vert.tex Vert.txt Vert.pl Vert.php

CWTupdate : CWTupdate.1 CWTupdate.html CWTupdate.tex CWTupdate.txt CWTupdate.pl CWTupdate.php

SRcompare : SRcompare.1 SRcompare.html SRcompare.tex SRcompare.txt SRcompare.pl SRcompare.php

Model8 : Model8.1 Model8.html Model8.tex Model8.txt Model8.pl Model8.php

Bootstrap : Bootstrap.1 Bootstrap.html Bootstrap.tex Bootstrap.txt Bootstrap.pl Bootstrap.php

Permute : Permute.1 Permute.html Permute.tex Permute.txt Permute.pl Permute.php

Dobasics : Dobasics.1 Dobasics.html Dobasics.tex Dobasics.txt Dobasics.pl Dobasics.php

Ztrim : Ztrim.1 Ztrim.html Ztrim.tex Ztrim.txt Ztrim.pl Ztrim.php

Ttransform : Ttransform.1 Ttransform.html Ttransform.tex Ttransform.txt Ttransform.pl Ttransform.php

TestExamples : TestExamples.1 TestExamples.html TestExamples.tex TestExamples.txt TestExamples.pl TestExamples.php

RCpermute : RCpermute.1 RCpermute.html RCpermute.tex RCpermute.txt RCpermute.pl RCpermute.php

SumLR : SumLR.1 SumLR.html SumLR.tex SumLR.txt SumLR.pl SumLR.php

CSSupdate : CSSupdate.1 CSSupdate.html CSSupdate.tex CSSupdate.txt CSSupdate.pl CSSupdate.php

clean:
	$(REMOVE)  *%   

veryclean: clean 
	$(REMOVE) $(HTMLS) $(TEXTS) $(PHPS)
	$(REMOVE) $(MAN) $(TEXS) $(PSPAGES)
	$(REMOVE) pod2htm*

# making dist propogates the files
dist: $(MAN) textfiles htmls bins
	$(INSTALL) $(MAN) $(UNIXBASE)/man

# making install sets the proper command line in the scripts and 
# copies them to the bindir, which should be set to where the QTL Cart 
# binaries reside.
install: $(SCRIPTS) $(PLS)
	for i in ${PLS}; do ( $(REMOVE) $(BINDIR)/$$i ); done
	for i in ${PLS}; do ( $(CAT) scripts.cfg $$i > $(BINDIR)/$$i); done 
	for i in ${SCRIPTS}; do ( sed s%INSERTCSHCOMMAND%$(CSH)% < $$i > $$i.seded ); done
	for i in ${SCRIPTS}; do ( $(MOVE) $$i.seded $(BINDIR)/$$i); done
	for i in ${PLS}; do ( chmod +x $(BINDIR)/$$i); done
	for i in ${SCRIPTS}; do ( chmod +x $(BINDIR)/$$i); done

#
#  sometimes, you need to touch the files before the make will work.
#
touch:
	$(TOUCH) $(SCRIPTS) $(PLS)

archive: veryclean
	cd .. ; $(REMOVE) scripts.tar ; tar cf scripts.tar scripts
	
###

#	for i in ${PLS}; do ( sed s%INSERTPERLCOMMAND%$(PERL)% < $$i > $$i.seded ); done
#	for i in ${PLS}; do ( $(MOVE) $$i.seded $(BINDIR)/$$i); done 
