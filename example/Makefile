#
#  Makefile to install the examples  for Unix systems
#  Chris Basten,  19 July 2002
#
#        Install the example files
#
FILES = COPYING README  mletest.map mletest.cro mletestq.inp\
	sample.raw sample.mps sample2.raw realdatm.inp realdatc.inp \
	weber519c.inp weber519m.inp weber701c.inp weber701m.inp \
	nuzhdinm.inp nuzhdinc.inp vieiram.inp vieirac.inp \
	lauriem.inp lauriem4c.inp lauriem6c.inp lauries4c.inp lauries6c.inp \
	rauhmap.inp rauhall.inp rauhratio.inp rauhrl.inp rauham.inp rauhrm.inp
#
MAKEFILE	= Makefile
#
#  Make sure that the REMOVE and INSTALL commands are correct,
#  and that the TESTDIR variable is set correctly.   
#  The base TESTDIR must exist prior to running make.
#
REMOVE = /bin/rm -f
INSTALL = /bin/cp -f
TESTDIR = ../qwork

all: prelim mletest sample realdat weber laurie rauh mackay

clean:
	test -d $(TESTDIR)/sim     && $(REMOVE) $(TESTDIR)/sim/*
	test -d $(TESTDIR)/mletest && $(REMOVE) $(TESTDIR)/mletest/*
	test -d $(TESTDIR)/sample  && $(REMOVE) $(TESTDIR)/sample/*
	test -d $(TESTDIR)/realdat && $(REMOVE) $(TESTDIR)/realdat/*
	test -d $(TESTDIR)/weber/519  && $(REMOVE) $(TESTDIR)/weber/519/*
	test -d $(TESTDIR)/weber/701  && $(REMOVE) $(TESTDIR)/weber/701/*
	test -d $(TESTDIR)/laurie/m4  && $(REMOVE) $(TESTDIR)/laurie/m4/*
	test -d $(TESTDIR)/laurie/m6  && $(REMOVE) $(TESTDIR)/laurie/m6/*
	test -d $(TESTDIR)/laurie/s4  && $(REMOVE) $(TESTDIR)/laurie/s4/*
	test -d $(TESTDIR)/laurie/s6  && $(REMOVE) $(TESTDIR)/laurie/s6/*
	test -d $(TESTDIR)/rauh/ratio && $(REMOVE) $(TESTDIR)/rauh/ratio/*
	test -d $(TESTDIR)/rauh/all   && $(REMOVE) $(TESTDIR)/rauh/all/*
	test -d $(TESTDIR)/rauh/rootmass   && $(REMOVE) $(TESTDIR)/rauh/rootmass/*
	test -d $(TESTDIR)/rauh/rootlength && $(REMOVE) $(TESTDIR)/rauh/rootlength/*
	test -d $(TESTDIR)/rauh/arialmass  && $(REMOVE) $(TESTDIR)/rauh/arialmass/*
	test -d $(TESTDIR)/mackay/nuzhdin  && $(REMOVE) $(TESTDIR)/mackay/nuzhdin/*
	test -d $(TESTDIR)/mackay/vieira   && $(REMOVE) $(TESTDIR)/mackay/vieira/*

prelim:
	test -d $(TESTDIR) || mkdir $(TESTDIR)
	test -d $(TESTDIR)/sim || mkdir $(TESTDIR)/sim

mletest: mletest.map mletestq.inp mletest.cro
	test -d $(TESTDIR)/mletest   || mkdir $(TESTDIR)/mletest
	$(REMOVE) $(TESTDIR)/mletest/*
	$(INSTALL) mletest.map   $(TESTDIR)/mletest/mletest.map       
	$(INSTALL) mletestq.inp   $(TESTDIR)/mletest/mletestq.inp       
	$(INSTALL) mletest.cro   $(TESTDIR)/mletest/mletest.cro       
	
sample: sample.mps sample.raw sample2.raw
	test -d $(TESTDIR)/sample || mkdir $(TESTDIR)/sample
	$(REMOVE) $(TESTDIR)/sample/*
	$(INSTALL) sample.mps   $(TESTDIR)/sample/sample.mps       
	$(INSTALL) sample2.raw   $(TESTDIR)/sample/sample2.raw       
	$(INSTALL) sample.raw   $(TESTDIR)/sample/sample.raw       

realdat: realdatm.inp realdatc.inp  
	test -d $(TESTDIR)/realdat || mkdir $(TESTDIR)/realdat
	$(REMOVE) $(TESTDIR)/realdat/*
	$(INSTALL) realdatc.inp   $(TESTDIR)/realdat/realdatc.inp       
	$(INSTALL) realdatm.inp   $(TESTDIR)/realdat/realdatm.inp       

weber: weber519m.inp weber519c.inp  weber701m.inp weber701c.inp
	test -d $(TESTDIR)/weber || mkdir $(TESTDIR)/weber
	test -d $(TESTDIR)/weber/519 || mkdir $(TESTDIR)/weber/519
	$(REMOVE) $(TESTDIR)/weber/519/*
	test -d $(TESTDIR)/weber/701 || mkdir $(TESTDIR)/weber/701
	$(REMOVE) $(TESTDIR)/weber/701/*
	$(INSTALL) weber519c.inp   $(TESTDIR)/weber/519/weber519c.inp       
	$(INSTALL) weber519m.inp   $(TESTDIR)/weber/519/weber519m.inp       
	$(INSTALL) weber701c.inp   $(TESTDIR)/weber/701/weber701c.inp       
	$(INSTALL) weber701m.inp   $(TESTDIR)/weber/701/weber701m.inp       

laurie: lauriem.inp lauriem4c.inp lauriem6c.inp lauries4c.inp lauries6c.inp
	test -d $(TESTDIR)/laurie || mkdir $(TESTDIR)/laurie
	test -d $(TESTDIR)/laurie/m4 || mkdir $(TESTDIR)/laurie/m4
	$(REMOVE) $(TESTDIR)/laurie/m4/*
	test -d $(TESTDIR)/laurie/m6 || mkdir $(TESTDIR)/laurie/m6
	$(REMOVE) $(TESTDIR)/laurie/m6/*
	test -d $(TESTDIR)/laurie/s4 || mkdir $(TESTDIR)/laurie/s4
	$(REMOVE) $(TESTDIR)/laurie/s4/*
	test -d $(TESTDIR)/laurie/s6 || mkdir $(TESTDIR)/laurie/s6
	$(REMOVE) $(TESTDIR)/laurie/s6/*
	$(INSTALL) lauriem.inp   $(TESTDIR)/laurie/m4/lauriem.inp       
	$(INSTALL) lauriem.inp   $(TESTDIR)/laurie/m6/lauriem.inp       
	$(INSTALL) lauriem.inp   $(TESTDIR)/laurie/s4/lauriem.inp       
	$(INSTALL) lauriem.inp   $(TESTDIR)/laurie/s6/lauriem.inp       
	$(INSTALL) lauriem4c.inp   $(TESTDIR)/laurie/m4/lauriem4c.inp       
	$(INSTALL) lauriem6c.inp   $(TESTDIR)/laurie/m6/lauriem6c.inp       
	$(INSTALL) lauries4c.inp   $(TESTDIR)/laurie/s4/lauries4c.inp       
	$(INSTALL) lauries6c.inp   $(TESTDIR)/laurie/s6/lauries6c.inp       

rauh: rauhmap.inp rauhall.inp rauhratio.inp rauhrl.inp rauham.inp rauhrm.inp
	test -d $(TESTDIR)/rauh || mkdir $(TESTDIR)/rauh
	test -d $(TESTDIR)/rauh/all || mkdir $(TESTDIR)/rauh/all
	$(REMOVE) $(TESTDIR)/rauh/all/*
	test -d $(TESTDIR)/rauh/rootmass || mkdir $(TESTDIR)/rauh/rootmass
	$(REMOVE) $(TESTDIR)/rauh/rootmass/*
	test -d $(TESTDIR)/rauh/rootlength || mkdir $(TESTDIR)/rauh/rootlength
	$(REMOVE) $(TESTDIR)/rauh/rootlength/*
	test -d $(TESTDIR)/rauh/ratio || mkdir $(TESTDIR)/rauh/ratio
	$(REMOVE) $(TESTDIR)/rauh/ratio/*
	test -d $(TESTDIR)/rauh/arialmass || mkdir $(TESTDIR)/rauh/arialmass
	$(REMOVE) $(TESTDIR)/rauh/arialmass/*
	$(INSTALL) rauhmap.inp   $(TESTDIR)/rauh/all/rauhmap.inp       
	$(INSTALL) rauhmap.inp   $(TESTDIR)/rauh/ratio/rauhmap.inp       
	$(INSTALL) rauhmap.inp   $(TESTDIR)/rauh/rootmass/rauhmap.inp       
	$(INSTALL) rauhmap.inp   $(TESTDIR)/rauh/rootlength/rauhmap.inp       
	$(INSTALL) rauhmap.inp   $(TESTDIR)/rauh/arialmass/rauhmap.inp       
	$(INSTALL) rauhall.inp   $(TESTDIR)/rauh/all/rauhall.inp       
	$(INSTALL) rauhratio.inp   $(TESTDIR)/rauh/ratio/rauhratio.inp       
	$(INSTALL) rauhrl.inp   $(TESTDIR)/rauh/rootlength/rauhrl.inp       
	$(INSTALL) rauham.inp   $(TESTDIR)/rauh/arialmass/rauham.inp       
	$(INSTALL) rauhrm.inp   $(TESTDIR)/rauh/rootmass/rauhrm.inp       

mackay: nuzhdinm.inp nuzhdinc.inp vieiram.inp vieirac.inp
	test -d $(TESTDIR)/mackay || mkdir $(TESTDIR)/mackay
	test -d $(TESTDIR)/mackay/nuzhdin || mkdir $(TESTDIR)/mackay/nuzhdin
	$(REMOVE) $(TESTDIR)/mackay/nuzhdin/*
	test -d $(TESTDIR)/mackay/vieira || mkdir $(TESTDIR)/mackay/vieira
	$(REMOVE) $(TESTDIR)/mackay/vieira/*
	$(INSTALL) nuzhdinm.inp   $(TESTDIR)/mackay/nuzhdin/nuzhdinm.inp      
	$(INSTALL) nuzhdinc.inp   $(TESTDIR)/mackay/nuzhdin/nuzhdinc.inp       
	$(INSTALL) vieiram.inp   $(TESTDIR)/mackay/vieira/vieiram.inp       
	$(INSTALL) vieirac.inp   $(TESTDIR)/mackay/vieira/vieirac.inp       

