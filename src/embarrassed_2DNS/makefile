#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

# Compiler options
# For the Intel F95 compiler use -c -w -O when compiling. In RH9 systems, 
# you also need libc-2.2.93 and link with -lc-2.2.93. For the Absoft F90 
# compiler, use -c -YEXT_SFX=__ -YEXT_NAMES=LCS -O to compile and set the 
# variable UNDER to 0. When UNDER is equal to 1 you can use the variable 
# APPEND to add characters to external functions before compiling.
# When EXTEN is equal to 1 you can use the variable F90EXT 
# to change the extension of the F90/95 files.

#
# Intel comp  jxb  
COMP     = mpiifort -c
LINK     = mpiifort

# REMOVE THE FLAG FOR SIZES LESS THAN OR EQUAL TO 2048
FLAGS    = -mcmodel=medium


UNDER    = 0
EXTEN    = 0
APPEND   = _
F90EXT   = f90
LMPI     = -limf  -lmpi
IDIR     =
LDIR     =
FFTWLDIR =  /linkhome/rech/genlpe01/uft62hz/LIBS/fftw/lib
FFTWLDIR2 = /obs/vballester/LIBS/fftw/lib
LFFTW    = -ldrfftw -ldfftw -L$(FFTWLDIR) -L$(FFTWLDIR2)
SHELL    = /bin/bash
MAKE     = /usr/bin/make
PERL     = /usr/bin/perl
ECHO     = echo

BINDIR   = .

default:
	@$(ECHO) "Usage: make <action>"
	@$(ECHO) "Possible actions are: hd2D, mhd2D, hall25, clean, dist"


pseudospec2D_inc.o:
	$(COMP) $(FLAGS) $(IDIR) pseudospec2D_mod.$(F90EXT) -o $(BINDIR)/pseudospec2D_mod.o
	$(COMP) $(FLAGS) $(IDIR) pseudospec2D_inc.$(F90EXT) -o $(BINDIR)/pseudospec2D_inc.o

pseudospec2D_adv.o: 
	$(COMP) $(IDIR) pseudospec2D_mod.$(F90EXT) -o $(BINDIR)/pseudospec2D_mod.o
	$(COMP) $(IDIR) pseudospec2D_adv.$(F90EXT) -o $(BINDIR)/pseudospec2D_adv.o


pseudospec25D_inc.o:
	$(COMP) $(IDIR) pseudospec25D_inc.$(F90EXT) -o $(BINDIR)/pseudospec25D_inc.o

hd2D:
	@$(MAKE) modify
	@$(MAKE) pseudospec2D_inc.o
	$(COMP) $(FLAGS) $(IDIR) hd2D.emb.$(F90EXT) -o $(BINDIR)/hd2D.emb.o
	$(LINK) $(FLAGS) $(BINDIR)/pseudospec2D_mod.o \
          $(BINDIR)/pseudospec2D_inc.o $(BINDIR)/hd2D.emb.o $(LDIR) $(LFFTW) -o $(BINDIR)/hd2D
	@$(MAKE) undo


modify:
	if [ ${UNDER} -eq 1 ]; then \
	$(PERL) -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 ; \
	for item in `cat external`; do \
	  `echo '$(PERL) -i -pe 's/$$item/$${item}${APPEND}/g' *.f90'`; \
	done; fi
	if [ ${EXTEN} -eq 1 ]; then \
        for file in *.f90; do \
          if [ ${UNDER} -eq 0 ]; then \
             cp $$file $${file}.bak; \
          fi; \
          mv $$file $${file%.f90}.$(F90EXT); \
        done; fi

undo:
	rm -f *.f
	if [ -a fftp.f90.bak ]; then \
        for file in *.f90.bak; do \
          mv $$file $${file%.bak}; \
        done; fi

clean:
	rm -f $(BINDIR)/*.o $(BINDIR)/*.d $(BINDIR)/*.s $(BINDIR)/*.int $(BINDIR)/*.inf $(BINDIR)/*.mod $(BINDIR)/*bal.txt $(BINDIR)/work.pc* $(BINDIR)/ifc*
	@$(MAKE) undo

dist:
	rm -f $(BINDIR)/hd2D $(BINDIR)/mhd2D $(BINDIR)/hall25 $(BINDIR)/adv2D
	@$(MAKE) clean

