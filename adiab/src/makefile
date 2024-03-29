# !! Define Environment variable JUPITER_ARCH to reflect your architecture
# for ex.: setenv JUPITER_ARCH LINUX (or export JUPITER_ARCH=LINUX depending on your shell)
# Possible choices : undefined, LINUX, INTEL, MACOSX, EV68, SPARC, OPTERON

# Generic IA32 Linux platform 
# for pentium 3 or higher (otherwise edit accordingly)
# JUPITER_ARCH must be set to LINUX
CC_LINUX  = gcc
OPT_LINUX = -O4 -march=pentium -Wall -ffast-math
OPTSEQ_LINUX = 
PARAOPT_LINUX =
PARACC_LINUX = mpicc

# JUPITER_ARCH must be set to OPTERON
CC_OPTERON = gcc
OPT_OPTERON = -O3 -m64 -Wall -ffast-math
#OPT_OPTERON = -g
OPTSEQ_OPTERON = 
PARAOPT_OPTERON =
#PARACC_OPTERON = mpicc
PARACC_OPTERON = /opt/mpich/gnu/bin/mpicc

# Mac Intel platform
# JUPITER_ARCH must be set to MacIntel
CC_MacIntel = clang
OPT_MacIntel = -O3 -Wall -arch x86_64
#OPT_MacIntel = -g -ggdb -Wall
OPTSEQ_MacIntel = 
PARAOPT_MacIntel =
PARACC_MacIntel = mpicc

# Intel Linux platform 
# optimized for pentium 4 (otherwise edit accordingly)
# JUPITER_ARCH must be set to INTEL
# (sequential only, MPI identical to generic IA32)
CC_INTEL  = icc
OPT_INTEL = 
OPTSEQ_INTEL = -O3 -tpp7 -ipo -unroll -w2 -Wall -wd1418,310,266,810,1419,981,161 
PARAOPT_INTEL = -O4 -march=pentium3 -Wall -ffast-math
PARACC_INTEL = mpicc

# Macintosh MAC OS X platform (with gcc)
# JUPITER_ARCH must be set to MACOSX
CC_MACOSX  = gcc
#OPT_MACOSX = -Wall -O3 -ffast-math -mcpu=G4 -fstrict-aliasing -malign-natural -fprofile-use
#OPT_MACOSX = -Wall -O3 -ffast-math -mcpu=G4 -fstrict-aliasing -malign-natural
#OPT_MACOSX = -Wall
OPT_MACOSX = -Wall -g
OPTSEQ_MACOSX = 
PARAOPT_MACOSX =
PARACC_MACOSX = mpicc

# Generic options valid for any platfom
# which has 'cc' in the path
# These are used if you leave JUPITER_ARCH undefined
CC_  = cc
#OPT_ = -g
#normal mode
OPT_ = -O3
#OPT_ = -O0 -g #debug mode
OPTSEQ_ = 
PARAOPT_ =
PARACC_ = mpicc
#PARACC_ = /softs/openmpi-1.6.3-mt-intel-13/bin/mpicc

# Setup for Compaq HP/SC (Digital EV68)
# with no OpenMP support
# Set the JUPITER_ARCH variable to EV68
CC_EV68  = cc
OPT_EV68 = -arch ev6 -fast -O3 -inline all -msg_disable ompdirignored,macroredef
OPTSEQ_EV68 = 
PARAOPT_EV68 =
PARACC_EV68 = cc

# Generic Sparc (v7) SunStation platform w/ gcc
# JUPITER_ARCH must be set to SPARC
CC_SPARC  = gcc
#OPT_SPARC = -O3 -Wall
OPT_SPARC = -g
OPTSEQ_SPARC = 
PARAOPT_SPARC =
PARACC_SPARC = mpicc
#
#
#
#--------------------No Changes Needed after this line------------------------
#
#
#
SHELL	     =  /bin/sh

MAINOBJ      =  main.o out.o calc_gamma.o testingGamma.o writetextfile.o secant_method.o in.o low.o var.o switch.o create.o iter.o\
		create_fluid.o actual.o diagnostics.o split.o desc.o\
		Output.o Init.o ExtPot.o flux.o diskdrift.o source.o slope.o\
		faces.o riemann.o plm.o correct.o beam.o gfo.o dumppot.o\
		alloc_beam.o sample.o kernel.o cfl.o stress.o visco.o\
		comm_init.o comm_strip.o comm_alloc.o comm_adapt.o\
		comm_exec.o comm_fill.o comm_check.o comm_flag.o\
	        boundary.o overlap.o monitor.o usage.o fpe.o dumpinit.o\
		gridfile.o gridbuild.o gridread.o gridabs.o\
		readfield.o merge_desc.o merge_field.o divergence.o\
		refine_desc.o refine_field.o keplerian.o definout.o\
	 	fillslope.o muscl.o multifluid.o\
		beambound.o initset.o initcode.o potset.o potcode.o\
		prescor.o Stockholm.o gfo_adiab.o plm_adiab.o redefoptions.o\
		subcycling.o log.o timespent.o paracheck.o globalinit.o\
		ray_tracing.o ray_integ.o ray_background.o\
		hstatprep.o hstatwrite.o hstatinit.o hstatscorr.o\
		hstatvalues.o meridcorrflux.o meridhighlev.o\
		enecorr.o fillsource_diff.o fillsource_pot.o fillsource_geom_col.o\
		filletot.o fillsource_geom_rad.o fillsource_predict.o\
                judit.o RT_accessories.o RT_main.o RT_solve.o RT_matrix.o\
		RT_boundaries.o RT_stellarirrad.o RT_viscousheating.o\
		muscl_adiab.o predictor_adiab.o predictor.o densfloor.o

DUMMY	    = mpi_dummy.o

COMP        = $(CC_$(JUPITER_ARCH))
OPT         = $(OPT_$(JUPITER_ARCH)) -DTEMPFLOOR
OPTSEQ      = $(OPTSEQ_$(JUPITER_ARCH))
PARAOPT     = $(PARAOPT_$(JUPITER_ARCH)) -D_PARALLEL
PARACOMP    = $(PARACC_$(JUPITER_ARCH))
LIBS        = -lm
#PARALIBS    = -L$(MPI_PREFIX)/lib -lmpich
PARALIBS    = -lmpi
AUTOINIT    = libinit.cx initcode.c includes/init.h dumpinit.c
AUTOINCL    = includes/var.h includes/extvar.h includes/global_switches.h
AUTOPOT	    = libpot.cx  potcode.c  includes/pot.h dumppot.c


include	.config

EXENAME        = ../jupiter

ARCHIVE		= $(EXENAME:../%=../.src.%.tar)
ARCHIVECOMP	= $(EXENAME:../%=../.src.%.tar.bz2)


SRC = *.c
INCLUDE = includes/*.h

ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ)
else
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(DUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(AUTOINIT) $(AUTOPOT) $(OBJ) $(EXENAME) $(ARCHIVECOMP)
	@echo "" 
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)."
else
	@echo "This built is SEQUENTIAL."
endif
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "make BUILD=sequential"
else
	@echo "make BUILD=parallel"
endif
	@echo ""

$(EXENAME): $(OBJ)
	$(COMPILER) $(OBJ) $(OPTIONS) -o $(EXENAME) $(LIBRARIES)

.PHONY: conditionalrebuild
ifneq ($(BUILD),$(OLDBUILD))
conditionalrebuild: clean
	@echo "OLDBUILD = $(BUILD)" > .config
	@echo "BUILD = $(BUILD)" >> .config
else
conditionalrebuild:
endif

.oldconfig:
.config:

para:
	@make BUILD=parallel

seq:
	@make BUILD=sequential

$(AUTOINCL) : var.c varsubst makefile includes/switches.h includes/jupiter.h
	@./varsubst var.c

$(AUTOINIT) : libinit.txt iniparser.pl makefile Init.c
	@./iniparser.pl

$(AUTOPOT) : libpot.txt potparser.pl makefile ExtPot.c
	@./potparser.pl

$(OBJ): makefile $(INCLUDE)

.PHONY: clean mrproper package

mrproper:
	rm -f *.o *~ *.s *.il $(AUTOINIT) $(AUTOINCL) $(AUTOPOT)\
	includes/*~ $(EXENAME) ../core.* *.tex *.dvi *.pdf *.ps *.log *.aux *.lint *.cx\
	*.gcno *.gcda $(ARCHIVE) $(ARCHIVECOMP)

clean:
	rm -f *.o *~ *.s *.il *.cx

# The Doxyfile provided here is relatively general and should allow
# you to build your own documentation if you modify the source. You
# should have the required executable files in your path
# (doxygen, dot, latex, epstopdf, dvips, gs). On Mac Os X I made soft
# links in /usr/bin to these executables (either in /sw/bin or
# /Applications/Doxygen.app/Content/Resources).
doc:
	doxygen Doxyfile

package: $(ARCHIVECOMP)

release:
	@echo "Creating archive jupiter.tar.gz for release"
	@cd ../..; tar -c\
	 --exclude jupiter/src/dumpinit.c\
	 --exclude jupiter/src/dumppot.c\
	 --exclude jupiter/src/initcode.c\
	 --exclude jupiter/src/potcode.c\
	 -f jupiter.tar jupiter/src/*.c
	@cd ../..; tar rf jupiter.tar jupiter/src/makefile
	@cd ../..; tar rf jupiter.tar jupiter/src/*.txt
	@cd ../..; tar rf jupiter.tar jupiter/src/*.pl
	@cd ../..; tar rf jupiter.tar jupiter/src/varsubst
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/def.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/*types.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/switches.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/*codes.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/mpi_dummy.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/jupiter.h
	@cd ../..; tar rf jupiter.tar jupiter/src/includes/prototypes.h
	@cd ../..; tar rf jupiter.tar jupiter/in/test/1d/*1d
	@cd ../..; tar rf jupiter.tar jupiter/in/test/1d/bifluid1
	@cd ../..; tar rf jupiter.tar jupiter/in/test/1d/*1dlat
	@cd ../..; tar rf jupiter.tar jupiter/in/test/2d/*2d
	@cd ../..; tar rf jupiter.tar jupiter/idl/*.pro
	@cd ../..; tar rf jupiter.tar jupiter/idl/test/*.pro
	@cd ../..; gzip -9 jupiter.tar

$(ARCHIVECOMP) : $(OBJ)
	@echo "Creating $(ARCHIVECOMP)"
	@tar cf $(ARCHIVE) *.c
	@tar rf $(ARCHIVE) includes/*.h
	@tar rf $(ARCHIVE) makefile
	@tar rf $(ARCHIVE) libinit.txt
	@tar rf $(ARCHIVE) iniparser.pl
	@tar rf $(ARCHIVE) libpot.txt
	@tar rf $(ARCHIVE) potparser.pl
	@tar rf $(ARCHIVE) varsubst
	@bzip2 -9 -f $(ARCHIVE)

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS) -Iincludes

#### Hereafter we run the different setups of the test suite

#VIEW	= cd ../idl/test; echo "showtest, '$@'" | idl
VIEW	= cd ../idl/test; echo "showtest, '$@'" | $(IDL_DIR)/bin/idl

# Issue 'make test' to run all the tests

test: test1d test2d

## 1D tests

sod1d sod1da aw1d aw1da daw1d awsph1dlat vkring1d bifluid1: all
	cd ..; ./jupiter in/test/1d/$@
	$(VIEW)

test1d: sod1d sod1da aw1d aw1da daw1d awsph1dlat vkring1d bifluid1

## 2D tests

tiltawsph2d tiltawsph2da sod2d sod2da vkoffring2d: all
	cd ..; ./jupiter in/test/2d/$@
	$(VIEW)

test2d: tiltawsph2d tiltawsph2da sod2d sod2da vkoffring2d
