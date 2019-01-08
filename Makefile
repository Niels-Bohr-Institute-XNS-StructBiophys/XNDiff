# GNU Make v3.81 supports else ifeq, v3.80 not !
#
# Intel AVX Optimization (has little influence) and OpenMP
ifeq ($(OLEVEL),DEBUG)
	OPT:=
else
	OPT:=-O$(OLEVEL)
endif

ifeq ($(CC),gcc)
	OPENMP=-fopenmp
	ifeq ($(OOPTIONS),AVX)
		OPT:=$(OPT) -march=corei7-avx -mtune=corei7-avx
	endif
endif
ifeq ($(CC),g++)
	OPENMP=-fopenmp
	ifeq ($(OOPTIONS),AVX)
		OPT:=$(OPT) -march=corei7-avx -mtune=corei7-avx
	endif
endif
ifeq ($(CC),icc)
	OPENMP=-openmp
	ifeq ($(OOPTIONS),-xP)
		OPT:=$(OPT) -xP
	endif
endif
ifeq ($(CC),icpc)
	OPENMP=-openmp
	ifeq ($(OOPTIONS),-xP)
		OPT:=$(OPT) -xP
	endif
endif

# compiler flags for Includes, Warnings, Definitons, Debug, etc: -I, -W, -D
CFLAGS:=-I/usr/local/include/gsl -Wall -Wno-maybe-uninitialized $(OPENMP) -DARMADILLO=$(ARMADILLO)
ifeq ($(OLEVEL),DEBUG)
	CFLAGS:=-g3 $(CFLAGS)
endif

# Math, GSL and Armadillo libraries, Linker
LFLAGS:=-lm -lgsl -lgslcblas
ifeq ($(ARMADILLO),1)
LFLAGS:=$(LFLAGS) -larmadillo
endif
LFLAGS:=$(LFLAGS) -lc

all: XNDiff
XNDiff: XNDiff.o mtrand.o
#	@echo $(CFLAGS)
#	@echo $(LFLAGS)
	$(CC) $(OPT) $(CFLAGS) -o XNDiff XNDiff.o mtrand.o $(LFLAGS)
mtrand.o: mtrand.cpp
	$(CC) $(OPT) $(CFLAGS) -c mtrand.cpp
XNDiff.o: XNDiff.cpp
	$(CC) $(OPT) $(CFLAGS) -c XNDiff.cpp
clean:
	rm -f *.o XNDiff
