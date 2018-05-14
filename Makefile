
ifeq ($(CC),cc)
  ICC:=$(shell which icc --tty-only 2>&1)
  #Can we find icc?
  ifeq (/icc,$(findstring /icc,${ICC}))
     CC = icc -vec_report0
     CXX = icpc
  else
     GCC:=$(shell which gcc --tty-only 2>&1)
     #Can we find gcc?
     ifeq (/gcc,$(findstring /gcc,${GCC}))
        CC = gcc
        CXX = g++
     endif
  endif
endif

LFLAGS += -lgsl -lgslcblas -lhdf5 -lboost_program_options -lfftw3 -lm
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O3 -g -c -w1 -openmp -xHost -parallel
  LINK +=${CXX} -openmp -parallel
else
  CFLAGS +=-O3 -g -c -Wall -fopenmp
  LINK +=${CXX} #$(PRO)
  LFLAGS += -lm -lgomp
endif
CXXFLAGS +=${CFLAGS}
CFLAGS+= -std=c++11

OBJS =  readfile.o accel.o integrate.o main.o
INCL   = main.h readfile.h


all:main

main: 	$(OBJS) 
		${LINK} $^ -o $@ ${LFLAGS}

%.o: %.cpp $(INCL) ${CXX} ${CFLAGS} $< -o $@

clean: 
	rm -f main $(OBJS)

