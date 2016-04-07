
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

LFLAGS += -lgsl -lgslcblas -lhdf5 -lm
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  #CFLAGS +=-O2 -g -c -w1 -openmp -I${GREAD}
  LINK +=${CXX} #-openmp
else
  #CFLAGS +=-O2 -g -c -Wall -fopenmp -I${GREAD}
  LINK +=${CXX} #-openmp $(PRO)
endif

objs =  readfile.o accel.o


all:main

main: 	main.o ${objs}
		${LINK} ${LFLAGS} $^ -o $@

readfile.o: readfile.cpp main.h
accel.o:	accel.cpp main.h
main.o: main.cpp main.h


clean: 
	rm -f main main.o readfile.o

