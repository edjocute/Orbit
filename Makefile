
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

LFLAGS += -lgsl -lgslcblas -lhdf5 
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  #CFLAGS +=-O3 -g -c -w1
  LINK +=${CXX} #-openmp
else
  #CFLAGS +=-O2 -g -c -Wall -fopenmp -I${GREAD}
  LINK +=${CXX} #-openmp $(PRO)
endif
CXXFLAGS +=${CFLAGS}

OBJS =  readfile.o accel.o
INCL   = main.h point_type.hpp


all:main

main: 	main.o $(OBJS)
		${LINK} ${LFLAGS} $^ -lm -o $@

#readfile.o: readfile.cpp main.h $(CFLAGS)
#accel.o:	accel.cpp main.h $(CFLAGS)
#main.o: main.cpp main.h
%.o: %.cpp $(INCL) $(CXX) ${CFLAGS} $< -o $@


clean: 
	rm -f main *.o

