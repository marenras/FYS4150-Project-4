# General makefile for c++ - choose PROG =   name of given program

# Here we define compiler option, libraries and the  target
CPPflags = c++ -O3

# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas

# Here we define the name of the executable
PROG = run.x


${PROG} :	  		main.o  ising.o
					${CPPflags} main.o ising.o -o ${PROG} ${LIB}

main.o :			main.cpp
		        	${CPPflags} -c main.cpp

ising.o :			ising.cpp
		        	${CPPflags} -c ising.cpp
