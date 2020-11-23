# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflagsmp= mpic++ -O3
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas

all: isingnopar.exe ising.exe
ising.exe :	   	main.o 
			${CPPflagsmp} main.o  ${LIB} -o ising.exe

main.o :		main.cpp 
		        	${CPPflagsmp} -c main.cpp

isingnopar.exe:		nopar.o 
			${CPPflags} nopar.o  ${LIB} -o isingnopar.exe

nopar.o :		nopar.cpp 
		        	${CPPflags} -c nopar.cpp