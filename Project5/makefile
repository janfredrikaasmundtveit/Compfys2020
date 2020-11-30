# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflagsmp= mpic++ -O3
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas

all: project5.exe
project5.exe :	   	main.o 
			${CPPflagsmp} main.o  ${LIB} -o project5.exe

main.o :		main.cpp 
		        	${CPPflagsmp} -c main.cpp

