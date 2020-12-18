# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas 

all: alpha.exe beta.exe alpha2.exe opptimizedaandb.exe
alpha.exe :	   	alpha.o 
			${CPPflags} alpha.o  ${LIB} -o alpha.exe

alpha.o :		alpha.cpp 
		        	${CPPflags} -c alpha.cpp

beta.exe :	   	beta.o 
			${CPPflags} beta.o  ${LIB} -o beta.exe

beta.o :		beta.cpp 
		        	${CPPflags} -c beta.cpp
alpha2.exe :	   	alpha2.o 
			${CPPflags} alpha2.o  ${LIB} -o alpha2.exe

alpha2.o :		alpha2.cpp 
		        	${CPPflags} -c alpha2.cpp

opptimizedaandb.exe : 	opptimizedaandb.o 
			${CPPflags} opptimizedaandb.o   ${LIB} -o opptimizedaandb.exe

opptimizedaandb.o :		opptimizedaandb.cpp 
		        	${CPPflags} -c opptimizedaandb.cpp









