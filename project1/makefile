# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= mpic++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= relativeerr.exe 
${PROG} :	    relativeerr.o 
			${CPPflags} relativeerr.o   ${LIB} -o ${PROG}

	
relativeerr.o :		 relativeerr.cpp 
		        	${CPPflags} -c  relativeerr.cpp



	