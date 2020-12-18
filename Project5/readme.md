this repository contains the nessesary sourecode for 4 programs used for calculating the groundstate enrgy of 2 electrons in a 3D harmonic oscillator potential.

The file qm.h contans 3 classes. The first class coord defines a 3D vector. 
The class psit1 defines an unnormilized wavefunction of the form e^(-alpha*w*(r1^2+r2^2)), where alpha and omega are defind when creating a wavefunction.
it contains member functions psi2(coord r1,coord r2), H0(coord r1,coord r2) and H(coord r1,coord r2). psi2 takes 2 coords as input and returns |psi|^2(r1,r2).
H0 takes  2 coords as input and returns the Eigenvalue of the noniteracting Hamiltonian H0(r1,r2)= $0.5*(\del_{r1}^2+\del_{r2}^2)+0.5*(\omega* r1^2+\omega* r2}^2).
H takes  2 coords as input and returns the Eigenvalue of the iteracting Hamiltonian H(r1,r2)= $0.5*(\del_{r1}^2+\del_{r2}^2)+0.5*(\omega* r1^2+\omega* r2}^2)+1/sqrt(r_1^2+r_2^2).

The class psit2 defines an unnormilized wavefunction of the form e^(-alpha*w*(r1^2+r2^2))e^(r12/2(1+beta* r12)).
in addition to the memberfunctions of psit1 it contain a function T returning the kinetic energy and V returning the potenial energy. 

The first program alpha.exe whose code is found in alpha.cpp takes 3 inputs first a string deining the name of the output file, second a double which defines the frequecy of the harmonic oscillator potential and finaly an integer defining how many monte carlo cycles  should be used. The program calculeates returns the number of accepted moves, the expectationvalue of the energy and the variance of the energy in the state psit1 as for choosen values of alpha.

The program beta.exe whose code is found in beta.cpp takes in addition to the inputs of alpha.exe a forth input which fixes the value of alpha and returns an output simmilar to alpha.exe except it uses the psit2, and varies beta. The idea is to use the alpha which produces a minimum from alpha.exe and input it here. 

the program alpha2.exe takes the same inputs as beta.exe, however the input referred to as alpha in beta.exe takes the role of beta. The program returns energies for different values of alpha. The idea is to use the value of beta which minnimizeds the energy from beta.exe as an input.

finaly the program opptimizedaandb.exe takes the same inputs as alpha.exe and returns energies as a function of the frequency. The progaram as given here uses psit2 with optimized parametres alpha and beta determinded by runing alpha.exe,beta.exe and alpha2.exe myself. In order to modify this you need to change the code opptimizedaandb.cpp.
In order to get the energies for psit1 with opptimized parameters as a function of the frequency comment out line 68 and comment in line 69.
