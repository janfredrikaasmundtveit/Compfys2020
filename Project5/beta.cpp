#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <string>
#include <time.h>
#include <random>
#include<armadillo>
#include "qm.h"

using namespace arma;
using namespace std;
ofstream ofile;


void output(double alpha,double beta,int accepted,double E,double E2,int MC){

	  		ofile << setw(15) << setprecision(8) << alpha;
	  		ofile << setw(15) << setprecision(8) << beta;
 		 	ofile << setw(15) << setprecision(8) << ((double)(accepted))/((double)(MC));
 		 	ofile << setw(15) << setprecision(8) <<	E/((double)(MC));
 			ofile << setw(15) << setprecision(8) <<	(E2/((double)(MC)))-(E/((double)(MC)))*(E/((double)(MC))) << "\n";

}

int main (int argc, char* argv[])
{ string filename;
 double w,beta,alpha;
 int MC;
// set up in master op

 if (argc < 2) {
    cout << "Bad Usage: " << argv[0] <<   " read output file, omega and MC cycles" << endl;
    exit(1);
  } 

  if (argc > 2) {
    filename=argv[1]; 
  	w=atof(argv[2]);
 	MC=pow(10,atoi(argv[3]));
 	alpha=atof(argv[4]);
}
       ofile.open(filename);
   ofile << setiosflags(ios::showpoint | ios::uppercase);
	std::random_device rd;
 	std::mt19937_64 gen(rd());
  	std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  

    
  	double sw=sqrt(w); //to avoid a sqrt in loop.
	beta=0.1;
	 
 	while(beta<1.0){ 

 			coord r1(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
  			coord r2(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
  			double delta=(0.704)/(sw*sqrt(alpha));
		  	double E=0.0;
		  	double E0=0.0;
		  	double E2=0.0;
  			int accepted=0;
  			psit2 pt1(alpha,beta,w);
			for(int i=0; i<MC; i++){
					
		  		
		  		coord rand1(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
		  		coord nr1=r1+rand1*delta; 
		  		coord rand2(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
		  		coord nr2=r2+rand2*delta;

		  		if(pt1.psi2(nr1,nr2)/pt1.psi2(r1,r2)>=RandomNumberGenerator(gen)){
		  			//accept move 
		  			r1=nr1;
		  			r2=nr2;
		  			accepted++;
		  			

		  		}	
		  		
		  			double locE=pt1.H(r1,r2);
		  		
		  			E+=locE;
		  			E2+=locE*locE;
		  	
		  	}
		   
  			
		  	output(alpha,beta,accepted,E,E2,MC);
		  	if (beta>0.3 && beta < 0.4) // after inital run we spesifi region we want hiher resoulusion in
		  	{
				beta+=0.001; 
		  	}else{
		  		beta+=0.1;
		  	}

			
 	} 

     ofile.close();
}  