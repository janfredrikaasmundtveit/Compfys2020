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


void output(double w,int accepted,double d,double E,double E2,double TV,int MC){

			ofile << setw(15) << setprecision(8) << w;
 		 	ofile << setw(15) << setprecision(8) << ((double)(accepted))/((double)(MC));
 		 	ofile << setw(15) << setprecision(8) <<	d/((double)(MC));
 		 	ofile << setw(15) << setprecision(8) <<	E/((double)(MC));
 			ofile << setw(15) << setprecision(8) <<	(E2/((double)(MC)))-(E/((double)(MC)))*(E/((double)(MC)));
 			ofile << setw(15) << setprecision(8) <<	TV/((double)(MC)) << "\n";

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
 
}
       ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
	std::random_device rd;
 	std::mt19937_64 gen(rd());
  	std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  

    
  	 //to avoid a sqrt in loop.
	alpha=0.713;
	beta=0.392;
	double TV=0.0;
	double dtot=0.0; 
 	while(w<=1.0){

 			coord r1(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
  			coord r2(2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0,2.0*RandomNumberGenerator(gen)-1.0);
  			double delta=(0.704)/(sqrt(w*alpha));
		  	double E=0.0;
		  	double E0=0.0;
		  	double E2=0.0;
  			int accepted=0;
  			psit2 pt1(alpha,beta,w); 
  			//psit1 pt1(0.61,w);
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
		  			coord d=r1-r2;
		  			double locE=pt1.H(r1,r2);
		  			dtot+=sqrt(d.r());
		  			TV+=pt1.T(r1,r2)/pt1.V(r1,r2);
		  			E+=locE;
		  			E2+=locE*locE;
		  	
		  	}
		   
  		
		  	output(w,accepted,dtot,E,E2,TV,MC);
		  	w+=0.01;
		  
		 }
			
 	 

     ofile.close();
}  