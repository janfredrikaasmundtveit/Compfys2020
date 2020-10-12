#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>
#include "planet.h"
#include "setup.h"

using namespace std;

ofstream ofile;
  

void out(solver solve,double t){
	ofile << setw(15) << setprecision(8) << t;
	for(int i=0;i<solve.s.size;i++){
  		ofile << setw(15) << setprecision(8) << solve.s.p[i].q.x;
		ofile << setw(15) << setprecision(8) << solve.s.p[i].q.y;
	}
	ofile << setw(15) << setprecision(8) << "\n";
}

void outputper(double tper,double xper,double yper,double theta,double v,double rmin){//writes the time at perehelion, xpos at perh, ypos, angle to positive x. axis velosity using previous perheilon only and distance at perh.
	ofile << setw(15) << setprecision(8) << tper;
	ofile << setw(15) << setprecision(8) << xper;
	ofile << setw(15) << setprecision(8) << yper;
	ofile << setw(15) << setprecision(8) << theta;
	ofile << setw(15) << setprecision(8) << v;
	ofile << setw(15) << setprecision(8) << sqrt(rmin);

}
 
int main(int argc, char *argv[]){  
	double beta=2.0;//exponetent of r in gravitational force
	int rel=0; //set to 1 to include relativistic correction.
	double pi=acos(-1.0);
	double G=4*pi*pi; //newtons constant in au/y
	string filename;
	clock_t tStart = clock();
 
	if( argc <= 4 ){ 
	    cout << "Bad Usage: " << argv[0] <<
	      " read also output file, number of steps, final time on same line" << endl;
	    exit(1);
	  }
	   filename=argv[1];
	    int exponetent=atoi(argv[2]); //number of timesteps
		double tfinal=atof(argv[3]); //time in years
		int option = atoi(argv[4]); //choose system 0 for full solarsystem, system 1 for sun-earth sytem, system 2 for sun-earth-jupiter system
		//, system 3 for sun mercury with GR correction, system 4 to vary beta and output E,l; option 5 to study escape velocities.
	  	int method = atoi(argv[5]); //choose method 0 for verlet, 1 for euler 
	//solarsystem s;
	double n=pow(10,exponetent);
	
	if(option==3){
		rel=1;
	}
	solarsystem s=setup(option,beta,rel);// contains all the nasty numbers and stuff
	solver solve(method,s);	

	double t=0.0;
	double h=tfinal/(n+1);
	int per=0;
	int perprev=0;
	int c=0;
	double rmax=0.0;
	double rmin=10.0;
	double theta,xper,yper,tper,v,tprev,thetaprev;

	//solve.s.p[1].q.x=0.0; solve.s.p[1].q.y=1.0; solve.s.p[1].v.x=2.0*pi; solve.s.p[1].v.y=0.0; //fix circular orbit
	
	double relerrr=0.0;
	ofile.open(filename);
	ofile << setiosflags(ios::showpoint | ios::uppercase);

	if(option!=5 && option!=4){
		while(t<=tfinal){
			solve.step(h,beta,rel);
			if(option!=3 && option!=4 && c==1000){
				out(solve,t); 
				c=0; //only include every 1000th data point to reduce filesize
			} 
			if(option==3){
				per=solve.s.p[1].findper();
				if(per==1){
					perprev=1;
					if(solve.s.p[1].q.r()<rmin){// comparing disctances^2
					 	rmin=solve.s.p[1].q.r();
					 	xper=solve.s.p[1].q.x;
					 	yper=solve.s.p[1].q.y;
					 	tper=t; 
					}
					
				}
				if(per==0 && perprev==1){//first datapoint not within  10^-4 of perrehelion. resetting to not compare prevous perhelion.
					 theta=atan(yper/xper);
					 v=rmin*(theta-thetaprev)/(t-tprev);
					 outputper(tper,xper,yper,theta,v,rmin);	
					 perprev=0;
					  rmin=10.0;
					  tprev=tper;
					  thetaprev=theta;

				}
			}
			
			t=t+h;
			c++;
		}
	}
	double a = 1.0;
	if(option==5){
		coord inpos=solve.s.p[1].q;
		while(solve.s.p[1].v.r()<4.0*G){// while the velocity being checked for escape < 2*theoretical prediction
			t=0;
			a=a+0.01;
			solve.s.p[1].v=solve.s.p[1].v*a;
			solve.s.p[1].q=inpos;// resetting the posision of the planet
			while(t<=tfinal){
				solve.step(h,beta,rel);
				t=t+h;
				if(sqrt(solve.s.p[1].q.r())>1.2){//defing escacpe by reacing distnace 1.1 AU within tfinal years, 
					cout << "escapevellocity=" << sqrt(solve.s.p[1].v.r())  << endl;
					cout << "theoretical value=" << sqrt(2.0*G)  << endl;
					return 0;
				}
			}
		}
	}
	if(option==4){
	s.p[1].q.x=0.0; s.p[1].q.y=1.0; s.p[1].v.x=5.0; s.p[1].v.y=0.0; //change velocity to 5 AU/y, note as the solarsystem is allready provided to the solver this does not afect the solver
		
		while(beta<=3){	
			t=0;
			double relerrE=0.0;
			double relerrl=0.0;

			s.p[1].updateE(beta); 	
			s.p[1].updatel();
			double lin=s.p[1].l; 
			double Ein=s.p[1].E;
			solve.s=s; // resetting the solarsytem 
			 while(t<=tfinal){
				solve.step(h,beta,rel);
				relerrE=relerrE+(fabs((solve.s.p[1].E-Ein)/Ein))/n;
  				relerrl=relerrl+(fabs((solve.s.p[1].l-lin)/lin))/n;
  				t=t+h;
			 }
			ofile << setw(15) << setprecision(8) << beta;
			ofile << setw(15) << setprecision(8) << relerrE;
			ofile << setw(15) << setprecision(8) << relerrl;
			ofile << setw(15) << setprecision(8) << "\n";
			
			beta=beta+0.01;	
			}
		}

	//cout << "avrage relatve deviaton from energy conservation: " << relerrE << "\n";
	//cout << "avrage relatve deviaton from angular momentum conservation: " << relerrl << "\n";
	//cout << "avrage relatve deviaton from circular orbit: " << relerrr << "\n";
	cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    ofile.close();
	return 0;
}


