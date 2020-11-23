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

using namespace arma;
using namespace std;
ofstream ofile;


void mcsampling(int MCC, double temp,int NSpin, mat &lattice,double *Energy,double *MagneticMoment,double *TotE,double *TotE2,double *TotM,double *TotM2,double *TotabsM, double *accepted){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  int totspin=NSpin*NSpin;
  vec EnergyDifference=zeros<vec>(17); 

   

  for( int de =-8; de <= 8; de+=4){ 
    EnergyDifference(de+8) = exp(-de/temp); // using actual energy not just above groundtate

  }   
  for (int i = 0; i < MCC; i++){
    for(int j=0; j<totspin; j++){
        int ix = (int) (RandomNumberGenerator(gen)*(NSpin));
        int iy = (int) (RandomNumberGenerator(gen)*(NSpin));
    
        int deltaE=2*lattice(ix,iy)*(lattice(ix,(iy+1)%NSpin)+lattice(ix,(iy+(NSpin)-1)%NSpin)+lattice((ix+1)%NSpin,iy)+lattice((ix+(NSpin)-1)%NSpin,iy));
           
        if( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8)){
            lattice(ix,iy) *=-1.0;
            *Energy +=((double)(deltaE));
            *MagneticMoment+=2.0*lattice(ix,iy); 
           // *accepted+=1.0;
          }
         
       }

     
        if(i>1000){
          *TotE += *Energy; *TotE2 += (*Energy)*(*Energy); *TotM += *MagneticMoment; *TotM2 +=(*MagneticMoment)*(*MagneticMoment); *TotabsM += fabs(*MagneticMoment);
         }


            /*
            ofile << setw(15) << setprecision(8) << i;
            ofile << setw(15) << setprecision(8) << *Energy/(NSpin*NSpin); 
            ofile << setw(15) << setprecision(8) << *MagneticMoment/(NSpin*NSpin)<< "\n";//for finding where equlibrium
            */


  }
  
}

void output(double temp, int NSpin, int MCC,double TotE[]){
//lattice size, cycles, temprature, other
 double norm = 1.0/((double) (MCC-1000));  // divided by  number of cycles
 double n=1.0/((double)(MCC)); 
   double AllSpins = 1.0/((double) NSpin*NSpin);
   double E_ExpectationValues = TotE[0]*norm; 
  double E2_ExpectationValues = TotE[1]*norm; 
  double M_ExpectationValues = TotE[2]*norm;
  double M2_ExpectationValues = TotE[3]*norm;
  double Mabs_ExpectationValues = TotE[4]*norm;  

  double HeatCapacity = (E2_ExpectationValues - E_ExpectationValues*E_ExpectationValues)/(temp*temp);
  double MagneticSusceptibility = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/(temp);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temp;
  ofile << setw(15) << setprecision(8) << (E_ExpectationValues)*AllSpins;
  ofile << setw(15) << setprecision(8) << HeatCapacity*AllSpins;
  ofile << setw(15) << setprecision(8) << M_ExpectationValues*AllSpins;
  ofile << setw(15) << setprecision(8) << MagneticSusceptibility*AllSpins;
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues*AllSpins << "\n";



}



int main (int argc, char* argv[])
{ string filename;
  int NSpin, MCC, mc;
  double ITemp, FTemp, TempStep;
 



 if (argc < 3) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins and MC cycles" << endl;
    exit(1);
  } 
  if (argc > 1) {
    filename=argv[1];
    NSpin = atoi(argv[2]);
    mc = atoi(argv[3]);  //mccycles  
    MCC=pow(10,mc);
    ITemp =2.2;
    FTemp = 2.3;
    TempStep = 0.01;
  }

    string fileout = filename;
    ofile.open(fileout);
  


  double  TimeStart, TimeEnd, TotalTime;
  //TimeStart = MPI_Wtime();
     std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
// for(int m=1000; m<MCC;m=10*m){ //loop for comparing accepeted configurations 
  mat lattice=ones<mat>(NSpin,NSpin); 
  for(int i=0;i<NSpin;i++){ //randomize lattice comment out to use groundstate
    for(int j=0;j<NSpin;j++){
        if( RandomNumberGenerator(gen) <=0.5){
          lattice(i,j)=-1;
        }
      }

 }
  double Energy= NSpin*NSpin*(-2); //E 
  double MagneticMoment=NSpin*NSpin; //M
	 
  for (double temp=ITemp; temp < FTemp; temp+=TempStep){ 
    double TotE = 0.0; 
    double TotE2 = 0.0;
    double TotM = 0.0;
    double TotM2 = 0.0;
    double  TotabsM = 0.0;
    double accepted=0.0; //number of accepted configurations 
    //double temp=1.0; //fix temprature coment out if looping over tempratures
   
   /*    mcsampling(m,temp,NSpin,lattice,&Energy,&MagneticMoment,&TotE,&TotE2,&TotM,&TotM2,&TotabsM,&accepted);//use if looping over number of mccycles for comparing accepeted configurations comment out otherwise

       ofile << setw(15) << setprecision(8) << m;
       ofile << setw(15) << setprecision(8) << (accepted)/((double)(m*NSpin*NSpin)) << "\n";*/
    mcsampling(MCC,temp,NSpin,lattice,&Energy,&MagneticMoment,&TotE,&TotE2,&TotM,&TotM2,&TotabsM,&accepted);//use if not looping over number of mccycles
  		
  	double LocE[5]; //E,E^2,M,M^2,|M|
    LocE[0]=TotE;
    LocE[1]=TotE2;
    LocE[2]=TotM;
    LocE[3]=TotM2;
    LocE[4]=TotabsM;


	 
  		output(temp,NSpin,MCC,LocE);	
     


  }

  
 
     ofile.close();  // close output file
  
  //  End MPI


return 0;
}

