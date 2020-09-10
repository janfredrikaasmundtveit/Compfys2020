//project 1
//project 1
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include<armadillo>
// use namespace for output and input
using namespace std;
 using namespace arma;
 ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}


 int main(int argc, char *argv[])
 {
  clock_t tStart = clock();
 	// creating a file containing number of intergrationpoints.
	int n; 
    string filename;
    int exponent;
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
}
		for (int i = 1; i <= exponent; i++){
      int n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);
      double h = 1.0/(n);
      double hh = h*h;
      // Set up arrays a lower off diagonal, b main diagonal, c upper of diagonal, u unknown, g modifid f(x)
     // double *b = new double [n+1]; double *c=new double [n+1]; double *a = new double [n+1]; double *u = new double [n+1];
     double *x = new double[n+1];// double *g =new double [n+1];
    
  
     vec b(n);
     mat A(n,n);
      vec u(n);
      for (int i = 0; i < n; i++)
      {
      b(i) = hh*f(i*h);
for(int j = 0; j < n-1; j++){   
     if (i==j)
      { A(i,j)=2.0;
    //B[i][j]=2.0;         
         }

     else if(i==j+1 || i==j-1){
   A(i,j)=-1.0;
   //B[i][j]=-1.0;
     } 
     else{
      A(i,j)=0.0;
    //  B[i][j]=0.0;
     }
      }
  		
        u=solve(A,b);
      
       	// gaussian eliminaton
	 
}
 	// making file (not makefile) 
  	 
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-u(i))/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << u(i);
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      } 
         //   ofile << setw(15) << setprecision(8) << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
      ofile.close();
	  	delete [] x; //delete [] a; delete [] b; delete [] c;  delete [] u; delete [] g;
 	}
  cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;

 	return 0;
 }
