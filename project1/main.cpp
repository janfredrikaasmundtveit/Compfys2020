
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
// use namespace for output and input
using namespace std;
 
 ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x);}

 int project_1_b(double *b, double *a, double *c, double *g,double *u,int n)
 {

	
      for (int i = 0; i <= n; i++){
        b[i]=2; // defining the matix (this is stupid)
        a[i]=-1;
        c[i]=-1;
      }
       	// gaussian eliminaton
	 for (int i = 1; i < n; i++) {
		b[i]=b[i]-((a[i]*c[i-1])/b[i-1]); // back sub
		g[i]=g[i]-((a[i]*g[i-1])/b[i-1]); //6n flops
      } 

       u[n-1] = g[n-1]/b[n-1];
      for (int i = n-2; i > 0; i--){
        u[i] = (g[i]-(c[i]*u[i+1]))/b[i]; // forward sub //3n flops
}
 	
 	;
 	return 0;
 }

 int project_1_c(double *b,double *g,double *u,int n)
 {
 
       	// gaussian eliminaton
	 for (int i = 1; i < n; i++) {
    double j=1.0*i;
		b[i]=(j+1.0)/j; // back sub
		g[i]=g[i]+(g[i-1]*(j-1.0)/j); //4 flops
	}
       u[n-1] = g[n-1]*n/(b[n-1]);

   


      for (int i = n-1; i > 0; i--){
       double j=1.0*i;
        u[i-1] = (g[i-1]+u[i])*(j-1.0)/j; // forward sub 2flops
}
 	return 0;
 }


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

 	clock_t tStart = clock();
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
      double *b = new double [n+1]; double *c=new double [n+1]; double *a = new double [n+1]; double *u = new double [n+1];
      double *x = new double[n+1]; double *g =new double [n+1];
    
      u[0] = u[n] = 0.0; // boundary conditions
      b[0] = b[n] = 2; // fixed elements of diagonal
     
      for (int i = 0; i <= n; i++){
		x[i]= i*h;
        g[i] = hh*f(i*h); //1*n flop
    }

    //project_1_b(b,a,c,g,u,n);
    project_1_c(b,g,u,n);
 	// making file (not makefile) 
  	 
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 0; i < n+2;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-u[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << u[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
       cout << log10(n) << ",";
       cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    //  ofile << setw(15) << setprecision(8) <<  (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
      ofile.close();
	  	delete [] a; delete [] b; delete [] c;  delete [] u; delete [] x; delete [] g;
 	}
 	 cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
 	return 0;
 }