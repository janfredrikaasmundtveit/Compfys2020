//project 1
#include <iostream>
#include <new>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
// use namespace for output and input
using namespace std;
#define   ZERO       1.0E-15
 ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double ** AllocateMatrix(int, int);
void DeallocateMatrix(double **, int, int); 
void MatrixInverse(double **, int);
void WriteMatrix(double **, int);
void MatrixMultiplication(double **, double **, int);
void LUDecomposition(double **, int, int *);
void LUBackwardSubstitution(double **, int, int *, double *);






void MatrixInverse(double **A, int n)
{        
  // allocate space in memory
  int *indx;  
  double *column;
  indx = new int[n];
  column  = new double[n];
  double **Y    = AllocateMatrix(n,n);
  // Perform the LU decomposition 
  LUDecomposition(A, n, indx);   // LU decompose  a[][] 
 // cout << "LU decomposed matrix  A:" << endl;
  // find inverse of a[][] by columns 
  for(int j = 0; j < n; j++) {
    // initialize right-side of linear equations 
    for(int i = 0; i < n; i++) column[i] = 0.0;
    column[j] = 1.0;
    LUBackwardSubstitution(A, n, indx, column);
    // save result in y[][] 
    for(int i = 0; i < n; i++) Y[i][j] = column[i];
  } 
  // return the inverse matrix in A[][] 
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) A[i][j] = Y[i][j];
  } 
  DeallocateMatrix(Y, n, n);     // release local memory 
  delete [] column;
  delete []indx;

}



double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}


void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}



void LUDecomposition(double **a, int n, int *indx)
{
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new double [n];
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a[i][j])) > big) big = temp;
      }
      if(big == ZERO) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }               
      vv[i] = 1.0/big;        
   } 

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a[i][j];    
	 for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
	 a[i][j] = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a[i][j];
	 for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
	 a[i][j] = sum;
	 if((dum = vv[i]*fabs(sum)) >= big) {
  	    big = dum;
	    imax = i;
	 }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
	    dum        = a[imax][k];
	    a[imax][k] = a[j][k];
	    a[j][k]    = dum;
	 }
	 vv[imax] = vv[j];         // also interchange scaling factor 
      }
      indx[j] = imax;
      if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;
      if(j < (n - 1)) {                   // divide by pivot element 
         dum = 1.0/a[j][j];
	 for(i=j+1;i < n; i++) a[i][j] *= dum;
      }
   } // end j-loop over columns
   delete [] vv;   // release local memory
}



void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
{
   int        i, ii = -1, ip, j;
   double     sum;

   for(i = 0; i< n; i++) {
      ip    = indx[i];
      sum   = b[ip];
      b[ip] = b[i];
      if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
      else if(sum) ii = i;
      b[i] = sum;
   }
   for(i = n - 1; i >= 0; i--) {
      sum = b[i];
      for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
      b[i] = sum/a[i][i];
   }
}

void WriteMatrix(double ** Matrix, int n){
  for(int i=0;i < n;i++){
    cout << endl;
     for (int j=0 ; j < n;j++){
        printf("  A[%2d][%2d] = %12.4E  ",i, j, Matrix[i][j]);
     }
  }
    cout << endl;
}

// Straightforward matrix-matrix multiplication

void MatrixMultiplication(double ** a, double **b, int n){
  double **c = AllocateMatrix(n, n);
  for(int i=0;i < n;i++){
     for (int j=0 ; j < n;j++){
       double sum = 0.0;
       for (int k = 0; k < n; k++) sum += a[i][k]*b[k][j];
       c[i][j] = sum;
     }
  }
  WriteMatrix(c,n);
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
      //double *u = new double [n+1]; 
    
    double **A = AllocateMatrix(n, n);
  	double *u=new double [n]; double *b=new double [n]; double *x=new double [n];
	double **B = AllocateMatrix(n, n);
 	

  
 
  // Setting general square matrices with random matrix elements
  for(int i = 0; i < n; i++) { 
  x[i]= i*h;
  b[i]=f(i*h);        

    for(int j = 0; j < n-1; j++){   
     if (i==j)
      { A[i][j]=2.0;
		//B[i][j]=2.0;      	 
		     }

     else if(i==j+1 || i==j-1){
	 A[i][j]=-1.0;
	 //B[i][j]=-1.0;
     } 
     else{
     	A[i][j]=0.0;
    // 	B[i][j]=0.0;
     }
    }
  
}
	
  MatrixInverse(A,n); 

  //WriteMatrix(A,n-1);
  //MatrixMultiplication(B,A,n-1);

    

  
  	for (int i = 0; i <n; i++)
  	{
  		
       double sum = 0.0;
       for (int k = 0; k < n; k++) {sum += A[i][k]*b[k];}
       
  		u[i]=sum;
  		//cout << u[i] << endl;
  	}




 	// making file (not makefile) 
  	 
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-u[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << u[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      
       cout << log10(n) << ",";
       cout << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;

           // ofile << setw(15) << setprecision(8) <<  (double)(clock() - tStart)/CLOCKS_PER_SEC;
      ofile.close();

		//	A.print(); 	
	   delete [] u; delete [] b;  delete [] x;
 	}

 	return 0;
 }
