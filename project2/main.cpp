#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>

// use namespace for output and input
using namespace arma;
using namespace std;
 
 ofstream ofile;
// Functions used
inline double V1(double x){x*x;
}
inline double V2(double x,double w){x*x*w*w+1.0/x;
}

 double normofdi(mat A,int n){
  //calculates the norm of all ellements of a matrix
  double sum=0.0;
      for(int i=0;i<n;i++){
         for(int j=0;j<n;j++){if(i!=j){
             sum=sum+A(i,j)*A(i,j);}
         }
        
      }


  return sum;
 }
 int maxofd(mat A,int *a,int *b,int n){
     //finds the largest ellement of a 
  double max=0.0;
      for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
          if(fabs(A(i,j))>max && i!=j){
              max=fabs(A(i,j));
                *a=i;*b=j;}
        }  
      }
      return 0;
 }

 void rot(mat A, int a, int b, int n, double *c,double *s){
       double tau,t;
        tau=(A(a,a)-A(b,b))/(2.0*A(a,b)); 
         if ( tau >= 0 ) {
           t = (-tau + sqrt(1.0 + tau*tau));
            } 
          else {
            t =(tau +sqrt(1.0 + tau*tau));
          }
        double cos = 1.0/(sqrt(1.0+t*t));
        *s = cos*t;
        *c =cos;

  return;
 }

 mat jacobirotate(mat A, int a, int b, int n, double c, double s){
      

    double aaa,abb,aia,aib;
    aaa = A(a,a);
    abb = A(b,b);
    A(a,a) = c*c*aaa - 2.0*c*s*A(a,b) + s*s*abb;
    A(b,b) = c*c*abb + 2.0*c*s*A(a,b) + s*s*aaa;        
       for(int i=0;i<n;i++){
         if(i!=a && i!=b){
           aia=A(i,a);
           aib=A(i,b);
           A(i,a)=c*aia-s*aib;
           A(i,b)=c*aib+s*aia;
           A(a,i)=A(i,a);
           A(b,i)=A(i,b);
         }

       }

   A(a,b)=A(a,b)*(c*c-s*s)+c*s*(aaa-abb);
   A(b,a)=A(a,b);
    //A(a,b)=0.0;
    //A(b,a)=0.0;


  return A;
 }
  mat eigenvectors(mat A, mat R, int a, int b, int n,double c, double s){
      
        double ria,rib;

        for(int i=0;i<n;i++){
          ria = R(i,a);
          rib = R(i,b);
          R(i,a) =c*ria - s*rib;
          R(i,b) = c*rib + s*ria;
        }


  return R;
 }

 int main(int argc, char *argv[])
 {
 //	clock_t tStart = clock();
 	// creating a file containing number of intergrationpoints.

    string filename;
    int exponent,option;
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
       exponent = atoi(argv[2]);
       option = atoi(argv[3]);
      
}
 	//clock_t tStart = clock();
      //int n = (int) pow(10.0,i);
      int n=exponent;
      // Declare new file name
     string fileout = filename;
      // Convert the power 10^i to a string
     // string argument = to_string(i);
      // Final filename as filename-i-
      //fileout.append(".dat");
     double h,hh;
       if(option==0){
        
        h = 1.0/(n);
        hh = h*h;}
      if(option==1 || option==2){
        
        h = 10.0/(n); // using infty=10
        hh = h*h;}
    
      dmat A(n, n,fill::zeros);
      dmat R(n, n, fill::zeros);

      int a,b;
      double c,s;
      for(int i=0; i<n; i++){
       if(option==0){
          A(i,i)=2.0;
       }
       if(option==1){
          A(i,i)=2.0+V1(i*h);
       }
       if(option==2){
        A(i,i)=2.0+V2(i*h,5.0); //write omega_r in here, this is not optimal.

       }
        R(i,i)=1.0;
          if(i!=n-1){
            A(i,i+1)=-1.0;
            A(i+1,i)=-1.0;
          }
      }

        mat B=A; //duplicating for tesing

      clock_t tStart = clock();

      double norm=1.0;
      int iter=0;
         while(norm>1E-14 && iter<1E7){
        maxofd(A,&a,&b,n);
        if(A(a,b)!=0){
        rot(A,a,b,n,&c,&s);
        R=eigenvectors(A,R,a,b,n,c,s);
        A=jacobirotate(A,a,b,n,c,s);
        norm=normofdi(A,n);
        iter++;
      }
        }


        double tjac = (double)(clock() - tStart)/CLOCKS_PER_SEC;
        cout << "jacobis metod: ";
        cout << tjac << endl;  
        tStart = clock();
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, B);
        double tarma = (double)(clock() - tStart)/CLOCKS_PER_SEC;
       cout << "armadillo: ";
        cout << tarma << endl;
                 cout << "jacobi/armadillo: ";
        cout << tjac/tarma << endl;

          //coment out for large matrecies
/*
          for(int i=0; i < n; i++){
           int j=0;
           while(fabs(A(i,i)-eigvec(j))<10E-8){ 
                if(j==n){
                  cout << "eigen value:" << A(i,i) << "not within presision." << "\n";
                  break;
                } 
                j++;
           }
          }*/


        //test preformed on 5*5 matrecies 
  //cout << "B" << endl;    
  //cout << B << endl;
  //cout << "A" << endl;
  //cout << A << endl;
  //cout << "RTBR-A" << endl;    
  //cout << R.t()*B*R-A << endl;
  //cout << "RTR" << endl;    
  //cout << R*R.t() << endl;
  //cout << norm << endl;
  cout << "number of transforms";      
  cout << iter << endl;
  	 
        if(option==1){

         ofile.open(fileout);
          ofile << setiosflags(ios::showpoint | ios::uppercase);
          ofile << "       x:  Eigenvalue:  err:  " << endl;
           for (int i = 0; i < n; i++){
              ofile << setw(15) << setprecision(8) << h*i;
              ofile << setw(15) << setprecision(8) << A(i,i);
              ofile << setw(15) << setprecision(8) << log10(fabs(A(i,i)-(3.0+i*4.0))) <<"\n" ; 
           }
          

        }
           
        
      
        if(option==2 || option==0){
            double min=1000.0;
      int imin;
      for (int i = 0; i < n; i++){
       if(A(i,i) < min){
          min=A(i,i);
          imin=i;
      }
        
      }

      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
          ofile << "       x:   eigenvector( excact, numeric if aplicapble)    " << endl;
           for(int j=0;j<n;j++){

              ofile << setw(15) << setprecision(8) << h*j;
              if(option==0){
                ofile << setw(15) << setprecision(8) << sin(h*j*3.14159265359);
                 ofile << setw(15) << setprecision(8) << R(j,imin) << "\n";
              }
              if(option==2){
                ofile << setw(15) << setprecision(8) << R(j,imin)*R(j,imin) << "\n";
              }
         } 
            }
    ofile.close();


 	return 0;
 }