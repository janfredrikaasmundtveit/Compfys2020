#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <mpi.h>
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



int main (int argc, char* argv[])
{ string filename;
 double w,alpha,beta;
// set up in master op

 if (argc < 2) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins and MC cycles" << endl;
    exit(1);
  } 
  if (argc > 2) {
    filename=argv[1];
    w=atof(argv[2]);
    
  }


}