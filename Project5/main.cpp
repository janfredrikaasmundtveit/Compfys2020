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

using namespace arma;
using namespace std;
ofstream ofile;

int main (int argc, char* argv[])
{ string filename;
  int NSpin, MCC, mc;

// set up in master op

 if (argc < 1) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins and MC cycles" << endl;
    exit(1);
  } 
  if (argc > 1) {
    filename=argv[1];
    
  }


}