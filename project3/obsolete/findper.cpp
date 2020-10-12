#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "planet.h"
#include "verlet.h"
#include "force.h"
#include <armadillo>
#include "findper.h"

int findper(planet p1,planet p2){

if(p1.r()<0.3076){// if planet is within 10^-4 AU of perreheilon.
	return 1;
}else{
	return 0;} 

}