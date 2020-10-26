#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>
#include "planet.h"
#include "verlet.h"
#include "force.h"
#include "totforce.h"
#include "findper.h"
#include "getplanet.h"


coord totforce(planet p, planet p1, planet p2,planet p3,planet p4,planet p5, planet p6, planet p7,planet p8,planet p9,double beta,int rel){
	coord F,F1,F2,F3,F4,F5,F6,F7,F8,F9;
	F1=force(p,p1,beta,rel);
	F2=force(p,p2,beta,rel);
	F3=force(p,p3,beta,rel);
	F4=force(p,p4,beta,rel);
	F5=force(p,p5,beta,rel);
	F6=force(p,p6,beta,rel);
	F7=force(p,p7,beta,rel);
	F8=force(p,p8,beta,rel);
	F9=force(p,p9,beta,rel);
	F.x=F1.x+F2.x+F3.x+F4.x+F5.x+F6.x+F7.x+F8.x+F9.x; 
	F.y=F1.y+F2.y+F3.y+F4.y+F5.y+F6.y+F7.y+F8.y+F9.y; 
	return F;
}