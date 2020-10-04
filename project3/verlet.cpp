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
#include "totforce.h"


planet verletstep(planet p, coord prevF,double h){
p.q.x=p.q.x+h*p.v.x-h*h*p.F.x/p.m;
p.q.y=p.q.y+h*p.v.y-h*h*p.F.y/p.m;
p.v.y=p.v.y-(p.F.y+prevF.y)*h/(2*p.m);
p.v.x=p.v.x-(p.F.x+prevF.x)*h/(2*p.m);

return p;
}


