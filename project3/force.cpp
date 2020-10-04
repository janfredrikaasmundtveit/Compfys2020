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


coord force(planet p1, planet p2){
double pi=acos(-1.0);
double c=4*pi*pi;
coord F;
coord d;
d.x=p1.q.x-p2.q.x;
d.y=p1.q.y-p2.q.y;
double R=sqrt(d.x*d.x+d.y*d.y); 
F.x=c*p1.m*p2.m*(p1.q.x-p2.q.x)/(R*R*R);
F.y=c*p1.m*p2.m*(p1.q.y-p2.q.y)/(R*R*R);

return F;
}
