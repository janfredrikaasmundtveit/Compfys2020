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


coord force(planet p1, planet p2,double beta,int rel){
	double pi=acos(-1.0);
	double a=4*pi*pi;
	double c=63197.8; //speed of light in AU/y
	coord F;
	coord d;
	d.x=p1.q.x-p2.q.x;
	d.y=p1.q.y-p2.q.y;
	double R=sqrt(d.x*d.x+d.y*d.y); 
	
		F.x=a*p1.m*p2.m*(p1.q.x-p2.q.x)/(pow(R,beta+1));
		F.y=a*p1.m*p2.m*(p1.q.y-p2.q.y)/(pow(R,beta+1));

	if(rel==1){//relatisvistic correction
		p1.l=p1.m*p1.v.x*p1.q.y-p1.m*p1.v.y*p1.q.x;// angular momentum (overall sign is irelevant)
		F.x=F.x*(1.0+3.0*p1.l*p1.l/(R*R*c*c));
		F.y=F.y*(1.0+3.0*p1.l*p1.l/(R*R*c*c));
	}
	return F;
}
