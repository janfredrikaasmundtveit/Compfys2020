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

using namespace std;

ofstream ofile;

void output(planet);

int main(int argc, char *argv[])
{
double pi=acos(-1.0);
double a=365;
planet mercury,venus,mars, earth, jupiter,saturn,uranus,neptune,pluto, sun;
//t=0 at 6/10-18
sun.m=1; sun.q.x=0.0; sun.q.y=0.0; sun.v.x=0.0; sun.v.y=0.0; sun.F.x=0.0; sun.F.y=0.0;
//sun.m=1; sun.q.x=-8.687531967304304E-05; sun.q.y=7.227908128550092E-03; sun.v.x=-7.565235600065487E-06*a; sun.v.y=2.663358189669651E-06*a;sun.F.x=0.0; sun.F.y=0.0;
mercury.m=1.5E-7; mercury.q.x=-3.009168510839890E-01; mercury.q.y=-3.298635039614782E-01; mercury.v.x=1.527008630750197E-02*a; mercury.v.y=-1.743504723628711E-02*a; mercury.F.x=0.0;mercury.F.y=0.0; 
venus.m=2.5E-6;venus.q.x=7.252714351031906E-01; venus.q.y= 7.465585461322467E-03; venus.v.x=-1.009455335419807E-04*a; venus.v.y=2.013651840696968E-02*a; venus.F.x=0.0; venus.F.y=0.0;
earth.m=3E-6; earth.q.x=9.763619062330592E-01; earth.q.y= 2.225327099640603E-01; earth.v.x=-3.988919278934853E-03 *a; earth.v.y=1.674541445029773E-02*a; earth.F.x=0.0;earth.F.y=0.0;
mars.m=3.3E-7; mars.q.x=1.355922774972073; mars.q.y=-2.677968997612383E-01; mars.v.x=3.307991206973E-03*a; mars.v.y=1.491281614814090E-02*a; mars.F.x=0.0; mars.F.y=0.0;
jupiter.m=(1.9/2)*1E-3; jupiter.q.x=-2.712032022234562; jupiter.q.y=-4.631713764473182E+00; jupiter.v.x=6.422390077545354E-03*a; jupiter.v.y=-3.452667922526868E-03*a;jupiter.F.x=0.0; jupiter.F.y=0.0;
saturn.m=2.8E-4; saturn.q.x=1.507501514277910; saturn.q.y=-9.941840797581150; saturn.v.x=5.209269080532190E-03*a; saturn.v.y= 8.181258887470869E-04*a; saturn.F.x=0.0; saturn.F.y=0.0;
uranus.m=4.4E-5; uranus.q.x=1.719195399355221E+01; uranus.q.y=9.971925134722831; uranus.v.x=-2.002240543902494E-03*a; uranus.v.y=3.218857618428033E-03*a;uranus.F.x=0.0; uranus.F.y=0.0;
neptune.m=0.5E-4; neptune.q.x=2.891397789182760E+01; neptune.q.y=-7.746949395325852; neptune.v.x=7.911747489361872E-04*a; neptune.v.y=3.050473747377198E-03*a; neptune.F.x=0.0; neptune.F.y=0.0;
pluto.m=0.7E-8; pluto.q.x=1.161977279419817E+01;pluto.q.y=-3.157851452292777E+01; pluto.v.x=3.023873754349508E-03*a; pluto.v.y=4.331006906136683E-04*a; pluto.F.x=0.0;pluto.F.y=0.0;



string filename;

if( argc <= 3 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file, number of steps, final time on same line" << endl;
    exit(1);
  }
  else{
    filename=argv[1];
   }
  
int n=atoi(argv[2]);
double tfinal=atof(argv[3]);
double t=0.0;
double h=tfinal/n;
//coord prevF0=totforce(sun,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,pluto);
coord prevF1=totforce(mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,pluto,sun);
coord prevF2=totforce(venus,mercury,earth,mars,jupiter,saturn,uranus,neptune,pluto,sun);
coord prevF3=totforce(earth,mercury,venus,mars,jupiter,saturn,uranus,neptune,pluto,sun);
coord prevF4=totforce(mars,mercury,venus,earth,jupiter,saturn,uranus,neptune,pluto,sun);
coord prevF5=totforce(jupiter,mercury,venus,earth,mars,saturn,uranus,neptune,pluto,sun);
coord prevF6=totforce(saturn,mercury,venus,earth,mars,jupiter,uranus,neptune,pluto,sun);
coord prevF7=totforce(uranus,mercury,venus,earth,mars,jupiter,saturn,neptune,pluto,sun);
coord prevF8=totforce(neptune,mercury,venus,earth,mars,jupiter,saturn,uranus,pluto,sun);
coord prevF9=totforce(pluto,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,sun);

ofile.open(filename);
ofile << setiosflags(ios::showpoint | ios::uppercase);
		 ofile << setw(15) << setprecision(8) <<t;
		 output(sun);
		 output(mercury);
		 output(venus);
		 output(earth);
		 output(mars);
	     output(jupiter);
	     output(saturn);
	     output(uranus);
	     output(neptune);
	     output(pluto);
	     ofile  << endl;
while(t<=tfinal){
//sun.F=totforce(sun,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,pluto);
mercury.F=totforce(mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,pluto,sun);
venus.F=totforce(venus,mercury,earth,mars,jupiter,saturn,uranus,neptune,pluto,sun);
earth.F=totforce(earth,mercury,venus,mars,jupiter,saturn,uranus,neptune,pluto,sun);
mars.F=totforce(mars,mercury,venus,earth,jupiter,saturn,uranus,neptune,pluto,sun);
jupiter.F=totforce(jupiter,mercury,venus,earth,mars,saturn,uranus,neptune,pluto,sun);
saturn.F=totforce(saturn,mercury,venus,earth,mars,jupiter,uranus,neptune,pluto,sun);
uranus.F=totforce(uranus,mercury,venus,earth,mars,jupiter,saturn,neptune,pluto,sun);
neptune.F=totforce(neptune,mercury,venus,earth,mars,jupiter,saturn,uranus,pluto,sun);
pluto.F=totforce(pluto,mercury,venus,earth,mars,jupiter,saturn,uranus,neptune,sun);

//sun=verletstep(sun,prevF0,h);
mercury=verletstep(mercury,prevF1,h);
venus=verletstep(venus,prevF2,h);
earth=verletstep(earth,prevF3,h);
mars=verletstep(mars,prevF4,h);
jupiter=verletstep(jupiter,prevF5,h);
saturn=verletstep(saturn,prevF6,h);
uranus=verletstep(uranus,prevF7,h);
neptune=verletstep(neptune,prevF8,h);
pluto=verletstep(pluto,prevF9,h);

//prevF0=sun.F;
prevF1=mercury.F;prevF2=venus.F;prevF3=earth.F; prevF4=mars.F; prevF5=jupiter.F;prevF6=saturn.F; prevF7=uranus.F;prevF8=neptune.F; prevF9=pluto.F;
t=t+h;  //c1 t:c2 sun x: c3 sun y: c4 mercury x: c5 mercury y: c6 venus x: c7 venus y: c8  earth x: c9 earth y
		//c10 mars x: c11 mars y: c12 jupiter x: c13 jupiter y: c14 saturn x: c15 saturn y: c16 uranus x: c17 uranus y
		//c18 neptune x: c19 neptune y: c20 pluto x: c21 pluto;
		 ofile << setw(15) << setprecision(8) <<t;
		 output(sun);
		 output(mercury);
		 output(venus);
		 output(earth);
		 output(mars);
	     output(jupiter);
	     output(saturn);
	     output(uranus);
	     output(neptune);
	     output(pluto);
	     ofile  << endl;

}

         ofile.close();
	return 0;
}



void output(planet p){

  	ofile << setw(15) << setprecision(8) <<p.q.x;
	ofile << setw(15) << setprecision(8) <<p.q.y; 

}