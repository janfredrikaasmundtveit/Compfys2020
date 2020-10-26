#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>
#include "planet.h"
#include "setup.h"


solarsystem setup(int option,double beta,int rel){
	solarsystem s;
	//initial positions(AU) and velocities (AU/day), also declare force but set to 0.
	coord sunqi(0.0,0.0); coord sunvi(0.0,0.0);	coord sunFi(0.0,0.0); 
	//coord sunqi(-1.550668714070811E-04,7.251508601760199E-03); coord sunvi(-7.588255259323373E-06, 2.581757885123079E-06);	coord sunFi(0.0,0.0); 
	coord mqi(-3.009168510839890E-01,7.465585461322467E-03); coord mvi(1.527008630750197E-02,-1.743504723628711E-02);	coord mfi(0.0,0.0); 
	coord vqi(7.252714351031906E-01,7.465585461322467E-03); coord vvi(-1.009455335419807E-04,2.013651840696968E-02); coord vfi(0.0,0.0);
	coord eqi(9.763619062330592E-01,2.225327099640603E-01); coord evi(-3.988919278934853E-03,1.674541445029773E-02); coord efi(0.0,0.0);
	coord maqi(1.355922774972073,-2.677968997612383E-01); coord mavi(3.307991206973E-03,1.491281614814090E-02); coord mafi(0.0,0.0);
	coord jqi(-2.712032022234562,-4.631713764473182E+00); coord jvi(6.422390077545354E-03,-3.452667922526868E-03); coord jfi(0.0,0.0);
	coord sqi(1.507501514277910,-9.941840797581150); coord svi(5.209269080532190E-03,8.181258887470869E-04); coord sfi(0.0,0.0);
	coord uqi(1.719195399355221E+01,9.971925134722831); coord uvi(-2.002240543902494E-03,3.218857618428033E-03); coord ufi(0.0,0.0);
	coord nqi(2.891397789182760E+01,-7.746949395325852); coord nvi(7.911747489361872E-04,3.050473747377198E-03); coord nfi(0.0,0.0);
	coord pqi(1.161977279419817E+01,-3.157851452292777E+01); coord pvi(3.023873754349508E-03,4.331006906136683E-04); coord pfi(0.0,0.0);

	planet sun(1.0,0.0,sunqi,sunvi,sunFi);
	planet mercury(1.5E-7,0.0,mqi,mvi,mfi);
	planet venus(2.5E-6,0.0,vqi,vvi,vfi);
	planet earth(3E-6,0.0,eqi,evi,efi);
	planet mars(3.3E-7,0.0,maqi,mavi,mafi);
	planet jupiter((1.9/2)*1E-3,0.0,jqi,jvi,jfi);
	planet saturn(2.8E-4,0.0,sqi,svi,sfi);
	planet uranus(4.4E-5,0.0,uqi,uvi,ufi);
	planet neptune(0.5E-4,0.0,nqi,nvi,nfi);
	planet pluto(0.7E-8,0.0,pqi,pvi,pfi); 

	if(option==0){
	 //full solarsytem including pluto
		s.setsize(10);
		s.addplanet(sun,0); s.addplanet(mercury,1); s.addplanet(venus,2); s.addplanet(earth,3); s.addplanet(mars,4);s.addplanet(jupiter,5);
		s.addplanet(saturn,6);s.addplanet(uranus,7);s.addplanet(neptune,8);s.addplanet(pluto,9);
	}

	if(option==1 || option==4 || option==5){
	 //sun-earth system
		s.setsize(2);
		s.addplanet(sun,0); s.addplanet(earth,1);
	}
	if(option==2){
	//sun-earth-jupiter system system
		s.setsize(3);
		s.addplanet(sun,0); s.addplanet(earth,1); s.addplanet(jupiter,2);
	}

	if(option==3){
	//sun-mercury system
		coord rmqi=(0.3075,0);
		coord rmvi=(0,12.44/365);
		planet relmerc(1.5E-7,0.0,rmqi,rmvi,mfi); // starting at perehelion x=0
		s.setsize(2);
		s.addplanet(sun,0); s.addplanet(relmerc,1);
	}

	for (int i = 0; i < s.size; i++)
	{	s.p[i].F=s.totforce(i,beta,rel); //computing intial force
		s.p[i].v=s.p[i].v*365.0;//converting from AU/day to Au/y
		s.p[i].updateE(beta); 	
		s.p[i].updatel();


	}

return s;

}



