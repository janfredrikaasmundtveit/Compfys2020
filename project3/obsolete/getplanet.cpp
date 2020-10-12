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

planet getplanet(planet p0,planet p1,planet p2,planet p3,planet p4,planet p5, planet p6,planet p7,planet p8,planet p9, int p){
	//picks out planet p.
	planet pl;
	if(p==0){
		pl=p0;
	}
	if(p==1){
		pl=p1;
	}
	if(p==2){
		pl=p2;
	}
	if(p==3){
		pl=p3;
	}
	if(p==4){
		pl=p4;
	}
	if(p==5){
		pl=p5;
	}
	if(p==6){
		pl=p6;
	}
	if(p==7){
		pl=p7;
	}
	if(p==8){
		pl=p8;
	}
	if(p==9){
		pl=p9;
	}

	return pl;
}