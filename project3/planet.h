#ifndef planet_H
#define planet_H

class coord{
	public:
	double x, y;	
};

class planet{
public:
	double m;
	coord q,v,F; 
	double r() {return sqrt(q.x*q.x+q.y*q.y);}
};

#endif