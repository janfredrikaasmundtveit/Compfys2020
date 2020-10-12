#ifndef planet_H
#define planet_H

class coord{
	public:
	mutable double x, y;
	coord(double x=0.0,double y=0.0):x(x),y(y){}
	void updatecoord(coord c){this->x=c.x;this->y=c.y;}		
	coord& operator= (const coord& c){
   		x = c.x;
   		y = c.y;
  	 return *this;
	}	
	inline coord operator+ (const coord& c){
  		 return coord(x+c.x,y+c.y);
	}
	inline coord operator- (const coord& c){
  		 return coord(x-c.x,y-c.y);
	}
	inline coord operator* (double s){
  		 return coord(s*x,s*y);
	}
	double r(){return x*x+y*y;} //returns vector(dot)vector


};

class planet{
	public:
	double m,l,E;// mass angular momentum and total energy. 
	coord q,v,F;
	planet(){m=0.0;l=0.0;
		coord temp(0.0,0.0);q=temp;v=temp;F=temp;}
	planet(double m,double l,coord q,coord v,coord F):m(m),l(l),q(q),v(v),F(F) {}
	planet& operator= (const planet& p){
   		m = p.m;
   		l = p.l;
   		q = p.q;
   		v = p.v;
   		F=p.F;
  	 return *this;
	}	
	void updateE(double beta){
		double pi=acos(-1.0);
		double G=4*pi*pi; //newtons constant in au/y
		E=0.5*m*v.r()+G/pow(sqrt(q.r()),beta-1.0);//kinetic+potential energy. here only the gravitional potential of the sun is included. M_s=1.
	}
	void updatel(){
		l=fabs(m*v.y*q.x-m*v.x*q.y); //sign will not be important 
	}
	
	int findper(){
		if(q.r()<0.3076*0.3076){// if planet is within 10^-4 AU of perreheilon spesificaly for mercury. comparing distances^2 to avoid sqrts.
			return 1;
		}else{
			return 0;} 
	}
};

class solarsystem{
public:
	int size; //number of planets including the sun
	mutable planet p[10];
	void setsize(int i){size=i;}
	void addplanet(planet pl,int j){this->p[j]=pl;}

	coord totforce(int i,double beta,int rel){// computethe force acting on planet number i,
 		coord ftemp(0.0,0.0);
 		coord fnew(0.0,0.0);
 		for(int j = 0; j < size; j++){
 			if(i!=j){
 			ftemp=force(p[i],p[j],beta,rel);
 			fnew=fnew+ftemp;
 			}
 		}
 		return fnew;
 	}

	coord force(planet p1, planet p2,double beta,int rel){//computes the force between from planet 2 on planet 1. beta is the exponent of the 1/r, empirical value 2.
		double pi=acos(-1.0);
		double G=4*pi*pi; //newtons constant in au/y
		double c=63197.8; //speed of light in AU/y
		coord d;
		d=p2.q-p1.q;
		double R=sqrt(d.x*d.x+d.y*d.y); 
		double s=G*p1.m*p2.m/(pow(R,beta+1));	
		coord F=d*s;

		if(rel==1){//relatisvistic correction
			F=F*(1.0+3.0*p1.l*p1.l/(R*R*c*c));
		}
		return F;
}


};

 class solver{
 public:
 	solarsystem s;
 	int method;//chooses whether to use euler(1) or verlet(0) method
	solver(int method,solarsystem s):method(method),s(s){}


 	void step(double h,double beta, int rel){
 		coord ftemp; 
	 
 	if(method==0){//verlet
	 	coord prevF;	
	 
			for(int i=1;i<s.size;i++){//incud 0 in loop if sun shoud be dynamic
				double minv=1.0/s.p[i].m;
				prevF=s.p[i].F;
				s.p[i].q=s.p[i].q+s.p[i].v*h+s.p[i].F*minv*0.5*h*h;
				s.p[i].F=s.totforce(i,beta,rel);
				s.p[i].v=s.p[i].v+(s.p[i].F+prevF)*h*0.5*minv;
				s.p[i].updateE(beta);
				s.p[i].updatel();

			}
		}
		if(method==1){//euler
			for(int i=1;i<s.size;i++){
				double minv=1.0/s.p[i].m;
				s.p[i].q=s.p[i].q+s.p[i].v*h;
				s.p[i].v=s.p[i].v+(s.p[i].F)*h*minv;
				s.p[i].F=s.totforce(i,beta,rel);
				s.p[i].updateE(beta);
				s.p[i].updatel();
			}
		}
	}
		
 };



#endif