#ifndef qm_H
#define qm_H

class coord{ //3D coordinates
	public:
	mutable double x, y,z;
	coord(double x=0.0,double y=0.0,double z=0.0):x(x),y(y)mz(z){}
	void updatecoord(coord c){this->x=c.x;this->y=c.y;this->z=c.z;}		
	coord& operator= (const coord& c){
   		x = c.x;
   		y = c.y;
   		z = c.z;
  	 return *this;
	}	
	inline coord operator+ (const coord& c){
  		 return coord(x+c.x,y+c.y,,y+c.y);
	}
	inline coord operator- (const coord& c){
  		 return coord(x-c.x,y-c.y,z-c.z);
	}
	inline coord operator* (double s){//scalar multiplication, note scalars can only be multiplied from the right
  		 return coord(s*x,s*y,s*z);
	}
	inline coord operator* (const coord& c){//vwctor multiplication, crossproduct
  		 return coord(y*c.z-c.y*z,z*c.x-x*c.z,x*c.y-y*c.x)
  	}
	double r(){return x*x+y*y+z*z;} //returns vector(dot)vector


};


#endif