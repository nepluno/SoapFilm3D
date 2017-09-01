#include <math.h>
#include <iostream> 
#include "adreal.h"
#include "cvec3t.h"


using namespace std;
typedef double baseclass;
typedef CVec3T<baseclass> Vec3;

//typedef adreal<3,1,double> baseclass;
//typedef CVec3T<baseclass> Vec3;



void instantiateAll(
		    Vec3 c0, Vec3 c1, Vec3 c2, Vec3 c3, 
		    const Vec3& cconst,
		    baseclass l, 
		    baseclass s,
		    baseclass d1, baseclass d2, baseclass d3,
		    baseclass  da[3]
		    ) { 

  bool b; 
  int i; 
  c0     = Vec3(c1);
  c0     = Vec3(d1,d2,d3);
  c0     = Vec3(da);
  const baseclass* dp = c1;
  c0     = c2;
  c0.x() = d1; 
  c0.y() = d2;
  c0.z() = d3;
  c0(0)  = d1; 
  c0(1)  = d2;
  c0(2)  = d3;  
  l  =  cconst.x();
  l  =  cconst.y();
  l  =  cconst.z();
  l  =  cconst(0);
  l  =  cconst(1);
  l  =  cconst(2);
        normalize(c1);
  l  =  lenSq    (c1);
  l  =  len      (c1);
  l  =  l2       (c1);
  l  =  l1       (c1);
  l  =  linfty   (c1);
  c0 =  -c1;
  c0 =  dir(c1); 
  // i  =  largestAbsComp (c1);
  // i  =  smallestAbsComp(c1);
  c0 += c1; 
  c0 -= c1;
  c0 *= s;
  c0 /= s;
  c0 =             c1+ c2;
  c0 =             c1+ c2;
  c0 = compMult   (c1, c2);
  c0 = compDiv    (c1, c2);
  c0 =             c1* s;
  c0 =              s* c1;
  c0 =             c1/ s;
  c0 = max        (c1,c2);
  c0 = min        (c1,c2);
  l  = dot        (c1,c2);
  l  = dist       (c1,c2);
  l  = angle      (c1,c2);
  c0 = cross      (c1,c2);
  c0 = project    (c1,c2);
  //  b  = isCollinear(c1,c2);
  c0 = lerp       (c1,c2,s);
  c0 = lerp       (c1,c2,s);
  c0 = rotate     (c1,c2,s);
  l =  tripleProd (c1,c2,c3);
}

void initializeInstantiate() { 
  Vec3 c0, c1(1,1,1),c2(2,2,2),c3(3,3,3);
  const Vec3 cconst(1,1,1);
  baseclass l, s = 2.0;
  baseclass d1 = 1.0,d2 = 2.0,d3 = 3.0;
  baseclass da[3] = { 1.0, 2.0, 3.0}; 
  bool b;
  int i; 

  instantiateAll(c0,c1,c2,c3,cconst,
		 l,s, d1,d2,d3, da);


}

int main() { 
  initializeInstantiate();
}
