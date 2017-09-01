#include "adreal.h"
#include "cvec3t.h"
#include "advec.h"

typedef adreal<12,1,double> ADScalar;
typedef CVec3T<ADScalar> ADVec3;
typedef CVec3T<double> Vec3;


void instantiateAll(
	       ADVec3 av1, ADVec3 av2, ADVec3 av3,
	       Vec3 v1, Vec3 v2, Vec3 v3,	       
	       ADScalar af1, double f1) {
  
  ADScalar af;
  ADVec3 av;

  av += v1;
  av -= v1;
  av *= f1;
  av /= f1;

  av =          av1  +  v1;
  av =           v1  + av1;
  av =          av1  -  v1;
  av =           v1  - av1;
  av = compMult(av1,    v1);
  av = compMult( v1,   av1);
  av = compDiv( av1,    v1);
  av = compDiv(  v1,   av1);
  av =     max( av1,    v1);
  av =     max(  v1,   av1);
  av =     min( av1,    v1);
  av =     min(  v1,   av1);
  av =   cross( av1,    v1);
  av =   cross(  v1,   av1);
  av = project( av1,    v1);
  av = project(  v1,   av1);

  av = av1* f1;
  av =  v1*af1;
  av =  f1*av1;
  av = af1* v1;
  av = av1/ f1;
  av =  v1/af1;

  af =   dot(av1, v2);
  af =   dot( v1,av2);
  af =  dist(av1, v2);
  af =  dist( v1,av2);
  af = angle(av1, v2);
  af = angle( v1,av2);

  av =    lerp( v1,av2,af1);
  av =    lerp(av1, v2,af1);
  av =    lerp(av1,av2, f1);
  av =    lerp(av1, v2, f1);
  av =    lerp( v1,av2, f1);
  av =    lerp( v1, v2,af1);

  av =    rotate( v1,av2,af1);
  av =    rotate(av1, v2,af1);
  av =    rotate(av1,av2, f1);
  av =    rotate(av1, v2, f1);
  av =    rotate( v1,av2, f1);
  av =    rotate( v1, v2,af1);


  af =   tripleProd( v1,av2,av3);
  af =   tripleProd(av1, v2,av3);
  af =   tripleProd(av1,av2, v3);
  af =   tripleProd(av1, v2, v3);
  af =   tripleProd( v1,av2, v3);
  af =   tripleProd( v1, v2,av3);
}


main() { 
  ADVec3 av1,av2,av3;
  Vec3 v1,v2,v3;

  set_independent( av1, v1,0);
  set_independent( av2, v2,3);
  set_independent( av3, v3,6);
  ADScalar af1;
  double f1 = 1;
  af1.set_independent(1.0,9);

  instantiateAll(av1, av2, av3,v1, v2, v3, af1,  f1);

}
