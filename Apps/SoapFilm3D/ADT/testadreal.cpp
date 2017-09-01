#include <iostream> 
#include "adreal.h"
#include "cvec3t.h"

using namespace std;

#define DO_HESS_TEST 0

typedef CVec3T<double> vec3;
typedef  adreal<3,DO_HESS_TEST,double>  adouble3;
typedef  CVec3T< adouble3 >  advec3;

  // automatically generated tests
  #include "testadreal_maple.cpp"




int main() { 
  adouble3 v1, v2, v3;
  HessianType<3,double> hh;
  v1.set_independent(0.5,0);
  v2.set_independent(3,1);

  // instantiate all 

  v1 = v2;
  v1 += v2;
  v1 -= v2;
  v1 *= v2;
  v1 /= v2;

  bool b; 
  b = v1 <  v2;
  b = v1 <= v2;
  b = v1 >  v2;
  b = v1 >= v2;
  b = v1 == v2;
  b = v1 != v2;


  b = 1.0 < v2;
  b = 1.0 <= v2;
  b = 1.0 > v2;
  b = 1.0 >= v2;
  b = 1.0 == v2;
  b = 1.0 != v2;

  b = v1 <  1.0;
  b = v1 <= 1.0;
  b = v1 >  1.0;
  b = v1 >= 1.0;
  b = v1 == 1.0;
  b = v1 != 1.0;


  v3 = v1 + v2;
  v3 = v1 - v2;
  v3 = v1*v2;
  v3 = v1/v2;
  v3 = atan2(v1,  v2);
  v3 = pow(v1,  v2);

  v3 = v1 + 1.0; 
  v3 = v1 - 1.0; 
  v3 = v1*1.0; 
  v3 = v1/1.0;
  v3 = atan2(1.0, v1);
  v3 = pow(1.0, v1);


  v3 = 1.0 + v1;
  v3 = 1.0 - v1;
  v3 = 1.0*v1;
  v3 = 1.0/v1;
  v3 = atan2(v1, 1.0);
  v3 = pow(v1, 1.0);

  v2 = -v1;
  v2 = sqrt(v1);
  v2 = cos(v1);
  v2 = sin(v1);
  v2 = tan(v1);
  v2 = exp(v1);

  v2 = asin(v1);
  v2 = acos(v1);
  v2 = atan(v1);
  v2 = log(v1);


  v1 = min(v2,v3);
  v1 = min(v2,1.0);
  v1 = min(1.0,v3);

  v1 = max(v2,v3);
  v1 = max(v2,1.0);
  v1 = max(1.0,v3);


  // instantiate vector ops
  advec3  V1, V2, V3;
  advec3  V4(v1,v2,v3);
  advec3  V5(V1); 
  V1 = V2;
  V1 = V2 + V3;
  V1 = V2 - V3;
  V1 = v1 * V2;
  V1 = V2 * v1;
  V1 = V2 / v1;

  V1 = -V1;


  v2 = dot(V1,V2);
  v1 = lenSq(V1);
  v1 = l2(V1);
  V1 = dir(V2);
  normalize(V1);
  V1 = cross(V2,V3);
  V1 = lerp(V2, V3,v1);


  test1(); 
  test2(); 
  test3();
  test4();
  test5();


}
