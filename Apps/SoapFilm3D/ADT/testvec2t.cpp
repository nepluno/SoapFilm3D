#include <math.h>
#include <float.h>
#include <iostream> 
#include "cvec2t.h"

#ifndef M_PI 
#define M_PI 3.141592653589793238462643
#endif

using namespace std;

typedef CVec2T<double> Vec2;

void testresult( const char* test, bool res) { 
  cout << "testing " << test << "..." << (res? "ok":"failed") << endl;
}

int main() { 

  bool ok; 
  Vec2 v1(1,1), v2(1,1), v3(1,2);	 
  Vec2 v4;
  double l;
  int r; 

  // test the cast to *Real
  const double* d = v1; 

  // test accessors 
  v4.x() = 1.0; v4.y() = 2.0;
  ok = (v4.x() - 1.0 == 0 && v4.y() - 2.0 == 0);
  testresult("accessors x(),y()", ok);

  v4(0) = 2.0; v4(1) = 3.0;
  ok = (v4(0) - 2.0 == 0 && v4(1) - 3.0 == 0);
  testresult("accessor (int)", ok);

  // test assignment,difference, length squared to use in the rest of tests
  v4 = v2;        ok = (v4.x() - v2.x() == 0 && v4.y() - v2.y() ==    0);  testresult("operator=", ok);
  v4 = v2 - v3;   ok = (v4.x()          == 0 && v4.y()          == -1.0);  testresult("operator-", ok);
  l  = lenSq(v2); ok = (l == 2.0);                                         testresult("lenSq",     ok);


  normalize(v1);   ok = lenSq(Vec2(1/sqrt(2.0),1/sqrt(2.0)) - v1)<DBL_EPSILON;     testresult("normalize", ok);
  l  = len   (v2); ok = l - sqrt(2.0) == 0;                                        testresult("len",       ok);
  l  = l2    (v2); ok = l - sqrt(2.0) == 0;                                        testresult("l2",        ok);
  l  = l1    (v2); ok = l - 2.0       == 0;                                        testresult("l1",        ok);
  l  = linfty(v2); ok = l - 1.0       == 0;                                        testresult("linfty",    ok);

  v4 =    -v3;     ok = lenSq(Vec2(-1,-2) - v4) == 0;                              testresult("unary-",    ok);
  v4 = dir(v2);    ok = lenSq(Vec2(1/sqrt(2.0),1/sqrt(2.0)) - v4) == 0;            testresult("dir   ",    ok);

  r  = largestAbsComp (v3);  ok = (r == 1);                                        testresult("largestAbsComp" ,ok);
  r  = smallestAbsComp(v3);  ok = (r == 0);                                        testresult("smallestAbsComp",ok);

  v4  = v2;
  v4 += v2;             ok = lenSq(Vec2(2,2)-v4) == 0.0;                         testresult("operator+=",ok);
  v4 -= v2;             ok = lenSq(Vec2(1,1)-v4) == 0.0;                         testresult("operator-=",ok);
  v4 *= 2.0;            ok = lenSq(Vec2(2,2)-v4) == 0.0;                         testresult("operator*=",ok);
  v4 /= 2.0;            ok = lenSq(Vec2(1,1)-v4) == 0.0;                         testresult("operator/=",ok); 
  v4 = v2 + v3;         ok = lenSq(Vec2(2,3)-v4) == 0.0;                         testresult("operator+" ,ok); 
  v4 = v2 - v3;         ok = lenSq(Vec2(0, -1)-v4) == 0.0;                       testresult("operator-" ,ok); 
  v4 = compMult(v2,v3); ok = lenSq(Vec2(1,  2)-v4) == 0.0;                       testresult("compMult"  ,ok); 
  v4 =  compDiv(v3,v3); ok = lenSq(Vec2(1,  1)-v4) == 0.0;                       testresult("compDiv"   ,ok); 
  v4 = v3*2.0;          ok = lenSq(Vec2(2,  4)-v4) == 0.0;                       testresult("operator* post" ,ok); 
  v4 = 2.0*v3;          ok = lenSq(Vec2(2,  4)-v4) == 0.0;                       testresult("operator* pre"  ,ok); 
  v4 = v3/2.0;          ok = lenSq(Vec2(0.5,1)-v4) == 0.0;                       testresult("operator/" ,ok); 

  v4 = max( Vec2(0,4), v3); ok = lenSq(Vec2(1,4)-v4) == 0;                       testresult("max"       ,ok);
  v4 = min( Vec2(0,4), v3); ok = lenSq(Vec2(0,2)-v4) == 0;                       testresult("min"       ,ok);

  l  = dot (v2,v3);     ok = l == 3;                                               testresult("dot"       ,ok);
  l  = dist(v2,v3);     ok = l == sqrt(10.0);                                      testresult("dist"      ,ok);

  l  = angle  (Vec2(0,1),Vec2(0, 1)); ok = l == 0.0;                           testresult("angle 1"   ,ok);
  l  = angle  (Vec2(1,0),Vec2(0, 1)); ok = (l - M_PI/2.0) < DBL_EPSILON;       testresult("angle 2"   ,ok);
  l  = angle  (Vec2(0,1),Vec2(0,-1)); ok = (l - M_PI)     < DBL_EPSILON;       testresult("angle 3"   ,ok);
  l  = angle  (Vec2(0,1),Vec2(0, 1)); ok = (l - M_PI/4.0) < DBL_EPSILON;       testresult("angle 4"   ,ok);

  v4 = project(Vec2(1,1), Vec2(1,0)); ok = lenSq(Vec2(1,0)-v4) == 0.0;      testresult("project 1" ,ok); 
  v4 = project(Vec2(1,1), Vec2(0,1)); ok = lenSq(Vec2(0,1)-v4) == 0.0;      testresult("project 2" ,ok); 
  v4 = project(Vec2(0,0), Vec2(1,1)); ok = lenSq(Vec2(1/3.0,1/3.0)-v4) == 0.0;  testresult("project 3" ,ok); 

  v4 = lerp   (v2,v3,0.  );               ok = lenSq(Vec2(1,1)-v4) == 0.0;      testresult("lerp 1"    ,ok); 
  v4 = lerp   (v2,v3,1.0 );               ok = lenSq(Vec2(1,2)-v4) == 0.0;      testresult("lerp 2"    ,ok); 
  v4 = lerp   (v2,v3,0.25);               ok = lenSq(Vec2(1,1.25)-v4) == 0.0;testresult("lerp 3"    ,ok); 
}
