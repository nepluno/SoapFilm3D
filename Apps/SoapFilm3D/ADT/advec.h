#ifndef __ADVEC_H__
#define __ADVEC_H__

#include "cvec3t.h"

// the lines in this file may be up to 160 characters long

#define AD_TEMPLATE template <int NUM_VARS, int DO_HESS, class constreal>
#define ADS  adreal<NUM_VARS,DO_HESS,constreal>
#define ADV CVec3T<ADS > 
#define CV CVec3T<constreal>
#define CS constreal


AD_TEMPLATE void set_independent( ADV& adv, const CV& v, int i) {
  adv(0).set_independent(v(0),i); 
  adv(1).set_independent(v(1),i+1); 
  adv(2).set_independent(v(2),i+2);
}

// assignments const-variable spcializations, only assignments variable =
// const need to be handled; discarding variable information is best done
// explicitly through value()  so const = variable assignments are not allowed

AD_TEMPLATE ADV operator+=(ADV& c1,  const  CV& c2) { return operator+=<ADS,ADS,constreal>(c1,c2);}
AD_TEMPLATE ADV operator-=(ADV& c1,  const  CV& c2) { return operator-=<ADS,ADS,constreal>(c1,c2);}
AD_TEMPLATE ADV operator*=(ADV& c1,  const  CS&  s) { return operator*=<ADS,ADS,constreal>(c1, s);}
AD_TEMPLATE ADV operator/=(ADV& c1,  const  CS&  s) { return operator/=<ADS,ADS,constreal>(c1, s);}

// bivariate const-variable specializations, two per function

AD_TEMPLATE ADV operator+(const ADV& c1,  const   CV&  c2) { return    opplus <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV operator+(const  CV& c1,  const  ADV&  c2) { return    opplus <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV operator-(const ADV& c1,  const   CV&  c2) { return operator- <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV operator-(const  CV& c1,  const  ADV&  c2) { return operator- <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV  compMult(const ADV& c1,  const   CV&  c2) { return  compMult <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV  compMult(const  CV& c1,  const  ADV&  c2) { return  compMult <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV   compDiv(const ADV& c1,  const   CV&  c2) { return   compDiv <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV   compDiv(const  CV& c1,  const  ADV&  c2) { return   compDiv <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV       max(const ADV& c1,  const   CV&  c2) { return       max <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV       max(const  CV& c1,  const  ADV&  c2) { return       max <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV       min(const ADV& c1,  const   CV&  c2) { return       min <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV       min(const  CV& c1,  const  ADV&  c2) { return       min <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV     cross(const ADV& c1,  const   CV&  c2) { return     cross <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV     cross(const  CV& c1,  const  ADV&  c2) { return     cross <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADV   project(const ADV& c1,  const   CV&  c2) { return   project <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADV   project(const  CV& c1,  const  ADV&  c2) { return   project <ADS,constreal,ADS>(c1, c2);}

AD_TEMPLATE ADV operator*(const ADV& c1,  const   CS&   s) { return operator* <ADS,ADS,constreal>(c1, s );}
AD_TEMPLATE ADV operator*(const  CV& c1,  const  ADS&   s) { return operator* <ADS,constreal,ADS>(c1, s );}
AD_TEMPLATE ADV operator*(const ADS&  s,  const   CV&  c1) { return operator* <ADS,ADS,constreal>( s, c1);}
AD_TEMPLATE ADV operator*(const  CS&  s,  const  ADV&  c1) { return operator* <ADS,constreal,ADS>( s, c1);}
AD_TEMPLATE ADV operator/(const ADV& c1,  const   CS&   s) { return operator/ <ADS,ADS,constreal>(c1, s );}
AD_TEMPLATE ADV operator/(const  CV& c1,  const  ADS&   s) { return operator/ <ADS,constreal,ADS>(c1, s );}

AD_TEMPLATE ADS       dot(const ADV& c1,  const   CV&  c2) { return       dot <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADS       dot(const  CV& c1,  const  ADV&  c2) { return       dot <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADS      dist(const ADV& c1,  const   CV&  c2) { return      dist <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADS      dist(const  CV& c1,  const  ADV&  c2) { return      dist <ADS,constreal,ADS>(c1, c2);}
AD_TEMPLATE ADS     angle(const ADV& c1,  const   CV&  c2) { return     angle <ADS,ADS,constreal>(c1, c2);}
AD_TEMPLATE ADS     angle(const  CV& c1,  const  ADV&  c2) { return     angle <ADS,constreal,ADS>(c1, c2);}

//trivariate, 6 versions per function, ouch

AD_TEMPLATE ADV      lerp(const  CV& c1,  const  ADV&  c2, const ADS& s ) { return     lerp <ADS,constreal,ADS,ADS>(c1, c2, s );}
AD_TEMPLATE ADV      lerp(const ADV& c1,  const   CV&  c2, const ADS& s ) { return     lerp <ADS,ADS,constreal,ADS>(c1, c2, s );}
AD_TEMPLATE ADV      lerp(const ADV& c1,  const  ADV&  c2, const  CS& s ) { return     lerp <ADS,ADS,ADS,constreal>(c1, c2, s );}
AD_TEMPLATE ADV      lerp(const ADV& c1,  const   CV&  c2, const  CS& s ) { return     lerp <ADS,ADS,constreal,constreal>(c1, c2, s );}
AD_TEMPLATE ADV      lerp(const  CV& c1,  const  ADV&  c2, const  CS& s ) { return     lerp <ADS,constreal,ADS,constreal>(c1, c2, s );}
AD_TEMPLATE ADV      lerp(const  CV& c1,  const   CV&  c2, const ADS& s ) { return     lerp <ADS,constreal,constreal,ADS>(c1, c2, s );}

AD_TEMPLATE ADV    rotate(const  CV& c1,  const  ADV&  c2, const ADS& s ) { return   rotate <ADS,constreal,ADS,ADS>(c1, c2, s );}
AD_TEMPLATE ADV    rotate(const ADV& c1,  const   CV&  c2, const ADS& s ) { return   rotate <ADS,ADS,constreal,ADS>(c1, c2, s );}
AD_TEMPLATE ADV    rotate(const ADV& c1,  const  ADV&  c2, const  CS& s ) { return   rotate <ADS,ADS,ADS,constreal>(c1, c2, s );}
AD_TEMPLATE ADV    rotate(const ADV& c1,  const   CV&  c2, const  CS& s ) { return   rotate <ADS,ADS,constreal,constreal>(c1, c2, s );}
AD_TEMPLATE ADV    rotate(const  CV& c1,  const  ADV&  c2, const  CS& s ) { return   rotate <ADS,constreal,ADS,constreal>(c1, c2, s );}
AD_TEMPLATE ADV    rotate(const  CV& c1,  const   CV&  c2, const ADS& s ) { return   rotate <ADS,constreal,constreal,ADS>(c1, c2, s );}

AD_TEMPLATE ADS tripleProd(const  CV& c1,  const  ADV&  c2, const ADV& c3) { return tripleProd <ADS,constreal,ADS,ADS>(c1, c2, c3);}
AD_TEMPLATE ADS tripleProd(const ADV& c1,  const   CV&  c2, const ADV& c3) { return tripleProd <ADS,ADS,constreal,ADS>(c1, c2, c3);}
AD_TEMPLATE ADS tripleProd(const ADV& c1,  const  ADV&  c2, const  CV& c3) { return tripleProd <ADS,ADS,ADS,constreal>(c1, c2, c3);}
AD_TEMPLATE ADS tripleProd(const ADV& c1,  const   CV&  c2, const  CV& c3) { return tripleProd <ADS,ADS,constreal,constreal>(c1, c2, c3);}
AD_TEMPLATE ADS tripleProd(const  CV& c1,  const  ADV&  c2, const  CV& c3) { return tripleProd <ADS,constreal,ADS,constreal>(c1, c2, c3);}
AD_TEMPLATE ADS tripleProd(const  CV& c1,  const   CV&  c2, const ADV& c3) { return tripleProd <ADS,constreal,constreal,ADS>(c1, c2, c3);}




#endif
