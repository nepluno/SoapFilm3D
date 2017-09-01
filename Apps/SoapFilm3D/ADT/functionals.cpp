#include "adreal.h"
#include "cvec3t.h"
#include "advec.h"


#define NDOF_DS 12
#define N_ORDER 1


// functions here are not tested and were intended primarily to 
// verify that the adreal.h and advec.h provide all 
// necessary functionality to implement them


typedef CVec3T<double> vec3;
typedef adreal<NDOF_DS,N_ORDER,double>  adrealDS;
typedef CVec3T<adreal<NDOF_DS,N_ORDER,double> > advecDS;

//    p2
//q1     q2
//    p1

adrealDS EDS(const vec3& p01, const vec3& p02, const vec3& q01, const vec3& q02) { 
  advecDS p1,p2;
  advecDS q1,q2;
  set_indepedent(p1,p01,0);   set_indepedent(p2,p02,3);
  set_indepedent(q1,q01,6);   set_indepedent(q2,q02,9); 

  advecDS n1 = cross(p2-p1,q1-p1), n2 = cross(q2-p1,p2-p1); 
  adrealDS n1len = len(n1), n2len = len(n2);
  n1 = n1/n1len; n2 = n2/n2len; 
  
  // correct for all cases
  //  adreal theta = atan2(dot( cross(n1,n2),dir(p2-p1)), dot(n1,n2));
  // ok for flat undeformed, as theta sign does not matter
  adrealDS theta = acos(dot(n1,n2));
  adrealDS l2 = lenSq(p2-p1);
  return 6.0*l2*sqr(theta)/(n1len + n2len);
}


#define NDOF_MEMBRANE 9

typedef adreal<NDOF_MEMBRANE,N_ORDER,double> adrealS;
typedef CVec3T<adreal<NDOF_MEMBRANE,N_ORDER,double> > advecS;

adrealS Emembrane(const vec3 p0[3], const vec3 pundef[3]) { 
  advecS p[3];
  set_indepedent(p[0],p0[0],0);   set_indepedent(p[1],p0[1],3);   set_indepedent(p[2],p0[2],6);

  
}


#define NDOF_MIDEDGE 12
typedef adreal<NDOF_MIDEDGE,N_ORDER,double> adrealME;
typedef CVec3T<adreal<NDOF_MIDEDGE,N_ORDER,double> > advecME;

adrealME Emidedge(const vec3 p0[3], const vec3 q0[3], const vec3& xi0, double nu) { 
  advecME p[3];
  advecME xi;
  set_indepedent(p[0],p0[0],0);   set_indepedent(p[1],p0[1],3);   set_indepedent(p[2],p0[2],6);
  set_indepedent(xi, xi0,9);

  // some work can be saved if these are precomputed, but probably not a
  // significant percentage
  vec3 na0[3] = { dir(cross (p0[2]-p0[0],q0[1]-p0[0])), dir(cross ( p0[0]-p0[1], q0[2]-p0[1])), dir(cross( p0[1]-p0[2], q0[0]-p0[2])) };
  vec3 n0 = dir(cross( p0[1]- p0[0], p0[2]-p0[0]));
  vec3 tau0[3] = { dir(na0[0]+n0), dir(na0[1]+n0), dir(na0[2]+n0)}; 

  advecME v[3] = { p[2]-p[1],p[0]-p[2],p[1]-p[0]};
  advecME n = cross(v[2],v[1]);
  adrealME nleninv = 1.0/len(n);
  n = n*nleninv;
  advecME t[3] = { cross(v[0],n), cross(v[1],n), cross(v[2],n) };
  
  adrealME T[6] = // 00, 11, 22, 01, 02, 12
    { sqr(lenSq(v[0])), sqr(lenSq(v[1])), sqr(len(v[2])), 
      nu*lenSq(v[0])*lenSq(v[1]) + (1-nu)*sqr(dot(v[0],v[1])), 
      nu*lenSq(v[0])*lenSq(v[2]) + (1-nu)*sqr(dot(v[0],v[2])), 
      nu*lenSq(v[1])*lenSq(v[2]) + (1-nu)*sqr(dot(v[1],v[2])) 
    };

  adrealME f[3] = { ( xi[0] - dot(n,tau0[0]))/dot(t[0],tau0[0]), ( xi[1] - dot(n,tau0[1]))/dot(t[1],tau0[1]), ( xi[2] - dot(n,tau0[2]))/dot(t[2],tau0[2])};

  return sqr(nleninv)*(sqr(f[0])*T[0] + sqr(f[1])*T[1] + sqr(f[2])*T[2] + 2.0*f[0]*f[1]*T[3] + 2.0*f[0]*f[2]*T[4] + 2.0*f[1]*f[2]*T[5]);
}

#define NDOF_AVG 18
typedef adreal<NDOF_MIDEDGE,N_ORDER,double> adrealAVG;
typedef CVec3T<adreal<NDOF_MIDEDGE,N_ORDER,double> > advecAVG;

adrealAVG Eaveraged(const vec3 p0[3], const vec3 q0[3], double nu) { 
  advecAVG p[3];
  advecAVG q[3];
  set_indepedent(p[0],p0[0],0);   set_indepedent(p[1],p0[1],3);   set_indepedent(p[2],p0[2],6);
  set_indepedent(q[0],q0[0],9);   set_indepedent(q[1],q0[1],12);   set_indepedent(q[2],q0[2],15);

  advecAVG na[3] = { dir(cross (p[2]-p[0],q[1]-p[0])), dir(cross ( p[0]-p[1], q[2]-p[1])), dir(cross( p[1]-p[2], q[0]-p[2])) };
  advecAVG v[3] = { p[2]-p[1],p[0]-p[2],p[1]-p[0]};
  adrealAVG linv[3] = {1.0/len(v[0]), 1.0/len(v[1]),1.0/len(v[2])};
  advecAVG n = cross(v[2],v[1]);
  adrealAVG nleninv = 1.0/len(n);
  n = n*nleninv;

  adrealAVG f[3] = {acos(dot(na[0],n))*linv[0], acos(dot(na[1],n))*linv[1], acos(dot(na[2],n))*linv[2]};
  
  adrealAVG T[6] = // 00, 11, 22, 01, 02, 12
    { sqr(lenSq(v[0])), sqr(lenSq(v[1])), sqr(lenSq(v[2])), 
      nu*lenSq(v[0])*lenSq(v[1]) + (1-nu)*sqr(dot(v[0],v[1])), 
      nu*lenSq(v[0])*lenSq(v[2]) + (1-nu)*sqr(dot(v[0],v[2])), 
      nu*lenSq(v[1])*lenSq(v[2]) + (1-nu)*sqr(dot(v[1],v[2])) 
    };

  return sqr(nleninv)*(sqr(f[0])*T[0] + sqr(f[1])*T[1] + sqr(f[2])*T[2] + 2.0*f[0]*f[1]*T[3] + 2.0*f[0]*f[2]*T[4] + 2.0*f[1]*f[2]*T[5]);
}


