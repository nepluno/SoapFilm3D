//
//  eigenheaders.h
//  MultiTracker
//
//  Created by Fang Da on 10/24/14.
//
//

#ifndef MultiTracker_eigenheaders_h
#define MultiTracker_eigenheaders_h

#include <Eigen/Core>
#include <Eigen/Dense>
#include "surftrack.h"

typedef Eigen::Matrix<double, 4, 4> Mat4d;
typedef Eigen::Matrix<double, 3, 3> Mat3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatXd;
typedef Eigen::Matrix<double, 4, 1> Vec4d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 2, 1> Vec2d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecXd;

typedef Eigen::Matrix<int, 3, 1> Vec3i;
typedef Eigen::Matrix<int, 2, 1> Vec2i;

Vec3d vc(const LosTopos::Vec3d & v);
LosTopos::Vec3d vc(const Vec3d & v);

class Vec2iComp
{
public:
    bool operator () (const Vec2i & v1, const Vec2i & v2) const { return v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]); }
};

#endif
