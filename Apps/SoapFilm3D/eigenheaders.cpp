//
//  eigenheaders.cpp
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#include "eigenheaders.h"

Vec3d vc(const LosTopos::Vec3d & v)
{
    return Vec3d(v[0], v[1], v[2]);
}

LosTopos::Vec3d vc(const Vec3d & v)
{
    return LosTopos::Vec3d(v[0], v[1], v[2]);
}

