//
//  MeshIO.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MeshIO__
#define __MeshIO__

#include <iostream>
#include <vector>

#include "VS3D.h"
#include "surftrack.h"

class MeshIO {
 public:
  static bool save(VS3D& vs, const std::string& filename, bool binary = true);
  static bool load(VS3D& vs, const std::string& filename, bool binary = true);

  static bool loadIntoRaw(std::vector<LosTopos::Vec3d>& vs,
                          std::vector<LosTopos::Vec3st>& fs,
                          std::vector<LosTopos::Vec2i>& ls,
                          const std::string& filename, bool binary = true);

  static bool saveOBJ(VS3D& vs, const std::string& filename);
};

#endif /* defined(__MeshIO__) */
