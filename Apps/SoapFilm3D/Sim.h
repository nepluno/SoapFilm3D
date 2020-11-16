//
//  Sim.h
//
//  Fang Da 2014
//
//

#ifndef __Sim__
#define __Sim__

#include <iostream>
#include <string>

#include "PRRenderer.h"
#include "Scenes.h"
#include "VS3D.h"

class Sim {
  friend class Scenes;

 public:
  Sim(bool verbose);
  ~Sim();

  VS3D* vs() { return m_vs; }

 public:
  bool init(const std::string& option_file, bool save_outputs, bool headless);

 public:
  void step();
  void stepOutput(bool headless);
  bool isFinished() const { return m_finished; }

  double dt() const { return m_dt; }
  double time() const { return m_time; }

  bool load(int inc = 1);

 public:
  enum RenderMode {
    RM_TRANSPARENT,
    RM_NONMANIFOLD,
    RM_OPAQUE_FLAT_SHADED,
    RM_OPAQUE_SMOOTH_SHADED,
    RM_PR,

    RM_COUNT
  };

  static const int SM_VERTEX = 0x01;
  static const int SM_EDGE = 0x02;
  static const int SM_FACE = 0x04;
  void render(RenderMode rm, const Vec2d& mousepos,
              int selection_mask = SM_VERTEX | SM_EDGE | SM_FACE);

  void showPrimitiveInfo();

 protected:
  bool m_verbose;

  std::string m_scene;
  std::string m_output_directory;
  std::string m_load_directory;

  VS3D* m_vs;

  double m_dt;
  double m_time;
  int m_frameid;
  bool m_finished;

  int m_nearest_vertex;
  int m_nearest_edge;
  int m_nearest_face;

  PRRenderer* m_prrenderer;
};

#endif
