//
//  Sim.cpp
//
//  Fang Da 2014
//
//

#include "Sim.h"

#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "MeshIO.h"
#include "SimOptions.h"
#include "eigenheaders.h"

#ifdef WIN32
#include <direct.h>
#endif

#include "Colormap.h"
#include "PRRenderer.h"
#include "YImage.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#ifdef WIN32
#include <Windows.h>
#include <glew.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

Sim::Sim(bool verbose)
    : m_verbose(verbose),
      m_scene("unspecified"),
      m_output_directory(""),
      m_vs(NULL),
      m_dt(0),
      m_time(0),
      m_frameid(0),
      m_finished(false),
      m_nearest_vertex(-1),
      m_nearest_edge(-1),
      m_nearest_face(-1),
      m_prrenderer(NULL) {}

Sim::~Sim() {}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  General initialization of a simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::init(const std::string &option_file, bool save_outputs,
               bool headless) {
  // declare and load the options
  Options::addStringOption("scene", "T1");
  Options::addStringOption("load-dir", "");
  Options::addDoubleOption("time-step", 0.01);
  Options::addDoubleOption("simulation-time", 1.0);
  Options::addBooleanOption("implicit-integration", false);
  Options::addBooleanOption("pbd-implicit", false);
  Options::addBooleanOption("RK4-velocity-integration", false);
  Options::addDoubleOption("smoothing-coef", 0.0);
  Options::addDoubleOption("damping-coef", 1.0);
  Options::addDoubleOption("sigma", 1.0);
  Options::addDoubleOption("gravity", 0.0);
  Options::addBooleanOption("fmmtl", true);
  Options::addBooleanOption("looped", true);
  Options::addDoubleOption("radius", 5e-5);
  Options::addDoubleOption("density", 1.32e3);
  Options::addDoubleOption("stretching", 5000.0);
  Options::addDoubleOption("bending", 250.0);

  Options::addBooleanOption("output-png", true);
  Options::addIntegerOption(
      "output-png-every-n-frames",
      0);  // 0 means synching with simulation frame rate (equivalent to 1).
  Options::addBooleanOption("output-mesh", false);
  Options::addIntegerOption(
      "output-mesh-every-n-frames",
      0);  // 0 means synching with simulation frame rate (equivalent to 1).
  Options::addBooleanOption("output-obj", false);
  Options::addIntegerOption(
      "output-obj-every-n-frames",
      0);  // 0 means synching with simulation frame rate (equivalent to 1).
  Options::addDoubleOption("remeshing-resolution", 0.1);
  Options::addIntegerOption("remeshing-iterations", 1);

  Options::addDoubleOption(
      "lostopos-collision-epsilon-fraction",
      1e-4);  // lostopos collision epsilon (fraction of mean edge length)
  Options::addDoubleOption(
      "lostopos-merge-proximity-epsilon-fraction",
      0.02);  // lostopos merge proximity epsilon (fraction of mean edge length)
  Options::addBooleanOption("lostopos-perform-smoothing",
                            false);  // whether or not to perform smoothing
  Options::addDoubleOption(
      "lostopos-max-volume-change-fraction",
      1e-4);  // maximum allowed volume change during a remeshing operation
              // (fraction of mean edge length cubed)
  Options::addDoubleOption("lostopos-min-triangle-angle",
                           3.0);  // min triangle angle (in degrees)
  Options::addDoubleOption("lostopos-max-triangle-angle",
                           177.0);  // max triangle angle (in degrees)
  Options::addDoubleOption("lostopos-large-triangle-angle-to-split",
                           160.0);  // threshold for large angles to be split
  Options::addDoubleOption("lostopos-min-triangle-area-fraction",
                           0.02);  // minimum allowed triangle area (fraction of
                                   // mean edge length squared)
  Options::addBooleanOption("lostopos-t1-transition-enabled",
                            true);  // whether t1 is enabled
  Options::addDoubleOption(
      "lostopos-t1-pull-apart-distance-fraction",
      0.1);  // t1 pull apart distance (fraction of mean edge legnth)
  Options::addBooleanOption(
      "lostopos-smooth-subdivision",
      false);  // whether to use smooth subdivision during remeshing
  Options::addBooleanOption(
      "lostopos-allow-non-manifold",
      true);  // whether to allow non-manifold geometry in the mesh
  Options::addBooleanOption("lostopos-allow-topology-changes",
                            true);  // whether to allow topology changes

  Options::addIntegerOption("mesh-size-n", 2);
  Options::addIntegerOption("mesh-size-m", 2);

  Options::addDoubleOption("foam-burst-interval", 20.0);
  Options::addDoubleOption("foam-burst-start", 10.0);

  Options::parseOptionFile(option_file, true);

  // select the scene
  m_scene = Options::strValue("scene");

  if (save_outputs) {
    std::stringstream output_dir_ss;
    output_dir_ss << "output_" << ::time(NULL);
    m_output_directory = output_dir_ss.str();
#ifdef WIN32
    _mkdir(m_output_directory.c_str());
#else
    mkdir(m_output_directory.c_str(), 0755);
#endif
    std::cout << "Outputing to directory: " << m_output_directory << std::endl;
  }

  if (m_scene == "load") {
    m_load_directory = Options::strValue("load-dir");
    assert(m_load_directory != "");

    VS3D *vs = new VS3D(std::vector<LosTopos::Vec3d>(),
                        std::vector<LosTopos::Vec3st>(),
                        std::vector<LosTopos::Vec2i>());
    MeshIO::load(*vs, m_load_directory + "/mesh000000.rec");

    std::vector<LosTopos::Vec3d> vertices = vs->m_st->pm_positions;
    std::vector<LosTopos::Vec3st> faces = vs->m_st->m_mesh.m_tris;
    std::vector<LosTopos::Vec2i> face_labels =
        vs->m_st->m_mesh.m_triangle_labels;

    m_vs = new VS3D(vertices, faces, face_labels);
    MeshIO::load(*m_vs, m_load_directory + "/mesh000000.rec");
  } else {
    std::vector<LosTopos::Vec3d> vertices;
    std::vector<LosTopos::Vec3st> faces;
    std::vector<LosTopos::Vec2i> face_labels;
    std::vector<size_t> constrained_vertices;
    std::vector<Vec3d> constrained_positions;

    if (m_scene == "sphere")
      m_vs = Scenes::sceneSphere(this, vertices, faces, face_labels,
                                 constrained_vertices, constrained_positions);
    else if (m_scene == "tet")
      m_vs = Scenes::sceneTet(this, vertices, faces, face_labels,
                              constrained_vertices, constrained_positions);
    else if (m_scene == "cube")
      m_vs = Scenes::sceneCube(this, vertices, faces, face_labels,
                               constrained_vertices, constrained_positions);
    else if (m_scene == "sheet")
      m_vs = Scenes::sceneSheet(this, vertices, faces, face_labels,
                                constrained_vertices, constrained_positions);
    else if (m_scene == "barrel")
      m_vs = Scenes::sceneBarrel(this, vertices, faces, face_labels,
                                 constrained_vertices, constrained_positions);
    else if (m_scene == "doublebubble")
      m_vs = Scenes::sceneDoubleBubble(this, vertices, faces, face_labels,
                                       constrained_vertices,
                                       constrained_positions);
    else if (m_scene == "twobubbles")
      m_vs =
          Scenes::sceneTwoBubbles(this, vertices, faces, face_labels,
                                  constrained_vertices, constrained_positions);
    else if (m_scene == "triplejunction")
      m_vs = Scenes::sceneTripleJunction(this, vertices, faces, face_labels,
                                         constrained_vertices,
                                         constrained_positions);
    else if (m_scene == "foaminit")
      m_vs = Scenes::sceneFoamInit(this, vertices, faces, face_labels,
                                   constrained_vertices, constrained_positions);
    else if (m_scene == "foam")
      m_vs = Scenes::sceneFoam(this, vertices, faces, face_labels,
                               constrained_vertices, constrained_positions);
    else if (m_scene == "quadjunction")
      m_vs = Scenes::sceneQuadJunction(this, vertices, faces, face_labels,
                                       constrained_vertices,
                                       constrained_positions);
    else if (m_scene == "constrainedsphere")
      m_vs = Scenes::sceneConstrainedSphere(this, vertices, faces, face_labels,
                                            constrained_vertices,
                                            constrained_positions);
    else if (m_scene == "bubblewand")
      m_vs =
          Scenes::sceneBubbleWand(this, vertices, faces, face_labels,
                                  constrained_vertices, constrained_positions);
    else if (m_scene == "tworingspinching")
      m_vs = Scenes::sceneTwoRingsPinching(this, vertices, faces, face_labels,
                                           constrained_vertices,
                                           constrained_positions);
    else if (m_scene == "pullingfoam")
      m_vs =
          Scenes::scenePullingFoam(this, vertices, faces, face_labels,
                                   constrained_vertices, constrained_positions);
    else if (m_scene == "peanutbubble")
      m_vs = Scenes::scenePeanutBubble(this, vertices, faces, face_labels,
                                       constrained_vertices,
                                       constrained_positions);
    else if (m_scene == "straw")
      m_vs = Scenes::sceneStraw(this, vertices, faces, face_labels,
                                constrained_vertices, constrained_positions);
    else if (m_scene == "carousel")
      m_vs = Scenes::sceneCarousel(this, vertices, faces, face_labels,
                                   constrained_vertices, constrained_positions);
    else if (m_scene == "octahedron")
      m_vs =
          Scenes::sceneOctahedron(this, vertices, faces, face_labels,
                                  constrained_vertices, constrained_positions);

    std::cout << "nv = " << vertices.size() << " nf = " << faces.size()
              << std::endl;
  }

  // prepare to start the simulation
  m_time = 0;
  m_dt = Options::doubleValue("time-step");
  m_finished = false;

  // output the initial frame
  if (m_output_directory != "" && Options::boolValue("output-mesh"))
    MeshIO::save(*m_vs, m_output_directory + "/mesh000000.rec");

  // PR rendering
  if (!headless) m_prrenderer = new PRRenderer(m_vs);

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::step() {
  assert(m_scene != "unspecified");
  assert(m_vs);
  assert(m_vs->surfTrack());
  if (m_verbose)
    std::cout << "Time stepping: t = " << m_time << ", dt = " << m_dt
              << std::endl;

  // scene-specific time stepping
  if (m_scene == "sphere")
    Scenes::stepSphere(m_dt, this, m_vs);
  else if (m_scene == "tet")
    Scenes::stepTet(m_dt, this, m_vs);
  else if (m_scene == "cube")
    Scenes::stepCube(m_dt, this, m_vs);
  else if (m_scene == "sheet")
    Scenes::stepSheet(m_dt, this, m_vs);
  else if (m_scene == "barrel")
    Scenes::stepBarrel(m_dt, this, m_vs);
  else if (m_scene == "doublebubble")
    Scenes::stepDoubleBubble(m_dt, this, m_vs);
  else if (m_scene == "twobubbles")
    Scenes::stepTwoBubbles(m_dt, this, m_vs);
  else if (m_scene == "triplejunction")
    Scenes::stepTripleJunction(m_dt, this, m_vs);
  else if (m_scene == "foaminit")
    Scenes::stepFoamInit(m_dt, this, m_vs);
  else if (m_scene == "foam")
    Scenes::stepFoam(m_dt, this, m_vs);
  else if (m_scene == "quadjunction")
    Scenes::stepQuadJunction(m_dt, this, m_vs);
  else if (m_scene == "constrainedsphere")
    Scenes::stepConstrainedSphere(m_dt, this, m_vs);
  else if (m_scene == "bubblewand")
    Scenes::stepBubbleWand(m_dt, this, m_vs);
  else if (m_scene == "tworingspinching")
    Scenes::stepTwoRingsPinching(m_dt, this, m_vs);
  else if (m_scene == "pullingfoam")
    Scenes::stepPullingFoam(m_dt, this, m_vs);
  else if (m_scene == "peanutbubble")
    Scenes::stepPeanutBubble(m_dt, this, m_vs);
  else if (m_scene == "straw")
    Scenes::stepStraw(m_dt, this, m_vs);
  else if (m_scene == "carousel")
    Scenes::stepCarousel(m_dt, this, m_vs);
  else if (m_scene == "octahedron")
    Scenes::stepOctahedron(m_dt, this, m_vs);

  // general time stepping
  double dt = m_vs->step(m_dt);

  // compute the volumes
#ifdef _DEBUG
  std::vector<double> volumes(m_vs->m_nregion, 0);
  Vec3d xref(0, 0, 0);
  for (size_t i = 0; i < m_vs->mesh().nt(); i++) {
    LosTopos::Vec3st t = m_vs->mesh().get_triangle(i);
    Vec3d x0 = m_vs->pos(t[0]);
    Vec3d x1 = m_vs->pos(t[1]);
    Vec3d x2 = m_vs->pos(t[2]);
    double v = (x0 - xref).cross(x1 - xref).dot(x2 - xref);

    LosTopos::Vec2i l = m_vs->mesh().get_triangle_label(i);
    volumes[l[0]] += v;
    volumes[l[1]] -= v;
  }

  std::cerr << m_time;
  for (size_t i = 0; i < volumes.size(); i++) std::cerr << " " << volumes[i];
  std::cerr << std::endl;
#endif
  // advance time
  m_frameid++;
  m_time += m_dt;
  if (m_time >= Options::doubleValue("simulation-time")) m_finished = true;
}

void Sim::stepOutput(bool headless) {
  if (m_output_directory != "") {
    int frameid = (int)(time() / dt() + 0.5);

    int pngfd = Options::intValue("output-png-every-n-frames");
    if ((pngfd == 0 || frameid % pngfd == 0) &&
        Options::boolValue("output-png") && !headless) {
      std::stringstream png_ss;
      png_ss << m_output_directory << "/frame" << std::setfill('0')
             << std::setw(6) << (frameid / std::max(pngfd, 1)) << ".png";

      int w, h;
      w = glutGet(GLUT_WINDOW_WIDTH);
      h = glutGet(GLUT_WINDOW_HEIGHT);

      YImage img;
      img.resize(w, h);
      glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE,
                   (unsigned char *)(img.data()));
      img.flip();
      img.save(png_ss.str().c_str());
    }

    int meshfd = Options::intValue("output-mesh-every-n-frames");
    bool outputmesh = Options::boolValue("output-mesh");
    if ((meshfd == 0 || frameid % meshfd == 0) && outputmesh) {
      std::stringstream mesh_ss;
      mesh_ss << m_output_directory << "/mesh" << std::setfill('0')
              << std::setw(6) << frameid << ".rec";
      MeshIO::save(*m_vs, mesh_ss.str());
    }

    int objfd = Options::intValue("output-obj-every-n-frames");
    if ((objfd == 0 || frameid % objfd == 0) &&
        Options::boolValue("output-obj")) {
      std::stringstream obj_ss;
      obj_ss << m_output_directory << "/mesh" << std::setfill('0')
             << std::setw(6) << frameid << ".obj";
      MeshIO::saveOBJ(*m_vs, obj_ss.str());
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Loading saved simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::load(int inc) {
  int current_frame = (int)((m_time + m_dt * 0.5) / m_dt);
  int next_frame = current_frame + inc;
  if (next_frame < 0) next_frame = 0;

  std::stringstream ss;
  ss << m_load_directory << "/mesh" << std::setfill('0') << std::setw(6)
     << next_frame << ".rec";
  if (!MeshIO::load(*m_vs, ss.str())) {
    std::cout << "Loading frame " << ss.str() << " unsuccessful." << std::endl;
    return false;
  }

  m_time = m_dt * next_frame;

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Rendering
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace {
bool is_edge_nonmanifold(const LosTopos::SurfTrack &st, size_t e) {
  return st.m_mesh.m_edge_to_triangle_map[e].size() != 2;
}

bool is_vertex_nonmanifold(const LosTopos::SurfTrack &st, size_t v) {
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
    if (is_edge_nonmanifold(st, st.m_mesh.m_vertex_to_edge_map[v][i]))
      return true;
  return false;
}

bool is_face_next_to_nonmanifold_vertices(const LosTopos::SurfTrack &st,
                                          size_t f) {
  return is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][0]) ||
         is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][1]) ||
         is_vertex_nonmanifold(st, st.m_mesh.m_tris[f][2]);
}

bool is_edge_next_to_nonmanifold_vertices(const LosTopos::SurfTrack &st,
                                          size_t e) {
  return is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][0]) ||
         is_vertex_nonmanifold(st, st.m_mesh.m_edges[e][1]);
}

bool is_vertex_next_to_nonmanifold_vertices(const LosTopos::SurfTrack &st,
                                            size_t v) {
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++)
    if (is_edge_next_to_nonmanifold_vertices(
            st, st.m_mesh.m_vertex_to_edge_map[v][i]))
      return true;
  return false;
}
}  // namespace

void Sim::render(RenderMode rm, const Vec2d &mousepos, int selection_mask) {
  if (rm == RM_PR) {
    m_prrenderer->render();
    return;
  }

  // find the primitive being picked by cursor
  Mat4d MV;
  Mat4d PJ;
  {
    float mv[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mv);
    float pj[16];
    glGetFloatv(GL_PROJECTION_MATRIX, pj);
    MV << mv[0], mv[4], mv[8], mv[12], mv[1], mv[5], mv[9], mv[13], mv[2],
        mv[6], mv[10], mv[14], mv[3], mv[7], mv[11], mv[15];
    PJ << pj[0], pj[4], pj[8], pj[12], pj[1], pj[5], pj[9], pj[13], pj[2],
        pj[6], pj[10], pj[14], pj[3], pj[7], pj[11], pj[15];
  }
  Mat4d MVP = PJ * MV;

  double mind = -1;
  m_nearest_vertex = -1;
  m_nearest_edge = -1;
  m_nearest_face = -1;
  if (selection_mask & SM_VERTEX) {
    for (size_t i = 0; i < m_vs->mesh().nv(); i++) {
      Vec3d pos = m_vs->pos(i);
      Vec4d scrpos_h = MVP * Vec4d(pos.x(), pos.y(), pos.z(), 1.0);
      Vec2d scrpos = Vec2d(scrpos_h.x(), scrpos_h.y()) / scrpos_h.w();

      double distance = (scrpos - mousepos).norm();
      if (distance < mind || mind < 0) {
        mind = distance;
        m_nearest_vertex = i;
      }
    }
  }

  if (selection_mask & SM_EDGE) {
    for (size_t i = 0; i < m_vs->mesh().ne(); i++) {
      Vec3d v0 = m_vs->pos(m_vs->mesh().m_edges[i][0]);
      Vec3d v1 = m_vs->pos(m_vs->mesh().m_edges[i][1]);

      Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
      Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
      Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
      Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();

      double distance = (mousepos - (scrv0 + scrv1) / 2).norm();
      //            double distance = (mousepos - scrv0 - (mousepos -
      //            scrv0).dot(scrv1 - scrv0) * (scrv1 - scrv0) / (scrv1 -
      //            scrv0).squaredNorm()).norm();
      if (distance < mind || mind < 0) {
        mind = distance;
        m_nearest_vertex = -1;
        m_nearest_edge = i;
      }
    }
  }

  if (selection_mask & SM_FACE) {
    for (size_t i = 0; i < m_vs->mesh().nt(); i++) {
      const LosTopos::Vec3st &t = m_vs->mesh().get_triangle(i);
      Vec3d v0 = m_vs->pos(t[0]);
      Vec3d v1 = m_vs->pos(t[1]);
      Vec3d v2 = m_vs->pos(t[2]);

      Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
      Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
      Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
      Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
      Vec4d scrv2_h = MVP * Vec4d(v2.x(), v2.y(), v2.z(), 1.0);
      Vec2d scrv2 = Vec2d(scrv2_h.x(), scrv2_h.y()) / scrv2_h.w();

      double distance = (mousepos - (scrv0 + scrv1 + scrv2) / 3).norm();
      if (distance < mind || mind < 0) {
        mind = distance;
        m_nearest_vertex = -1;
        m_nearest_edge = -1;
        m_nearest_face = i;
      }
    }
  }

  //    assert(mind >= 0);
  //    assert(m_nearest_vertex >= 0 || m_nearest_edge >= 0 || m_nearest_face >=
  //    0);

  //    Colormap cm(Colormap::MATLAB_JET, 256);
  //    double gamma_max = gamma(0);
  //    double gamma_min = gamma(0);
  //    for (size_t i = 0; i < m_st->m_mesh.nt(); i++)
  //    {
  //        if (gamma(i) < gamma_min)
  //            gamma_min = gamma(i);
  //        if (gamma(i) > gamma_max)
  //            gamma_max = gamma(i);
  //    }
  //    if (gamma_max == gamma_min)
  //    {
  //        gamma_max += 1e-4;
  //        gamma_min -= 1e-4;
  //    }

  //    gamma_max = 0.05;
  //    gamma_min = -0.05;

  bool truncate = false;

  glEnable(GL_DEPTH_TEST);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  if (rm != RM_TRANSPARENT && rm != RM_NONMANIFOLD) {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    GLfloat mat_diffuse[] = {0.3, 0.6, 1.0, 0.95};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    GLfloat mat_specular[] = {1.0, 1.0, 1.0, 0.95};
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    GLfloat mat_shininess[] = {50.0};
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

    GLfloat light_ambient[] = {0.05, 0.05, 0.05, 1.0};
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    GLfloat light_direction[] = {1.0, 1.0, 1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
  }

  // pre-compute vertex normals (area-weighted face normals), for
  // RM_OPAQUE_SMOOTH_SHADED mode
  std::vector<Vec3d> vn(m_vs->mesh().nv(), Vec3d(0, 0, 0));
  for (size_t i = 0; i < m_vs->mesh().nt(); i++) {
    LosTopos::Vec3st t = m_vs->surfTrack()->m_mesh.get_triangle(i);
    Vec3d x0 = m_vs->pos(t[0]);
    Vec3d x1 = m_vs->pos(t[1]);
    Vec3d x2 = m_vs->pos(t[2]);

    Vec3d nt = (x1 - x0).cross(x2 - x0);
    if (m_vs->surfTrack()->m_mesh.get_triangle_label(i)[0] <
        m_vs->surfTrack()->m_mesh.get_triangle_label(i)[1])
      nt = -nt;

    vn[t[0]] += nt;
    vn[t[1]] += nt;
    vn[t[2]] += nt;
  }
  for (size_t i = 0; i < m_vs->mesh().nv(); i++) vn[i].normalize();

  if (false) {
    glLineWidth(5);
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_vs->mesh().nt(); i++) {
      LosTopos::Vec3st t = m_vs->surfTrack()->m_mesh.get_triangle(i);
      Vec3d x0 = m_vs->pos(t[0]);
      Vec3d x1 = m_vs->pos(t[1]);
      Vec3d x2 = m_vs->pos(t[2]);
      Vec3d c = (x0 + x1 + x2) / 3;
      Vec3d n = (x1 - x0).cross(x2 - x0).normalized();
      Vec3d eout = c + n * 0.03;
      Vec3d ein = c - n * 0.03;

      LosTopos::Vec2i l = m_vs->mesh().get_triangle_label(i);

      if (l[1] == 0)
        glColor3d(1, 0, 0);
      else if (l[1] == 1)
        glColor3d(0, 1, 0);
      else if (l[1] == 2)
        glColor3d(0, 0, 1);
      else if (l[1] == 3)
        glColor3d(0.8, 0.8, 0);
      else if (l[1] == 4)
        glColor3d(0.8, 0, 0.8);
      else if (l[1] == 5)
        glColor3d(0, 0.8, 0.8);
      else if (l[1] == 6)
        glColor3d(0.4, 0.4, 1);
      else if (l[1] == 7)
        glColor3d(0.4, 1, 0.4);
      else if (l[1] == 8)
        glColor3d(1, 0.4, 0.4);
      else
        glColor3d(0, 0, 0);
      glVertex3d(c[0], c[1], c[2]);
      glVertex3d(eout[0], eout[1], eout[2]);

      if (l[0] == 0)
        glColor3d(1, 0, 0);
      else if (l[0] == 1)
        glColor3d(0, 1, 0);
      else if (l[0] == 2)
        glColor3d(0, 0, 1);
      else if (l[0] == 3)
        glColor3d(0.8, 0.8, 0);
      else if (l[0] == 4)
        glColor3d(0.8, 0, 0.8);
      else if (l[0] == 5)
        glColor3d(0, 0.8, 0.8);
      else if (l[0] == 6)
        glColor3d(0.4, 0.4, 1);
      else if (l[0] == 7)
        glColor3d(0.4, 1, 0.4);
      else if (l[0] == 8)
        glColor3d(1, 0.4, 0.4);
      else
        glColor3d(0, 0, 0);
      glVertex3d(c[0], c[1], c[2]);
      glVertex3d(ein[0], ein[1], ein[2]);
    }
    glEnd();
    glLineWidth(1);
  }

  glBegin(GL_TRIANGLES);
  for (size_t i = 0; i < m_vs->surfTrack()->m_mesh.nt(); i++) {
    if (rm == RM_NONMANIFOLD) {
      if (!is_face_next_to_nonmanifold_vertices(*m_vs->surfTrack(), i))
        continue;
    }

    if (m_vs->m_st->triangle_is_all_solid(i)) continue;

    //        double r, g, b;
    //        cm.getColorByDensity((gamma(i) - gamma_min) / (gamma_max -
    //        gamma_min), r, g, b); glColor3d(r, g, b);

    LosTopos::Vec3st t = m_vs->surfTrack()->m_mesh.get_triangle(i);
    LosTopos::Vec3d x0 = m_vs->surfTrack()->pm_positions[t[0]];
    LosTopos::Vec3d x1 = m_vs->surfTrack()->pm_positions[t[1]];
    LosTopos::Vec3d x2 = m_vs->surfTrack()->pm_positions[t[2]];
    LosTopos::Vec3d c = (x0 + x1 + x2) / 3;

    if (truncate)
      if (c[0] > 0.5 || c[0] < -0.5) continue;

    double shrink = (rm == RM_TRANSPARENT ? 0.05 : 0);
    x0 += (c - x0) * shrink;
    x1 += (c - x1) * shrink;
    x2 += (c - x2) * shrink;
    Vec3d n0, n1, n2;

    if (rm == RM_OPAQUE_FLAT_SHADED) {
      n0 = vc(cross(x1 - x0, x2 - x0));
      if (m_vs->surfTrack()->m_mesh.get_triangle_label(i)[0] <
          m_vs->surfTrack()->m_mesh.get_triangle_label(i)[1])
        n0 = -n0;
      n0.normalize();
      n1 = n0;
      n2 = n0;
    } else {
      n0 = vn[t[0]];
      n1 = vn[t[1]];
      n2 = vn[t[2]];
    }

    if (m_nearest_face == i)
      glColor4d(0.4, 0.5, 0.6, 0.5);
    else
      glColor4d(0.7, 0.8, 0.9, 0.2);

    glNormal3d(n0[0], n0[1], n0[2]);
    glVertex3d(x0[0], x0[1], x0[2]);
    glNormal3d(n1[0], n1[1], n1[2]);
    glVertex3d(x1[0], x1[1], x1[2]);
    glNormal3d(n2[0], n2[1], n2[2]);
    glVertex3d(x2[0], x2[1], x2[2]);
  }
  glEnd();

  glDisable(GL_BLEND);
  glDepthMask(GL_TRUE);

  if (rm != RM_TRANSPARENT && rm != RM_NONMANIFOLD) {
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
  }

  // render edges
  if (rm != RM_OPAQUE_SMOOTH_SHADED) {
    glLineWidth(1);
    glColor3d(0, 0, 0);
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_vs->mesh().ne(); i++) {
      if (rm == RM_NONMANIFOLD) {
        //                if
        //                (!is_edge_next_to_nonmanifold_vertices(*m_vs->surfTrack(),
        //                i))
        if (!is_edge_nonmanifold(*m_vs->surfTrack(), i)) continue;
      }

      Vec3d x0 = m_vs->pos(m_vs->mesh().m_edges[i][0]);
      Vec3d x1 = m_vs->pos(m_vs->mesh().m_edges[i][1]);

      if (truncate)
        if (x0.x() + x1.x() > 1 || x0.x() + x1.x() < -1) continue;

      glVertex3d(x0[0], x0[1], x0[2]);
      glVertex3d(x1[0], x1[1], x1[2]);
    }
    glEnd();

    glLineWidth(5);
    glColor3d(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_vs->mesh().ne(); i++) {
      if (m_vs->isVertexConstrained(m_vs->mesh().m_edges[i][0]) &&
          m_vs->isVertexConstrained(m_vs->mesh().m_edges[i][1])) {
        if (rm == RM_NONMANIFOLD) {
          //                if
          //                (!is_edge_next_to_nonmanifold_vertices(*m_vs->surfTrack(),
          //                i))
          if (!is_edge_nonmanifold(*m_vs->surfTrack(), i)) continue;
        }

        Vec3d x0 = m_vs->pos(m_vs->mesh().m_edges[i][0]);
        Vec3d x1 = m_vs->pos(m_vs->mesh().m_edges[i][1]);

        if (truncate)
          if (x0.x() + x1.x() > 1 || x0.x() + x1.x() < -1) continue;

        glVertex3d(x0[0], x0[1], x0[2]);
        glVertex3d(x1[0], x1[1], x1[2]);
      }
    }
    glEnd();
    glLineWidth(1);
  }

  // render vertices, with mean curvature coloring
  if (rm != RM_OPAQUE_SMOOTH_SHADED) {
    glPointSize(2);
    glColor3f(0, 0, 0);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < m_vs->surfTrack()->m_mesh.nv(); i++) {
      if (rm == RM_NONMANIFOLD) {
        if (!is_vertex_nonmanifold(*m_vs->surfTrack(), i)) continue;
      }

      const LosTopos::Vec3d &x = m_vs->surfTrack()->pm_positions[i];

      if (truncate)
        if (x[0] > 0.5 || x[0] < -0.5) continue;

      glVertex3d(x[0], x[1], x[2]);
    }
    glEnd();
    glPointSize(1);
  }

  // render constrained vertices
  if (rm != RM_OPAQUE_SMOOTH_SHADED) {
    glPointSize(8);
    glColor3f(0, 0.7, 1);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < m_vs->constrainedPositions().size(); i++) {
      Vec3d x = m_vs->constrainedPositions()[i];
      glVertex3d(x[0], x[1], x[2]);
    }
    glEnd();
    glPointSize(10);
    glColor3f(1, 0.7, 0);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < m_vs->mesh().nv(); i++) {
      Vec3d x = m_vs->pos(i);
      if (m_vs->m_st->vertex_is_any_solid(i)) glVertex3d(x[0], x[1], x[2]);
    }
    glEnd();
    glPointSize(1);
  }

  if (m_nearest_vertex >= 0) {
    glPointSize(12);
    glColor3f(0, 0, 0);
    glBegin(GL_POINTS);
    Vec3d x = m_vs->pos(m_nearest_vertex);
    glVertex3d(x[0], x[1], x[2]);
    glEnd();
  }

  // render gamma
  if (rm != RM_OPAQUE_SMOOTH_SHADED) {
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_vs->mesh().nt(); i++) {
      if (rm == RM_NONMANIFOLD) {
        if (!is_face_next_to_nonmanifold_vertices(*m_vs->surfTrack(), i))
          continue;
      }

      LosTopos::Vec3st t = m_vs->mesh().get_triangle(i);
      if (m_vs->m_st->vertex_is_any_solid(t[0]) &&
          m_vs->m_st->vertex_is_any_solid(t[1]) &&
          m_vs->m_st->vertex_is_any_solid(t[2]))
        continue;  // all-solid faces don't contribute vorticity.

      LosTopos::Vec2i l = m_vs->mesh().get_triangle_label(i);
      Vec3d x0 = m_vs->pos(t[0]);
      Vec3d x1 = m_vs->pos(t[1]);
      Vec3d x2 = m_vs->pos(t[2]);

      if (truncate)
        if (x0.x() + x1.x() + x2.x() > 1.5 || x0.x() + x1.x() + x2.x() < -1.5)
          continue;

      Vec3d e01 = x1 - x0;
      Vec3d e12 = x2 - x1;
      Vec3d e20 = x0 - x2;

      Vec3d gamma = -(e01 * (*m_vs->m_Gamma)[t[2]].get(l) +
                      e12 * (*m_vs->m_Gamma)[t[0]].get(l) +
                      e20 * (*m_vs->m_Gamma)[t[1]].get(l));

      double s = 0.2;
      Vec3d b = (x0 + x1 + x2) / 3;
      Vec3d e = b + s * gamma;

      glColor3d(1, 0, 0);
      glVertex3d(b[0], b[1], b[2]);
      glColor3d(0, 0, 1);
      glVertex3d(e[0], e[1], e[2]);
    }
    glEnd();
  }

  if (m_nearest_edge >= 0) {
    glColor3d(0, 0, 0);
    glLineWidth(3);
    glBegin(GL_LINES);
    Vec3d x0 = m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][0]);
    Vec3d x1 = m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][1]);
    glVertex3d(x0[0], x0[1], x0[2]);
    glVertex3d(x1[0], x1[1], x1[2]);
    glEnd();
    glLineWidth(1);
  }
}

void Sim::showPrimitiveInfo() {
  if (m_nearest_vertex >= 0) {
    std::cout << "Vertex of Interest: " << m_nearest_vertex << " ("
              << m_vs->pos(m_nearest_vertex).transpose() << ")" << std::endl;
    std::cout << "  circulation = " << std::endl
              << m_vs->Gamma(m_nearest_vertex).values << std::endl;
    std::cout << "  incident edges:";
    for (size_t i = 0;
         i < m_vs->mesh().m_vertex_to_edge_map[m_nearest_vertex].size(); i++)
      std::cout << " "
                << m_vs->mesh().m_vertex_to_edge_map[m_nearest_vertex][i];
    std::cout << std::endl;
    std::cout << "  incident faces:";
    for (size_t i = 0;
         i < m_vs->mesh().m_vertex_to_triangle_map[m_nearest_vertex].size();
         i++)
      std::cout << " "
                << m_vs->mesh().m_vertex_to_triangle_map[m_nearest_vertex][i];
    std::cout << std::endl;
    //        std::cout << "  ====== " << (*m_vs->m_Gamma)[m_nearest_vertex][0]
    //        - (*m_vs->m_Gamma)[m_nearest_vertex][1] << std::endl; if
    //        (m_vs->m_dbg_v1.size() == m_vs->mesh().nv())
    //            std::cout << "  ====== " <<
    //            m_vs->m_dbg_v1[m_nearest_vertex][0] << std::endl;
  }

  if (m_nearest_edge >= 0) {
    std::cout << "Edge of Interest: " << m_nearest_edge << ": "
              << m_vs->mesh().m_edges[m_nearest_edge][0] << " ("
              << m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][0]).transpose()
              << ") - " << m_vs->mesh().m_edges[m_nearest_edge][1] << " ("
              << m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][1]).transpose()
              << ") length = "
              << (m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][1]) -
                  m_vs->pos(m_vs->mesh().m_edges[m_nearest_edge][0]))
                     .norm()
              << std::endl;
    // std::cout << "  dbg e1 = " << m_vs->m_dbg_e1[m_nearest_edge][0] << " " <<
    // m_vs->m_dbg_e1[m_nearest_edge][1] << " " <<
    // m_vs->m_dbg_e1[m_nearest_edge][2] << std::endl;
    std::cout << "  incident faces:";
    for (size_t i = 0;
         i < m_vs->mesh().m_edge_to_triangle_map[m_nearest_edge].size(); i++)
      std::cout << " "
                << m_vs->mesh().m_edge_to_triangle_map[m_nearest_edge][i];
    std::cout << std::endl;
  }

  if (m_nearest_face >= 0) {
    std::cout << "Face of Interest: " << m_nearest_face << ": "
              << m_vs->mesh().m_tris[m_nearest_face][0] << " ("
              << m_vs->pos(m_vs->mesh().m_tris[m_nearest_face][0]).transpose()
              << "), " << m_vs->mesh().m_tris[m_nearest_face][1] << " ("
              << m_vs->pos(m_vs->mesh().m_tris[m_nearest_face][1]).transpose()
              << "), " << m_vs->mesh().m_tris[m_nearest_face][2] << " ("
              << m_vs->pos(m_vs->mesh().m_tris[m_nearest_face][2]).transpose()
              << ")" << std::endl;
    std::cout << "  labels: " << m_vs->mesh().get_triangle_label(m_nearest_face)
              << std::endl;
    std::cout << "  incident edges:";
    for (size_t i = 0; i < 3; i++)
      std::cout << " "
                << m_vs->mesh().m_triangle_to_edge_map[m_nearest_face][i];
    std::cout << std::endl;
  }
}
