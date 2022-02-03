//
//  VS3D.cpp
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#include "VS3D.h"

#ifdef PRINT_TIMING
#include <chrono>
#endif

#include "LinearBendingForce.h"
#include "LinearizedImplicitEuler.h"
#include "LosTopos/LosTopos3D/subdivisionscheme.h"
#include "SimOptions.h"
#include "SimpleGravityForce.h"
#include "SpringForce.h"
#include "VertexAreaForce.h"

VecXd BiotSavartFunc(VS3D& vs, const VecXd& dx);

bool VS3D::isVertexConstrained(size_t vert) {
  return m_constrained_mapping.find(vert) != m_constrained_mapping.end();
}

VS3D::VS3D(const std::vector<LosTopos::Vec3d>& vs,
           const std::vector<LosTopos::Vec3st>& fs,
           const std::vector<LosTopos::Vec2i>& ls,
           const std::vector<size_t>& constrained_vertices,
           const std::vector<Vec3d>& constrained_positions,
           const std::vector<Vec3d>& constrained_velocities,
           const std::vector<unsigned char>& constrained_fixed) {
  // load sim options
  m_sim_options.implicit = Options::boolValue("implicit-integration");
  m_sim_options.pbd = Options::boolValue("pbd-implicit");
  m_sim_options.smoothing_coef = Options::doubleValue("smoothing-coef");
  m_sim_options.damping_coef = Options::doubleValue("damping-coef");
  m_sim_options.sigma = Options::doubleValue("sigma");
  m_sim_options.gravity = Options::doubleValue("gravity");
  m_sim_options.looped = Options::boolValue("looped");
  m_sim_options.radius = Options::doubleValue("radius");
  m_sim_options.density = Options::doubleValue("density");
  m_sim_options.stretching = Options::doubleValue("stretching");
  m_sim_options.bending = Options::doubleValue("bending");

  // construct the surface tracker
  double mean_edge_len = Options::doubleValue("remeshing-resolution");
  if (mean_edge_len == 0) {
    for (size_t i = 0; i < fs.size(); i++) {
      mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
      mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
      mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
    }
    mean_edge_len /= (fs.size() * 3);
  }
  double min_edge_len = mean_edge_len * 0.5;
  double max_edge_len = mean_edge_len * 1.5;
#ifdef _DEBUG
  std::cout << "mean edge length = " << mean_edge_len
            << " min edge length = " << min_edge_len
            << " max edge length = " << max_edge_len << std::endl;
#endif
  LosTopos::SurfTrackInitializationParameters params;
  params.m_proximity_epsilon =
      Options::doubleValue("lostopos-collision-epsilon-fraction") *
      mean_edge_len;
  params.m_merge_proximity_epsilon =
      Options::doubleValue("lostopos-merge-proximity-epsilon-fraction") *
      mean_edge_len;
  params.m_allow_vertex_movement_during_collapse = true;
  params.m_perform_smoothing = Options::boolValue("lostopos-perform-smoothing");
  params.m_min_edge_length = min_edge_len;
  params.m_max_edge_length = max_edge_len;
  params.m_max_volume_change =
      Options::doubleValue("lostopos-max-volume-change-fraction") *
      pow(mean_edge_len, 3);
  params.m_min_triangle_angle =
      Options::doubleValue("lostopos-min-triangle-angle");
  params.m_max_triangle_angle =
      Options::doubleValue("lostopos-max-triangle-angle");
  params.m_large_triangle_angle_to_split =
      Options::doubleValue("lostopos-large-triangle-angle-to-split");
  params.m_min_triangle_area =
      Options::doubleValue("lostopos-min-triangle-area-fraction") *
      pow(mean_edge_len, 2);
  params.m_verbose = false;
  params.m_allow_non_manifold =
      Options::boolValue("lostopos-allow-non-manifold");
  params.m_allow_topology_changes =
      Options::boolValue("lostopos-allow-topology-changes");
  params.m_collision_safety = true;
  params.m_remesh_boundaries = true;
  params.m_t1_transition_enabled =
      Options::boolValue("lostopos-t1-transition-enabled");
  params.m_pull_apart_distance =
      Options::doubleValue("lostopos-t1-pull-apart-distance-fraction") *
      mean_edge_len;

  params.m_velocity_field_callback = NULL;

  if (Options::boolValue("lostopos-smooth-subdivision"))
    params.m_subdivision_scheme = new LosTopos::ModifiedButterflyScheme();
  else
    params.m_subdivision_scheme = new LosTopos::MidpointScheme();

  params.m_use_curvature_when_collapsing = false;
  params.m_use_curvature_when_splitting = false;

  m_constrained_vertices = constrained_vertices;
  m_constrained_positions = constrained_positions;
  m_constrained_mapping.insert(constrained_vertices.begin(),
                               constrained_vertices.end());

  if (m_constrained_positions.size() > 0) {
    m_constrained_velocities = constrained_velocities;
    if (m_constrained_velocities.size() != m_constrained_positions.size()) {
      m_constrained_velocities.resize(m_constrained_positions.size(),
                                      Vec3d(0, 0, 0));
    }

    m_constrained_fixed = constrained_fixed;
    if (m_constrained_fixed.size() != m_constrained_positions.size()) {
      m_constrained_fixed.resize(m_constrained_positions.size(), false);
    }

    m_constrained_mass.resize(m_constrained_positions.size() * 3);
    const size_t np = m_constrained_positions.size();
    if (m_sim_options.looped) {
      const size_t ne = np;
      VecXd rest_length(ne);
      for (size_t i = 0; i < ne; ++i) {
        int iforw = (i == np - 1) ? 0 : (i + 1);
        rest_length(i) =
            (m_constrained_positions[iforw] - m_constrained_positions[i])
                .norm();
      }

      for (size_t i = 0; i < np; ++i) {
        int ieforw = i;
        int ieback = (i == 0) ? (ne - 1) : (i - 1);
        double len = (rest_length(ieforw) + rest_length(ieback)) * 0.5;
        double mass = M_PI * m_sim_options.radius * m_sim_options.radius *
                      m_sim_options.density * len;
        for (size_t r = 0; r < 3; ++r) m_constrained_mass[i * 3 + r] = mass;
      }

      // add spring forces
      for (size_t i = 0; i < ne; ++i) {
        int iforw = (i == np - 1) ? 0 : (i + 1);
        m_forces.push_back(new SpringForce(std::pair<int, int>(i, iforw),
                                           m_sim_options.stretching,
                                           rest_length(i)));
      }

      // add bending forces
      for (size_t i = 0; i < np; ++i) {
        int iforw = (i == np - 1) ? 0 : (i + 1);
        int iback = (i == 0) ? (np - 1) : (i - 1);

        int ieforw = i;
        int ieback = (i == 0) ? (ne - 1) : (i - 1);

        m_forces.push_back(new LinearBendingForce(
            iback, i, iforw, m_sim_options.bending,
            m_sim_options.bending * m_sim_options.damping_coef * 0.0001,
            Vector2s::Zero(), rest_length(ieback), rest_length(ieforw)));
      }
    } else {
      const size_t ne = np - 1;
      VecXd rest_length(ne);
      for (size_t i = 0; i < ne; ++i) {
        rest_length(i) =
            (m_constrained_positions[i + 1] - m_constrained_positions[i])
                .norm();
      }
      for (size_t i = 0; i < np; ++i) {
        double len;
        if (i == 0) {
          len = rest_length(0) * 0.5;
        } else if (i == np - 1) {
          len = rest_length(ne - 1) * 0.5;
        } else {
          int ieforw = i;
          int ieback = i - 1;
          len = (rest_length(ieforw) + rest_length(ieback)) * 0.5;
        }

        double mass = M_PI * m_sim_options.radius * m_sim_options.radius *
                      m_sim_options.density * len;
        for (size_t r = 0; r < 3; ++r) m_constrained_mass[i * 3 + r] = mass;
      }

      for (size_t i = 0; i < ne; ++i) {
        m_forces.push_back(new SpringForce(std::pair<int, int>(i, i + 1),
                                           m_sim_options.stretching,
                                           rest_length(i)));
      }

      for (size_t i = 1; i < np - 1; ++i) {
        int iforw = i + 1;
        int iback = i - 1;

        int ieback = i - 1;
        int ieforw = i;

        m_forces.push_back(new LinearBendingForce(
            iback, i, iforw, m_sim_options.bending,
            m_sim_options.bending * m_sim_options.damping_coef,
            Vector2s::Zero(), rest_length(ieback), rest_length(ieforw)));
      }
    }
  }

  std::vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
  for (size_t i = 0; i < m_constrained_vertices.size(); i++)
    masses[m_constrained_vertices[i]] *=
        std::numeric_limits<double>::infinity();
  m_st = new LosTopos::SurfTrack(vs, fs, ls, masses, params);
  m_st->m_solid_vertices_callback = this;
  m_st->m_mesheventcallback = this;

  // find out the number of regions
  m_nregion = 0;
  for (size_t i = 0; i < mesh().nt(); i++) {
    LosTopos::Vec2i l = mesh().get_triangle_label(i);
    m_nregion = std::max(m_nregion, l[0] + 1);
    m_nregion = std::max(m_nregion, l[1] + 1);
  }

  // initialize the sheet quantities
  m_Gamma = new LosTopos::NonDestructiveTriMesh::VertexData<GammaType>(
      &(m_st->m_mesh));
  for (size_t i = 0; i < mesh().nv(); i++) (*m_Gamma)[i] = GammaType(m_nregion);

  // Biot-Savart kernel regularization parameter
  m_delta = max_edge_len * 0.5;

  m_forces.push_back(
      new SimpleGravityForce(Vector3s(0, 0, -m_sim_options.gravity)));
  m_forces.push_back(new VertexAreaForce(this, m_sim_options.sigma));
  // initialize constraint stepper
  m_constraint_stepper = new LinearizedImplicitEuler();
}

VS3D::~VS3D() {
  if (m_st) delete m_st;
  for (Force* f : m_forces) delete f;

  if (m_constraint_stepper) delete m_constraint_stepper;
}

namespace {
Mat3d skewSymmetric(const Vec3d& v) {
  Mat3d ss = Mat3d::Zero();
  ss(0, 1) = -v(2);
  ss(1, 0) = v(2);
  ss(0, 2) = v(1);
  ss(2, 0) = -v(1);
  ss(1, 2) = -v(0);
  ss(2, 1) = v(0);
  return ss;
}
}  // namespace

void VS3D::stepConstrainted(const scalar& dt) {
  m_constraint_stepper->stepScene(*this, dt);
}

double VS3D::step(double dt) {
  static int counter = 0;
  counter++;

#ifdef PRINT_TIMING
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;

  TimePoint last_time_point;
#endif
  if (counter % 2 == 0) {
    // mesh improvement
#ifdef PRINT_TIMING
    last_time_point = Clock::now();
#endif
    for (int i = 0; i < Options::intValue("remeshing-iterations"); i++) {
      m_st->topology_changes();
      m_st->improve_mesh();
    }
#ifdef PRINT_TIMING
    std::cout << "[mesh improvement] "
              << static_cast<scalar>(
                     std::chrono::duration_cast<std::chrono::nanoseconds>(
                         Clock::now() - last_time_point)
                         .count()) *
                     1e-6
              << " ms" << std::endl;
#endif
    // defrag the mesh in the end, to ensure the next step starts with a clean
    // mesh
#ifdef PRINT_TIMING
    last_time_point = Clock::now();
#endif
    m_st->defrag_mesh_from_scratch(m_constrained_vertices);
#ifdef PRINT_TIMING
    std::cout << "[defrag_mesh_from_scratch] "
              << static_cast<scalar>(
                     std::chrono::duration_cast<std::chrono::nanoseconds>(
                         Clock::now() - last_time_point)
                         .count()) *
                     1e-6
              << " ms" << std::endl;
#endif
#ifdef _DEBUG
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
      assert(m_constrained_vertices[i] < mesh().nv());
#endif

  } else {
    // update gamma due to external forces (surface tension, gravity, etc), and
    // update velocity from gamma
    if (simOptions().implicit) {
      if (simOptions().pbd)
        step_PBD_implicit(dt);
      else
        step_implicit(dt);
    } else {
#ifdef PRINT_TIMING
      last_time_point = Clock::now();
#endif
      step_explicit(dt);
#ifdef PRINT_TIMING
      std::cout << "[step_explicit] "
                << static_cast<scalar>(
                       std::chrono::duration_cast<std::chrono::nanoseconds>(
                           Clock::now() - last_time_point)
                           .count()) *
                       1e-6
                << " ms" << std::endl;
#endif
    }
  }
#ifdef PRINT_TIMING
  last_time_point = Clock::now();
#endif
  // shift all Gammas for each region pair to their global mean
  MatXd means(m_nregion, m_nregion);
  means.setZero();
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> counters(m_nregion,
                                                              m_nregion);
  counters.setZero();

  for (size_t i = 0; i < mesh().nv(); i++) {
    std::set<Vec2i, Vec2iComp> region_pairs;
    for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++) {
      LosTopos::Vec2i l =
          mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
      region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
    }

    for (std::set<Vec2i, Vec2iComp>::iterator j = region_pairs.begin();
         j != region_pairs.end(); j++) {
      Vec2i rp = *j;
      means(rp[0], rp[1]) += (*m_Gamma)[i].get(rp);
      counters(rp[0], rp[1])++;
    }
  }

  for (int i = 0; i < m_nregion; i++)
    for (int j = i + 1; j < m_nregion; j++)
      if (counters(i, j) != 0) means(i, j) /= counters(i, j);

  for (size_t i = 0; i < mesh().nv(); i++) {
    std::set<Vec2i, Vec2iComp> region_pairs;
    for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++) {
      LosTopos::Vec2i l =
          mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
      region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
    }

    for (std::set<Vec2i, Vec2iComp>::iterator j = region_pairs.begin();
         j != region_pairs.end(); j++) {
      Vec2i rp = *j;
      (*m_Gamma)[i].set(rp, (*m_Gamma)[i].get(rp) - means(rp[0], rp[1]));
    }
  }

  std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> >
      incident_region_pairs(mesh().nv());
  for (size_t i = 0; i < mesh().nv(); i++)
    incident_region_pairs[i].setZero(m_nregion, m_nregion);

  for (size_t i = 0; i < mesh().nt(); i++) {
    LosTopos::Vec3st t = mesh().get_triangle(i);
    LosTopos::Vec2i l = mesh().get_triangle_label(i);
    incident_region_pairs[t[0]](l[0], l[1]) =
        incident_region_pairs[t[0]](l[1], l[0]) = true;
    incident_region_pairs[t[1]](l[0], l[1]) =
        incident_region_pairs[t[1]](l[1], l[0]) = true;
    incident_region_pairs[t[2]](l[0], l[1]) =
        incident_region_pairs[t[2]](l[1], l[0]) = true;
  }

  std::vector<GammaType> newGamma = m_Gamma->m_data;
  for (size_t i = 0; i < mesh().nv(); i++) {
    for (int j = 0; j < m_nregion; j++) {
      for (int k = j + 1; k < m_nregion; k++) {
        if (incident_region_pairs[i](j, k)) {
          double neighborhood_mean = 0;
          int neighborhood_counter = 0;
          for (size_t l = 0; l < mesh().m_vertex_to_edge_map[i].size(); l++) {
            LosTopos::Vec2st e =
                mesh().m_edges[mesh().m_vertex_to_edge_map[i][l]];
            size_t vother = (e[0] == i ? e[1] : e[0]);
            if (incident_region_pairs[vother](j, k)) {
              neighborhood_mean += (*m_Gamma)[vother].get(j, k);
              neighborhood_counter++;
            }
          }
          if (neighborhood_counter != 0)
            neighborhood_mean /= neighborhood_counter;
          newGamma[i].set(j, k,
                          (*m_Gamma)[i].get(j, k) +
                              (neighborhood_mean - (*m_Gamma)[i].get(j, k)) *
                                  simOptions().smoothing_coef * dt);
        } else {
          newGamma[i].set(j, k, 0);
        }
      }
    }
  }

  for (size_t i = 0; i < mesh().nv(); i++) (*m_Gamma)[i] = newGamma[i];
#ifdef PRINT_TIMING
  std::cout << "[update Gamma] "
            << static_cast<scalar>(
                   std::chrono::duration_cast<std::chrono::nanoseconds>(
                       Clock::now() - last_time_point)
                       .count()) *
                   1e-6
            << " ms" << std::endl;

  last_time_point = Clock::now();
#endif
  // before enforcing constraints, first scan through the mesh to find any solid
  // vertices not registered as constraints. they can appear due to remeshing
  // (splitting an all-solid edge)
  std::vector<int> constrained_vertices_map(mesh().nv(), -1);
  for (size_t i = 0; i < m_constrained_vertices.size(); i++)
    constrained_vertices_map[m_constrained_vertices[i]] = i;

  // contruct the open boudnary extra faces
  m_obefc.clear();
  m_obefe.clear();
  m_obefn.clear();
  std::vector<size_t> ob;    // open boundary edges
  std::set<size_t> obv_set;  // open boundary vertices
  for (size_t i = 0; i < mesh().ne(); i++) {
    if (surfTrack()->edge_is_all_solid(i))  // this is a constrained edge
    {
      // assert(mesh().m_edge_to_triangle_map[i].size() == 2);
      if (mesh().m_edge_to_triangle_map[i].size() == 2) {
        size_t f0 = mesh().m_edge_to_triangle_map[i][0];
        size_t f1 = mesh().m_edge_to_triangle_map[i][1];
        bool s0 = surfTrack()->triangle_is_all_solid(f0);
        bool s1 = surfTrack()->triangle_is_all_solid(f1);

        if ((s0 && !s1) || (s1 && !s0))  // this is an open boundary edge
        {
          size_t f = (s0 ? f1 : f0);

          size_t v0 = mesh().m_edges[i][0];
          size_t v1 = mesh().m_edges[i][1];
          size_t v2 = mesh().get_third_vertex(i, f);

          Vec3d x0 = pos(v0);
          Vec3d x1 = pos(v1);
          Vec3d x2 = pos(v2);

          m_obefc.push_back((x0 + x1) / 2);
          m_obefe.push_back((x1 - x0).normalized());
          m_obefn.push_back((x1 - x0).cross(x2 - x0).normalized());
          ob.push_back(i);
          obv_set.insert(v0);
          obv_set.insert(v1);
        }
      } else if (mesh().m_edge_to_triangle_map[i].size() == 1) {
        size_t f = mesh().m_edge_to_triangle_map[i][0];

        size_t v0 = mesh().m_edges[i][0];
        size_t v1 = mesh().m_edges[i][1];
        size_t v2 = mesh().get_third_vertex(i, f);

        Vec3d x0 = pos(v0);
        Vec3d x1 = pos(v1);
        Vec3d x2 = pos(v2);

        m_obefc.push_back((x0 + x1) / 2);
        m_obefe.push_back((x1 - x0).normalized());
        m_obefn.push_back((x1 - x0).cross(x2 - x0).normalized());
        ob.push_back(i);
        obv_set.insert(v0);
        obv_set.insert(v1);
      }
    }
  }

  std::vector<size_t> obv;
  obv.assign(obv_set.begin(), obv_set.end());
  assert(obv.size() ==
         ob.size());  // assume the open boundary has simple topology
  size_t nob = ob.size();

  if (nob > 0) {
    std::vector<Vec3d> obvn(nob);  // open boundary vertex normals
    for (size_t i = 0; i < nob; i++) {
      Vec3d n(0, 0, 0);
      for (size_t j = 0; j < nob; j++)
        if (mesh().m_edges[ob[j]][0] == obv[i] ||
            mesh().m_edges[ob[j]][1] == obv[i])
          n += m_obefn[j];
      obvn[i] = n.normalized();
    }

    // solve for the obef vorticities such that the open boundary does not move
    MatXd obA = MatXd::Zero(nob, nob);
    for (size_t i = 0; i < nob; i++) {
      Vec3d x = pos(obv[i]);

      for (size_t j = 0; j < nob; j++) {
        Vec3d xp = m_obefc[j];

        Vec3d dx = x - xp;
        double dxn = sqrt(dx.squaredNorm() + m_delta * m_delta);

        obA(i, j) =
            -(skewSymmetric(dx) * m_obefe[j]).dot(obvn[i]) / (dxn * dxn * dxn);
      }
    }

    obA *= dt / (4 * M_PI);

    VecXd obrhs = VecXd::Zero(nob);
    for (size_t i = 0; i < nob; i++)
      obrhs[i] = (m_constrained_positions[constrained_vertices_map[obv[i]]] -
                  vc(surfTrack()->pm_newpositions[obv[i]]))
                     .dot(obvn[i]);

    //    VecXd obefv = obA.partialPivLu().solve(obrhs);
    double oblambda = 0.1;
    VecXd obefv =
        (obA.transpose() * obA +
         oblambda * oblambda * MatXd::Identity(nob, nob))
            .partialPivLu()
            .solve(obA.transpose() *
                   obrhs);  // regularized solve to avoid blowing up in presence
                            // of near-dependent constraints

    m_obefv.resize(nob, 0);
    for (size_t i = 0; i < nob; i++) m_obefv[i] = obefv[i];
  }
  // following the open boundary solve, recompute the velocities
  VecXd newv = BiotSavartFunc(*this, VecXd::Zero(mesh().nv() * 3));
  for (size_t i = 0; i < mesh().nv(); i++)
    m_st->pm_newpositions[i] =
        m_st->pm_positions[i] + vc(newv.segment<3>(i * 3)) * dt;
#ifdef PRINT_TIMING
  std::cout << "[solve open boundary] "
            << static_cast<scalar>(
                   std::chrono::duration_cast<std::chrono::nanoseconds>(
                       Clock::now() - last_time_point)
                       .count()) *
                   1e-6
            << " ms" << std::endl;
#endif
  // project to remove motion on the constrained vertices

  // assume that constrained vertices can only be manifold for now
  // first of all, exclude from the solve those constrained vertices being
  // surrounded by all constrained vertices (i.e. with no incident face that is
  // not fully constrained). These vertices don't contribute vorticity because
  //  their incident faces don't contribute vorticity, and they don't desire
  //  displacement correction because they just passively move to wherever they
  //  are prescribed to go. So they don't appear in either columns or rows in
  //  the solve.
#ifdef PRINT_TIMING
  last_time_point = Clock::now();
#endif
  std::vector<int> constrained_vertex_nonconstrained_neighbors(
      m_constrained_vertices.size(),
      0);  // how many non-fully-constrained incident faces each constrained
           // vertex has
  for (size_t i = 0; i < m_constrained_vertices.size(); i++)
    for (size_t j = 0;
         j < mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]].size();
         j++) {
      LosTopos::Vec3st t = mesh().get_triangle(
          mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
      if (!(m_st->vertex_is_any_solid(t[0]) &&
            m_st->vertex_is_any_solid(t[1]) && m_st->vertex_is_any_solid(t[2])))
        constrained_vertex_nonconstrained_neighbors[i]++;
    }

  std::vector<size_t>
      relevant_constrained_vertices;  // constrained vertices that are not fully
                                      // constrained (circulations on those
                                      // vertices don't matter) and are not on
                                      // open boundary (the constraint solve
                                      // there is different and they have
                                      // already been handled by the OB solve
                                      // above)
  for (size_t i = 0; i < m_constrained_vertices.size(); i++) {
    if (constrained_vertex_nonconstrained_neighbors[i] == 0)
      continue;  // ignore constrained vertices who don't have any
                 // non-fully-constrianed incident faces

    std::set<Vec2i, Vec2iComp> rps;
    for (size_t j = 0;
         j < mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]].size();
         j++) {
      LosTopos::Vec3st t = mesh().get_triangle(
          mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
      if (m_st->vertex_is_any_solid(t[0]) && m_st->vertex_is_any_solid(t[1]) &&
          m_st->vertex_is_any_solid(t[2]))
        continue;
      LosTopos::Vec2i l = mesh().get_triangle_label(
          mesh().m_vertex_to_triangle_map[m_constrained_vertices[i]][j]);
      Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
      rps.insert(rp);
    }
    if (rps.size() > 1)
      continue;  // ignore constrained vertices who are incident to more than
                 // one region pair (through non-fully-constrained incident
                 // faces), because the code below can't handle them

    bool ob = false;
    for (size_t j = 0; j < nob; j++)
      if (obv[j] == m_constrained_vertices[i]) ob = true;
    if (ob) continue;  // ignore open boundary vertices

    relevant_constrained_vertices.push_back(i);
  }

  size_t nc = relevant_constrained_vertices.size();

  if (nc > 0) {
    std::vector<Vec2i> constrained_vertex_region_pair(
        nc);  // one region pair for each constrained vertex
    std::vector<Vec3d> constrained_vertex_normal(nc);  // surface normals
    for (size_t i = 0; i < nc; i++) {
      size_t cv = m_constrained_vertices[relevant_constrained_vertices[i]];
      LosTopos::Vec2i l(-1, -1);
      for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[cv].size(); j++) {
        LosTopos::Vec3st t =
            mesh().get_triangle(mesh().m_vertex_to_triangle_map[cv][j]);
        if (m_st->vertex_is_any_solid(t[0]) &&
            m_st->vertex_is_any_solid(t[1]) && m_st->vertex_is_any_solid(t[2]))
          continue;
        l = mesh().get_triangle_label(
            mesh().m_vertex_to_triangle_map[cv][j]);  // grab any triangle on
                                                      // the vertex, because
                                                      // it's manifold.
      }
      assert(l[0] >= 0 && l[1] >= 0);  // this vertex can't be incident to no
                                       // unconstrained triangle.
      constrained_vertex_region_pair[i] =
          (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));

      Vec3d normal(0, 0, 0);
      for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[cv].size(); j++) {
        LosTopos::Vec3st t =
            mesh().get_triangle(mesh().m_vertex_to_triangle_map[cv][j]);
        if (m_st->vertex_is_any_solid(t[0]) &&
            m_st->vertex_is_any_solid(t[1]) && m_st->vertex_is_any_solid(t[2]))
          continue;
        LosTopos::Vec2i ll =
            mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[cv][j]);
        assert(
            (l[0] == ll[0] && l[1] == ll[1]) ||
            (l[0] == ll[1] && l[1] == ll[0]));  // assume constrained vertices
                                                // are all manifold for now
        Vec3d x0 = pos(t[0]);
        Vec3d x1 = pos(t[1]);
        Vec3d x2 = pos(t[2]);
        if (!(m_st->vertex_is_any_solid(t[0]) &&
              m_st->vertex_is_any_solid(t[1]) &&
              m_st->vertex_is_any_solid(t[2])))
          normal += (x1 - x0).cross(x2 - x0) * (l[0] == ll[0] ? 1 : -1);
      }
      assert(normal.norm() != 0);
      constrained_vertex_normal[i] = normal.normalized();
    }

    MatXd A = MatXd::Zero(
        nc,
        nc);  // the transformation from circulations (the additional amount to
              // be added to existing circulations to enforce constraints) on
              // constrained vertices to normal displacements (additional amount
              // as well) on constrained vertices
    for (size_t ii = 0; ii < nc; ii++) {
      // Biot-Savart for constrained vertex i
      size_t i = m_constrained_vertices[relevant_constrained_vertices[ii]];

      Vec3d v(0, 0, 0);
      Vec3d x = pos(i);

      for (size_t jj = 0; jj < nc; jj++) {
        size_t j = m_constrained_vertices[relevant_constrained_vertices[jj]];
        for (size_t k = 0; k < mesh().m_vertex_to_triangle_map[j].size(); k++) {
          LosTopos::Vec3st t =
              mesh().get_triangle(mesh().m_vertex_to_triangle_map[j][k]);
          if (m_st->vertex_is_any_solid(t[0]) &&
              m_st->vertex_is_any_solid(t[1]) &&
              m_st->vertex_is_any_solid(t[2]))
            continue;  // all-solid faces don't contribute vorticity.

          LosTopos::Vec2i l =
              mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[j][k]);
          Vec3d xp = (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;

          Vec3d e01 = pos(t[1]) - pos(t[0]);
          Vec3d e12 = pos(t[2]) - pos(t[1]);
          Vec3d e20 = pos(t[0]) - pos(t[2]);

          Vec3d e_opposite;
          if (t[0] == j)
            e_opposite = e12;
          else if (t[1] == j)
            e_opposite = e20;
          else if (t[2] == j)
            e_opposite = e01;
          assert(t[0] == j || t[1] == j || t[2] == j);

          Vec3d gamma =
              -(e01 * (*m_Gamma)[t[2]].get(l) + e12 * (*m_Gamma)[t[0]].get(l) +
                e20 * (*m_Gamma)[t[1]].get(l));

          Vec3d dx = x - xp;
          //                double dxn = dx.norm();
          double dxn = sqrt(dx.squaredNorm() + m_delta * m_delta);

          v += gamma.cross(dx) / (dxn * dxn * dxn);
          //                v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 -
          //                exp(-dxn / m_delta));

          A(ii, jj) += (l[0] < l[1] ? 1 : -1) *
                       (skewSymmetric(dx) * e_opposite / (dxn * dxn * dxn))
                           .dot(constrained_vertex_normal[ii]);
        }
      }
    }
    A *= dt / (4 * M_PI);

    VecXd rhs = VecXd::Zero(
        nc);  // rhs = the additional normal displacements needed on constrained
              // vertices to correct the current displacement to satisfiy the
              // constraints (in the normal direction)
    for (size_t ii = 0; ii < nc; ii++) {
      size_t i = m_constrained_vertices[relevant_constrained_vertices[ii]];
      rhs(ii) = (m_constrained_positions[relevant_constrained_vertices[ii]] -
                 vc(m_st->pm_newpositions[i]))
                    .dot(constrained_vertex_normal[ii]);
    }

    double lambda =
        (m_st->m_min_edge_length + m_st->m_max_edge_length) / 2 * 0.1;
    VecXd result =
        (A.transpose() * A + lambda * lambda * MatXd::Identity(nc, nc))
            .partialPivLu()
            .solve(A.transpose() *
                   rhs);  // regularized solve to avoid blowing up in presence
                          // of near-dependent constraints

    for (size_t ii = 0; ii < nc; ii++) {
      size_t i = m_constrained_vertices[relevant_constrained_vertices[ii]];
      Vec2i rp = constrained_vertex_region_pair[ii];
      (*m_Gamma)[i].set(rp, (*m_Gamma)[i].get(rp) + result[ii]);
    }
  }
  // recompute the velocity after the constraint projection
  newv = BiotSavartFunc(*this, VecXd::Zero(mesh().nv() * 3));
  for (size_t i = 0; i < mesh().nv(); i++)
    m_st->pm_newpositions[i] =
        m_st->pm_positions[i] + vc(newv.segment<3>(i * 3)) * dt;

  // enforce the constraints exactly, potentially sliding them tangentially.
  // TODO: resample Gammas?
  for (size_t ii = 0; ii < m_constrained_vertices.size(); ii++) {
    size_t i = m_constrained_vertices[ii];
    m_st->pm_newpositions[i] = vc(m_constrained_positions[ii]);
  }
#ifdef PRINT_TIMING
  std::cout << "[enforce the constraints] "
            << static_cast<scalar>(
                   std::chrono::duration_cast<std::chrono::nanoseconds>(
                       Clock::now() - last_time_point)
                       .count()) *
                   1e-6
            << " ms" << std::endl;
#endif
  // move the mesh
  if (counter % 2 != 0) {
#ifdef PRINT_TIMING
    last_time_point = Clock::now();
#endif
    double actual_dt;
    m_st->integrate(dt, actual_dt);
    if (actual_dt != dt)
      std::cout << "Warning: SurfTrack::integrate() failed to step the full "
                   "length of the time step!"
                << std::endl;
#ifdef PRINT_TIMING
    std::cout << "[integrate] "
              << static_cast<scalar>(
                     std::chrono::duration_cast<std::chrono::nanoseconds>(
                         Clock::now() - last_time_point)
                         .count()) *
                     1e-6
              << " ms" << std::endl;
#endif
  }

  return (counter % 2 == 0 ? 0 : dt);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool VS3D::generate_collapsed_position(LosTopos::SurfTrack& st, size_t v0,
                                       size_t v1, LosTopos::Vec3d& pos) {
  if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1)) {
    return false;
  } else if (st.vertex_is_any_solid(v0)) {
    pos = st.pm_positions[v0];
    return true;
  } else if (st.vertex_is_any_solid(v1)) {
    pos = st.pm_positions[v1];
    return true;
  } else {
    pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    return true;
  }
}

bool VS3D::generate_split_position(LosTopos::SurfTrack& st, size_t v0,
                                   size_t v1, LosTopos::Vec3d& pos) {
#ifdef _DEBUG
  std::cout << "solid callback: generate split position: " << v0 << " " << v1
            << " " << (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
            << std::endl;
#endif
  pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
  if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    return false;
  else
    return true;
}

LosTopos::Vec3c VS3D::generate_collapsed_solid_label(
    LosTopos::SurfTrack& st, size_t v0, size_t v1,
    const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1) {
  return LosTopos::Vec3c(1, 1, 1);
}

LosTopos::Vec3c VS3D::generate_split_solid_label(
    LosTopos::SurfTrack& st, size_t v0, size_t v1,
    const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1) {
  if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    return LosTopos::Vec3c(1, 1, 1);
  else if (st.vertex_is_any_solid(v0))
    return LosTopos::Vec3c(0, 0, 0);
  else if (st.vertex_is_any_solid(v1))
    return LosTopos::Vec3c(0, 0, 0);
  else
    return LosTopos::Vec3c(0, 0, 0);
}

bool VS3D::generate_edge_popped_positions(LosTopos::SurfTrack& st, size_t oldv,
                                          const LosTopos::Vec2i& cut,
                                          LosTopos::Vec3d& pos_upper,
                                          LosTopos::Vec3d& pos_lower) {
  return false;
}

bool VS3D::generate_vertex_popped_positions(LosTopos::SurfTrack& st,
                                            size_t oldv, int A, int B,
                                            LosTopos::Vec3d& pos_a,
                                            LosTopos::Vec3d& pos_b) {
  return false;
}

bool VS3D::solid_edge_is_feature(const LosTopos::SurfTrack& st, size_t e) {
  return false;
}

LosTopos::Vec3d VS3D::sampleVelocity(LosTopos::Vec3d& pos) {
  return LosTopos::Vec3d(0, 0, 0);
}

bool VS3D::sampleDirectionalDivergence(const LosTopos::Vec3d& pos,
                                       const LosTopos::Vec3d& dir,
                                       double& output) {
  return false;
}

struct CollapseTempData {
  size_t v0;
  size_t v1;

  Vec3d old_x0;
  Vec3d old_x1;

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void VS3D::pre_collapse(const LosTopos::SurfTrack& st, size_t e, void** data) {
  CollapseTempData* td = new CollapseTempData;
  td->v0 = st.m_mesh.m_edges[e][0];
  td->v1 = st.m_mesh.m_edges[e][1];

  td->old_x0 = vc(st.pm_positions[td->v0]);
  td->old_x1 = vc(st.pm_positions[td->v1]);

  td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
    td->v0_incident_region_pairs(l[0], l[1]) =
        td->v0_incident_region_pairs(l[1], l[0]) = true;
  }
  td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
    td->v1_incident_region_pairs(l[0], l[1]) =
        td->v1_incident_region_pairs(l[1], l[0]) = true;
  }

  *data = (void*)td;
#ifdef _DEBUG
  std::cout << "pre collapse: " << e << ": " << td->v0 << " " << td->v1
            << std::endl;
#endif
}

void VS3D::post_collapse(const LosTopos::SurfTrack& st, size_t e,
                         size_t merged_vertex, void* data) {
  CollapseTempData* td = (CollapseTempData*)data;
#ifdef _DEBUG
  std::cout << "post collapse: " << e << ": " << td->v0 << " " << td->v1
            << " => " << merged_vertex << std::endl;
#endif
  assert((st.m_mesh.vertex_is_deleted(td->v0) && merged_vertex == td->v1) ||
         (st.m_mesh.vertex_is_deleted(td->v1) && merged_vertex == td->v0));

  Vec3d merged_x = vc(st.pm_positions[merged_vertex]);
  double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) /
             (td->old_x1 - td->old_x0).squaredNorm();

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>
      merged_vertex_incident_region_pairs;
  merged_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0;
       i < st.m_mesh.m_vertex_to_triangle_map[merged_vertex].size(); i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[merged_vertex][i]);
    merged_vertex_incident_region_pairs(l[0], l[1]) =
        merged_vertex_incident_region_pairs(l[1], l[0]) = true;
  }

  // for region pairs not existing in original v0 and v1, we cannot simply
  // assume 0 circulation because the neighbors that do have those region pairs
  //  may have accumulated some amount of circulations, and filling in 0 will
  //  cause vorticity spikes. The information that should be filled in here can
  //  only come from neighbors.
  GammaType newGamma(m_nregion);
  newGamma.setZero();
  for (int i = 0; i < m_nregion; i++) {
    for (int j = i + 1; j < m_nregion; j++) {
      if (merged_vertex_incident_region_pairs(i, j)) {
        if (!td->v0_incident_region_pairs(i, j) &&
            !td->v1_incident_region_pairs(i, j)) {
          double neighborhood_mean = 0;
          int neighborhood_counter = 0;
          for (size_t k = 0;
               k < st.m_mesh.m_vertex_to_edge_map[merged_vertex].size(); k++) {
            LosTopos::Vec2st e =
                st.m_mesh
                    .m_edges[st.m_mesh.m_vertex_to_edge_map[merged_vertex][k]];
            size_t vother = (e[0] == merged_vertex ? e[1] : e[0]);
            bool incident_to_this_region_pair = false;
            for (size_t l = 0;
                 l < st.m_mesh.m_vertex_to_triangle_map[vother].size(); l++) {
              LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                  st.m_mesh.m_vertex_to_triangle_map[vother][l]);
              if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i)) {
                incident_to_this_region_pair = true;
                break;
              }
            }

            if (incident_to_this_region_pair) {
              neighborhood_mean += (*m_Gamma)[vother].get(i, j);
              neighborhood_counter++;
            }
          }
          if (neighborhood_counter != 0)
            neighborhood_mean /= neighborhood_counter;

          newGamma.set(i, j, neighborhood_mean);
        } else if (!td->v0_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
        } else if (!td->v1_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
        } else {
          newGamma.set(i, j,
                       (*m_Gamma)[td->v0].get(i, j) * (1 - s) +
                           (*m_Gamma)[td->v1].get(i, j) * s);
        }
      }
    }
  }
  (*m_Gamma)[merged_vertex] = newGamma;
}

struct SplitTempData {
  size_t v0;
  size_t v1;

  Vec3d old_x0;
  Vec3d old_x1;

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void VS3D::pre_split(const LosTopos::SurfTrack& st, size_t e, void** data) {
  SplitTempData* td = new SplitTempData;
  td->v0 = st.m_mesh.m_edges[e][0];
  td->v1 = st.m_mesh.m_edges[e][1];

  td->old_x0 = vc(st.pm_positions[td->v0]);
  td->old_x1 = vc(st.pm_positions[td->v1]);

  td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
    td->v0_incident_region_pairs(l[0], l[1]) =
        td->v0_incident_region_pairs(l[1], l[0]) = true;
  }
  td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
    td->v1_incident_region_pairs(l[0], l[1]) =
        td->v1_incident_region_pairs(l[1], l[0]) = true;
  }

  *data = (void*)td;
#ifdef _DEBUG
  std::cout << "pre split: " << e << ": " << td->v0 << " " << td->v1
            << std::endl;
#endif
}

void VS3D::post_split(const LosTopos::SurfTrack& st, size_t e,
                      size_t new_vertex, void* data) {
  SplitTempData* td = (SplitTempData*)data;
#ifdef _DEBUG
  std::cout << "post split: " << e << ": " << td->v0 << " " << td->v1 << " => "
            << new_vertex << std::endl;
#endif
  Vec3d midpoint_x = vc(st.pm_positions[new_vertex]);
  double s = (midpoint_x - td->old_x0).dot(td->old_x1 - td->old_x0) /
             (td->old_x1 - td->old_x0).squaredNorm();

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>
      new_vertex_incident_region_pairs;
  new_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[new_vertex].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[new_vertex][i]);
    new_vertex_incident_region_pairs(l[0], l[1]) =
        new_vertex_incident_region_pairs(l[1], l[0]) = true;
  }

  // for region pairs not existing in original v0 and v1, we cannot simply
  // assume 0 circulation because the neighbors that do have those region pairs
  //  may have accumulated some amount of circulations, and filling in 0 will
  //  cause vorticity spikes. The information that should be filled in here can
  //  only come from neighbors.
  GammaType newGamma(m_nregion);
  newGamma.setZero();
  for (int i = 0; i < m_nregion; i++) {
    for (int j = i + 1; j < m_nregion; j++) {
      if (new_vertex_incident_region_pairs(i, j)) {
        if (!td->v0_incident_region_pairs(i, j) &&
            !td->v1_incident_region_pairs(i, j)) {
          double neighborhood_mean = 0;
          int neighborhood_counter = 0;
          for (size_t k = 0;
               k < st.m_mesh.m_vertex_to_edge_map[new_vertex].size(); k++) {
            LosTopos::Vec2st e =
                st.m_mesh
                    .m_edges[st.m_mesh.m_vertex_to_edge_map[new_vertex][k]];
            size_t vother = (e[0] == new_vertex ? e[1] : e[0]);

            bool incident_to_this_region_pair = false;
            for (size_t l = 0;
                 l < st.m_mesh.m_vertex_to_triangle_map[vother].size(); l++) {
              LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                  st.m_mesh.m_vertex_to_triangle_map[vother][l]);
              if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i)) {
                incident_to_this_region_pair = true;
                break;
              }
            }

            if (incident_to_this_region_pair) {
              neighborhood_mean += (*m_Gamma)[vother].get(i, j);
              neighborhood_counter++;
            }
          }
          if (neighborhood_counter != 0)
            neighborhood_mean /= neighborhood_counter;

          newGamma.set(i, j, neighborhood_mean);
        } else if (!td->v0_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
        } else if (!td->v1_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
        } else {
          newGamma.set(i, j,
                       (*m_Gamma)[td->v0].get(i, j) * (1 - s) +
                           (*m_Gamma)[td->v1].get(i, j) * s);
        }
      }
    }
  }
  (*m_Gamma)[new_vertex] = newGamma;
}

void VS3D::pre_flip(const LosTopos::SurfTrack& st, size_t e, void** data) {}

void VS3D::post_flip(const LosTopos::SurfTrack& st, size_t e, void* data) {}

struct T1TempData {
  std::vector<size_t> neighbor_verts;
  std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> >
      neighbor_region_pairs;
};

void VS3D::pre_t1(const LosTopos::SurfTrack& st, size_t v, void** data) {
  T1TempData* td = new T1TempData;

  for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v].size(); i++) {
    LosTopos::Vec2st e =
        st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[v][i]];
    size_t vother = (e[0] == v ? e[1] : e[0]);
    td->neighbor_verts.push_back(vother);

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> rps;
    rps.setZero(m_nregion, m_nregion);
    for (size_t j = 0; j < st.m_mesh.m_vertex_to_triangle_map[vother].size();
         j++) {
      LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
          st.m_mesh.m_vertex_to_triangle_map[vother][j]);
      rps(l[0], l[1]) = rps(l[1], l[0]) = true;
    }
    td->neighbor_region_pairs.push_back(rps);
  }

  *data = (void*)td;
}

void VS3D::post_t1(const LosTopos::SurfTrack& st, size_t v, size_t a, size_t b,
                   void* data) {
#ifdef _DEBUG
  std::cout << "v = " << v << " -> " << a << " " << b << std::endl;
#endif
  T1TempData* td = (T1TempData*)data;

  (*m_Gamma)[a] = (*m_Gamma)[v];
  (*m_Gamma)[b] = (*m_Gamma)[v];
#ifdef _DEBUG
  std::cout << "v Gammas: " << std::endl << (*m_Gamma)[v].values << std::endl;
#endif
  for (size_t i = 0; i < td->neighbor_verts.size(); i++) {
    size_t vother = td->neighbor_verts[i];

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> rps;
    rps.setZero(m_nregion, m_nregion);
    for (size_t j = 0; j < st.m_mesh.m_vertex_to_triangle_map[vother].size();
         j++) {
      LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
          st.m_mesh.m_vertex_to_triangle_map[vother][j]);
      rps(l[0], l[1]) = rps(l[1], l[0]) = true;
    }
#ifdef _DEBUG
    std::cout << vother << std::endl
              << td->neighbor_region_pairs[i] << std::endl
              << "-> " << std::endl
              << rps << std::endl;
#endif
    for (int j = 0; j < m_nregion; j++)
      for (int k = j + 1; k < m_nregion; k++)
        if (rps(j, k) && !td->neighbor_region_pairs[i](j, k))
          (*m_Gamma)[vother].set(j, k, (*m_Gamma)[v].get(j, k));
#ifdef _DEBUG
    std::cout << "Gammas: " << std::endl
              << (*m_Gamma)[vother].values << std::endl;
#endif
  }
}

struct FaceSplitTempData {
  size_t v0;
  size_t v1;
  size_t v2;
};

void VS3D::pre_facesplit(const LosTopos::SurfTrack& st, size_t f, void** data) {
  FaceSplitTempData* td = new FaceSplitTempData;

  td->v0 = st.m_mesh.get_triangle(f)[0];
  td->v1 = st.m_mesh.get_triangle(f)[1];
  td->v2 = st.m_mesh.get_triangle(f)[2];

  *data = (void*)td;
}

void VS3D::post_facesplit(const LosTopos::SurfTrack& st, size_t f,
                          size_t new_vertex, void* data) {
  FaceSplitTempData* td = (FaceSplitTempData*)data;

  (*m_Gamma)[new_vertex] = GammaType(m_nregion);
  (*m_Gamma)[new_vertex].values =
      ((*m_Gamma)[td->v0].values + (*m_Gamma)[td->v1].values +
       (*m_Gamma)[td->v2].values) /
      3;
}

struct SnapTempData {
  size_t v0;
  size_t v1;

  Vec3d old_x0;
  Vec3d old_x1;

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v0_incident_region_pairs;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> v1_incident_region_pairs;
};

void VS3D::pre_snap(const LosTopos::SurfTrack& st, size_t v0, size_t v1,
                    void** data) {
  SnapTempData* td = new SnapTempData;
  td->v0 = v0;
  td->v1 = v1;

  td->old_x0 = vc(st.pm_positions[v0]);
  td->old_x1 = vc(st.pm_positions[v1]);

  td->v0_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v0].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v0][i]);
    td->v0_incident_region_pairs(l[0], l[1]) =
        td->v0_incident_region_pairs(l[1], l[0]) = true;
  }
  td->v1_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[td->v1].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[td->v1][i]);
    td->v1_incident_region_pairs(l[0], l[1]) =
        td->v1_incident_region_pairs(l[1], l[0]) = true;
  }

  *data = (void*)td;
#ifdef _DEBUG
  std::cout << "pre snap: " << v0 << " " << v1 << std::endl;
#endif
}

void VS3D::post_snap(const LosTopos::SurfTrack& st, size_t v_kept,
                     size_t v_deleted, void* data) {
  SnapTempData* td = (SnapTempData*)data;
#ifdef _DEBUG
  std::cout << "post snap: " << td->v0 << " " << td->v1 << " => " << v_kept
            << std::endl;
#endif
  assert((td->v0 == v_kept && td->v1 == v_deleted) ||
         (td->v1 == v_kept && td->v0 == v_deleted));
  assert(v_kept != v_deleted);
  assert(st.m_mesh.vertex_is_deleted(v_deleted));
  assert(!st.m_mesh.vertex_is_deleted(v_kept));

  Vec3d merged_x = vc(st.pm_positions[v_kept]);
  double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) /
             (td->old_x1 - td->old_x0).squaredNorm();

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>
      merged_vertex_incident_region_pairs;
  merged_vertex_incident_region_pairs.setZero(m_nregion, m_nregion);
  for (size_t i = 0; i < st.m_mesh.m_vertex_to_triangle_map[v_kept].size();
       i++) {
    LosTopos::Vec2i l = st.m_mesh.get_triangle_label(
        st.m_mesh.m_vertex_to_triangle_map[v_kept][i]);
#ifdef _DEBUG
    std::cout << "triangle " << st.m_mesh.m_vertex_to_triangle_map[v_kept][i]
              << " label = " << l << std::endl;
#endif
    merged_vertex_incident_region_pairs(l[0], l[1]) =
        merged_vertex_incident_region_pairs(l[1], l[0]) = true;
  }
#ifdef _DEBUG
  std::cout << "v0 incident region pairs: " << std::endl
            << td->v0_incident_region_pairs << std::endl;
  std::cout << "v1 incident region pairs: " << std::endl
            << td->v1_incident_region_pairs << std::endl;
  std::cout << "merged vertex incident region pairs: " << std::endl
            << merged_vertex_incident_region_pairs << std::endl;
#endif

  // for region pairs not existing in original v0 and v1, we cannot simply
  // assume 0 circulation because the neighbors that do have those region pairs
  //  may have accumulated some amount of circulations, and filling in 0 will
  //  cause vorticity spikes. The information that should be filled in here can
  //  only come from neighbors.
  GammaType newGamma(m_nregion);
  newGamma.setZero();
  for (int i = 0; i < m_nregion; i++) {
    for (int j = i + 1; j < m_nregion; j++) {
      if (merged_vertex_incident_region_pairs(i, j)) {
        if (!td->v0_incident_region_pairs(i, j) &&
            !td->v1_incident_region_pairs(i, j)) {
#ifdef _DEBUG
          std::cout << "region pair " << i << " " << j
                    << " is being computed from 1-ring neighbors" << std::endl;
#endif
          double neighborhood_mean = 0;
          int neighborhood_counter = 0;
          for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[v_kept].size();
               k++) {
            LosTopos::Vec2st e =
                st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[v_kept][k]];
            size_t vother = (e[0] == v_kept ? e[1] : e[0]);

            bool incident_to_this_region_pair = false;
            for (size_t l = 0;
                 l < st.m_mesh.m_vertex_to_triangle_map[vother].size(); l++) {
              LosTopos::Vec2i ll = st.m_mesh.get_triangle_label(
                  st.m_mesh.m_vertex_to_triangle_map[vother][l]);
              if ((ll[0] == i && ll[1] == j) || (ll[0] == j && ll[1] == i)) {
                incident_to_this_region_pair = true;
                break;
              }
            }
#ifdef _DEBUG
            std::cout << "vother = " << vother
                      << " incident = " << incident_to_this_region_pair
                      << std::endl;
#endif
            if (incident_to_this_region_pair) {
              neighborhood_mean += (*m_Gamma)[vother].get(i, j);
              neighborhood_counter++;
            }
          }
          if (neighborhood_counter != 0)
            neighborhood_mean /= neighborhood_counter;

          newGamma.set(i, j, neighborhood_mean);
        } else if (!td->v0_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v1].get(i, j));
        } else if (!td->v1_incident_region_pairs(i, j)) {
          newGamma.set(i, j, (*m_Gamma)[td->v0].get(i, j));
        } else {
          newGamma.set(i, j,
                       (*m_Gamma)[td->v0].get(i, j) * (1 - s) +
                           (*m_Gamma)[td->v1].get(i, j) * s);
        }
      }
    }
  }
#ifdef _DEBUG
  std::cout << "v0 Gamma = " << std::endl
            << (*m_Gamma)[td->v0].values << std::endl;
  std::cout << "v1 Gamma = " << std::endl
            << (*m_Gamma)[td->v1].values << std::endl;
  std::cout << "new Gamma = " << std::endl << newGamma.values << std::endl;
#endif
  (*m_Gamma)[v_kept] = newGamma;
}

void VS3D::accumulateGradU(VectorXs& F, const VectorXs& dx,
                           const VectorXs& dv) {
  const int ndof = m_constrained_positions.size() * 3;

  // Accumulate all energy gradients
  if (dx.size() == 0) {
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addGradEToTotal(
          Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), F);
  } else {
    VectorXs nx =
        Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
    VectorXs nv =
        Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addGradEToTotal(
          nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
          F);
  }
}

void VS3D::accumulateddUdxdx(TripletXs& A, const VectorXs& dx,
                             const VectorXs& dv) {
  const int ndof = m_constrained_positions.size() * 3;
  if (dx.size() == 0) {
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addHessXToTotal(
          Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), A);
  } else {
    VectorXs nx =
        Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
    VectorXs nv =
        Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addHessXToTotal(
          nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
          A);
  }
}

// Kind of a misnomer.
void VS3D::accumulateddUdxdv(TripletXs& A, const VectorXs& dx,
                             const VectorXs& dv) {
  const int ndof = m_constrained_positions.size() * 3;
  if (dx.size() == 0) {
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addHessVToTotal(
          Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), A);
  } else {
    VectorXs nx =
        Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
    VectorXs nv =
        Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->addHessVToTotal(
          nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
          A);
  }
}

void VS3D::preCompute(const VectorXs& dx, const VectorXs& dv,
                      const scalar& dt) {
  const int ndof = m_constrained_positions.size() * 3;
  if (dx.size() == 0) {
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->preCompute(
          Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof),
          Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof), dt);
  } else {
    VectorXs nx =
        Eigen::Map<VectorXs>((double*)&m_constrained_positions[0], ndof) + dx;
    VectorXs nv =
        Eigen::Map<VectorXs>((double*)&m_constrained_velocities[0], ndof) + dv;
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
      m_forces[i]->preCompute(
          nx, nv, Eigen::Map<VectorXs>((double*)&m_constrained_mass[0], ndof),
          dt);
  }
}

void VS3D::pre_smoothing(const LosTopos::SurfTrack& st, void** data) {}

void VS3D::post_smoothing(const LosTopos::SurfTrack& st, void* data) {}
