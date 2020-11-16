// ---------------------------------------------------------
//
//  surftrack.cpp
//  Tyson Brochu 2008
//  Christopher Batty, Fang Da 2014
//
//  Implementation of the SurfTrack class: a dynamic mesh with
//  topological changes and mesh maintenance operations.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <array3.h>
#include <broadphase.h>
#include <ccd_wrapper.h>
#include <collisionpipeline.h>
#include <collisionqueries.h>
#include <edgeflipper.h>
#include <impactzonesolver.h>
#include <lapack_wrapper.h>
#include <nondestructivetrimesh.h>
#include <runstats.h>
#include <stdio.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>
#include <vec.h>
#include <wallclocktime.h>

#include <cassert>
#include <queue>
#include <vector>

// ---------------------------------------------------------
//  Global externs
// ---------------------------------------------------------
namespace LosTopos {

double G_EIGENVALUE_RANK_RATIO = 0.03;

extern RunStats g_stats;

// ---------------------------------------------------------
//  Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Default initialization parameters
///
// ---------------------------------------------------------

SurfTrackInitializationParameters::SurfTrackInitializationParameters()
    : m_proximity_epsilon(1e-4),
      m_friction_coefficient(0.0),
      m_min_triangle_area(1e-7),
      m_t1_transition_enabled(false),
      m_velocity_field_callback(NULL),
      m_improve_collision_epsilon(2e-6),
      m_use_fraction(false),
      m_min_edge_length(
          UNINITIALIZED_DOUBLE),  // <- Don't allow instantiation without
                                  // setting these parameters
      m_max_edge_length(UNINITIALIZED_DOUBLE),    // <-
      m_max_volume_change(UNINITIALIZED_DOUBLE),  // <-
      m_min_triangle_angle(2.0),
      m_max_triangle_angle(178.0),
      m_large_triangle_angle_to_split(135.0),
      m_use_curvature_when_splitting(false),
      m_use_curvature_when_collapsing(false),
      m_min_curvature_multiplier(1.0),
      m_max_curvature_multiplier(1.0),
      m_allow_vertex_movement_during_collapse(true),
      m_perform_smoothing(true),
      m_merge_proximity_epsilon(1e-5),
      m_subdivision_scheme(NULL),
      m_collision_safety(true),
      m_allow_topology_changes(true),
      m_allow_non_manifold(true),
      m_perform_improvement(true),
      m_remesh_boundaries(true),
      m_pull_apart_distance(0.1),
      m_verbose(false) {}

// ---------------------------------------------------------
///
/// Create a SurfTrack object from a set of vertices and triangles using the
/// specified parameters
///
// ---------------------------------------------------------

#ifdef _MSC_VER
#pragma warning(disable : 4355)
#endif

SurfTrack::SurfTrack(
    const std::vector<Vec3d>& vs, const std::vector<Vec3st>& ts,
    const std::vector<Vec2i>& labels, const std::vector<Vec3d>& masses,
    const SurfTrackInitializationParameters& initial_parameters)
    :

      DynamicSurface(
          vs, ts, labels, masses, initial_parameters.m_proximity_epsilon,
          initial_parameters.m_friction_coefficient,
          initial_parameters.m_collision_safety, initial_parameters.m_verbose),

      m_collapser(*this, initial_parameters.m_use_curvature_when_collapsing,
                  initial_parameters.m_remesh_boundaries,
                  initial_parameters.m_min_curvature_multiplier),
      m_splitter(*this, initial_parameters.m_use_curvature_when_splitting,
                 initial_parameters.m_remesh_boundaries,
                 initial_parameters.m_max_curvature_multiplier),
      m_flipper(*this),
      m_smoother(*this),
      m_merger(*this),
      m_pincher(*this),
      m_cutter(*this),
      m_snapper(*this, initial_parameters.m_use_curvature_when_splitting,
                initial_parameters.m_remesh_boundaries,
                initial_parameters.m_max_curvature_multiplier),
      m_t1transition(*this, initial_parameters.m_velocity_field_callback,
                     initial_parameters.m_remesh_boundaries),
      m_t1_transition_enabled(initial_parameters.m_t1_transition_enabled),
      m_improve_collision_epsilon(
          initial_parameters.m_improve_collision_epsilon),
      m_max_volume_change(UNINITIALIZED_DOUBLE),
      m_min_edge_length(UNINITIALIZED_DOUBLE),
      m_max_edge_length(UNINITIALIZED_DOUBLE),
      m_merge_proximity_epsilon(initial_parameters.m_merge_proximity_epsilon),
      m_min_triangle_area(initial_parameters.m_min_triangle_area),
      m_min_triangle_angle(initial_parameters.m_min_triangle_angle),
      m_max_triangle_angle(initial_parameters.m_max_triangle_angle),
      m_hard_min_edge_len(0.05 * initial_parameters.m_min_edge_length),
      m_hard_max_edge_len(10.0 * initial_parameters.m_max_edge_length),
      m_large_triangle_angle_to_split(
          initial_parameters.m_large_triangle_angle_to_split),
      m_subdivision_scheme(initial_parameters.m_subdivision_scheme),
      should_delete_subdivision_scheme_object(
          m_subdivision_scheme == NULL ? true : false),
      m_dirty_triangles(0),
      m_allow_topology_changes(initial_parameters.m_allow_topology_changes),
      m_allow_non_manifold(initial_parameters.m_allow_non_manifold),
      m_perform_improvement(initial_parameters.m_perform_improvement),
      m_remesh_boundaries(initial_parameters.m_remesh_boundaries),
      m_aggressive_mode(false),
      m_allow_vertex_movement_during_collapse(
          initial_parameters.m_allow_vertex_movement_during_collapse),
      m_perform_smoothing(initial_parameters.m_perform_smoothing),
      m_mesheventcallback(NULL),
      m_solid_vertices_callback(NULL),
      m_vertex_change_history(),
      m_triangle_change_history(),
      m_defragged_triangle_map(),
      m_defragged_vertex_map() {
  if (m_verbose) {
    std::cout << " ======== SurfTrack ======== " << std::endl;
    std::cout << "m_allow_topology_changes: " << m_allow_topology_changes
              << std::endl;
    std::cout << "m_perform_improvement: " << m_perform_improvement
              << std::endl;
    std::cout << "m_min_triangle_area: " << m_min_triangle_area << std::endl;
    std::cout << "initial_parameters.m_use_fraction: "
              << initial_parameters.m_use_fraction << std::endl;
  }

  if (m_collision_safety) {
    rebuild_static_broad_phase();
  }

  assert(initial_parameters.m_min_edge_length != UNINITIALIZED_DOUBLE);
  assert(initial_parameters.m_max_edge_length != UNINITIALIZED_DOUBLE);
  assert(initial_parameters.m_max_volume_change != UNINITIALIZED_DOUBLE);

  if (initial_parameters.m_use_fraction) {
    double avg_length = DynamicSurface::get_average_non_solid_edge_length();
    m_collapser.m_min_edge_length =
        initial_parameters.m_min_edge_length * avg_length;
    m_collapser.m_max_edge_length =
        initial_parameters.m_max_edge_length * avg_length;

    m_splitter.m_max_edge_length =
        initial_parameters.m_max_edge_length * avg_length;
    m_splitter.m_min_edge_length =
        initial_parameters.m_min_edge_length * avg_length;

    m_min_edge_length = initial_parameters.m_min_edge_length * avg_length;
    m_max_edge_length = initial_parameters.m_max_edge_length * avg_length;
    m_max_volume_change = initial_parameters.m_max_volume_change * avg_length *
                          avg_length * avg_length;

    m_t1transition.m_pull_apart_distance =
        avg_length * initial_parameters.m_pull_apart_distance * 2;
    m_collapser.m_t1_pull_apart_distance =
        avg_length * initial_parameters.m_pull_apart_distance * 2;

  } else {
    m_collapser.m_min_edge_length = initial_parameters.m_min_edge_length;
    m_collapser.m_max_edge_length = initial_parameters.m_max_edge_length;

    m_splitter.m_max_edge_length = initial_parameters.m_max_edge_length;
    m_splitter.m_min_edge_length = initial_parameters.m_min_edge_length;

    m_min_edge_length = initial_parameters.m_min_edge_length;
    m_max_edge_length = initial_parameters.m_max_edge_length;
    m_max_volume_change = initial_parameters.m_max_volume_change;

    m_t1transition.m_pull_apart_distance =
        initial_parameters.m_pull_apart_distance * 2;
    m_collapser.m_t1_pull_apart_distance =
        initial_parameters.m_pull_apart_distance * 2;
  }

  if (m_verbose) {
    std::cout << "m_min_edge_length: " << m_min_edge_length << std::endl;
    std::cout << "m_max_edge_length: " << m_max_edge_length << std::endl;
    std::cout << "m_max_volume_change: " << m_max_volume_change << std::endl;
  }

  if (m_subdivision_scheme == NULL) {
    m_subdivision_scheme = new MidpointScheme();
    should_delete_subdivision_scheme_object = true;
  } else {
    should_delete_subdivision_scheme_object = false;
  }

  if (false == m_allow_topology_changes) {
    m_allow_non_manifold = false;
  }

  assert_no_bad_labels();
}

// ---------------------------------------------------------
///
/// Destructor.  Deallocates the subdivision scheme object if we created one.
///
// ---------------------------------------------------------

SurfTrack::~SurfTrack() {
  if (should_delete_subdivision_scheme_object) {
    delete m_subdivision_scheme;
  }
}

// ---------------------------------------------------------
///
/// Add a triangle to the surface.  Update the underlying TriMesh and
/// acceleration grid.
///
// ---------------------------------------------------------

size_t SurfTrack::add_triangle(const Vec3st& t, const Vec2i& label) {
  size_t new_triangle_index = m_mesh.nondestructive_add_triangle(t, label);

  assert(t[0] < get_num_vertices());
  assert(t[1] < get_num_vertices());
  assert(t[2] < get_num_vertices());

  if (m_collision_safety) {
    // Add to the triangle grid
    Vec3d low, high;
    triangle_static_bounds(new_triangle_index, low, high);
    m_broad_phase->add_triangle(new_triangle_index, low, high,
                                triangle_is_all_solid(new_triangle_index));

    // Add edges to grid as well
    size_t new_edge_index = m_mesh.get_edge_index(t[0], t[1]);
    assert(new_edge_index != m_mesh.m_edges.size());
    edge_static_bounds(new_edge_index, low, high);
    m_broad_phase->add_edge(new_edge_index, low, high,
                            edge_is_all_solid(new_edge_index));

    new_edge_index = m_mesh.get_edge_index(t[1], t[2]);
    assert(new_edge_index != m_mesh.m_edges.size());
    edge_static_bounds(new_edge_index, low, high);
    m_broad_phase->add_edge(new_edge_index, low, high,
                            edge_is_all_solid(new_edge_index));

    new_edge_index = m_mesh.get_edge_index(t[2], t[0]);
    assert(new_edge_index != m_mesh.m_edges.size());
    edge_static_bounds(new_edge_index, low, high);
    m_broad_phase->add_edge(new_edge_index, low, high,
                            edge_is_all_solid(new_edge_index));
  }

  m_triangle_change_history.push_back(TriangleUpdateEvent(
      TriangleUpdateEvent::TRIANGLE_ADD, new_triangle_index, t));

  return new_triangle_index;
}

// ---------------------------------------------------------
///
/// Remove a triangle from the surface.  Update the underlying TriMesh and
/// acceleration grid.
///
// ---------------------------------------------------------

void SurfTrack::remove_triangle(size_t t) {
  m_mesh.nondestructive_remove_triangle(t);
  if (m_collision_safety) {
    m_broad_phase->remove_triangle(t);
  }

  m_triangle_change_history.push_back(
      TriangleUpdateEvent(TriangleUpdateEvent::TRIANGLE_REMOVE, t, Vec3st(0)));
}

// ---------------------------------------------------------
///
/// Efficiently renumber a triangle (replace its vertices) for defragging the
/// mesh. Assume that the underlying geometry doesn't change, so we can keep the
/// same broadphase data unchanged.
///
// ---------------------------------------------------------

void SurfTrack::renumber_triangle(size_t tri, const Vec3st& verts) {
  assert(verts[0] < get_num_vertices());
  assert(verts[1] < get_num_vertices());
  assert(verts[2] < get_num_vertices());

  m_mesh.nondestructive_renumber_triangle(tri, verts);

  // Assume the geometry doesn't change, so we leave all the broad phase
  // collision data alone...

  // Kind of a hack for the history, for now.
  m_triangle_change_history.push_back(TriangleUpdateEvent(
      TriangleUpdateEvent::TRIANGLE_REMOVE, tri, Vec3st(0)));
  m_triangle_change_history.push_back(
      TriangleUpdateEvent(TriangleUpdateEvent::TRIANGLE_ADD, tri, verts));
}

// ---------------------------------------------------------
///
/// Add a vertex to the surface.  Update the acceleration grid.
///
// ---------------------------------------------------------

size_t SurfTrack::add_vertex(const Vec3d& new_vertex_position,
                             const Vec3d& new_vertex_mass) {
  size_t new_vertex_index = m_mesh.nondestructive_add_vertex();

  if (new_vertex_index > get_num_vertices() - 1) {
    pm_positions.resize(new_vertex_index + 1);
    pm_newpositions.resize(new_vertex_index + 1);
    m_masses.resize(new_vertex_index + 1);

    pm_velocities.resize(new_vertex_index + 1);
  }

  pm_positions[new_vertex_index] = new_vertex_position;
  pm_newpositions[new_vertex_index] = new_vertex_position;
  m_masses[new_vertex_index] = new_vertex_mass;

  pm_velocities[new_vertex_index] = Vec3d(0);

  ////////////////////////////////////////////////////////////

  if (m_collision_safety) {
    m_broad_phase->add_vertex(new_vertex_index, get_position(new_vertex_index),
                              get_position(new_vertex_index),
                              vertex_is_all_solid(new_vertex_index));
  }

  return new_vertex_index;
}

// ---------------------------------------------------------
///
/// Remove a vertex from the surface.  Update the acceleration grid.
///
// ---------------------------------------------------------

void SurfTrack::remove_vertex(size_t vertex_index) {
  m_mesh.nondestructive_remove_vertex(vertex_index);

  if (m_collision_safety) {
    m_broad_phase->remove_vertex(vertex_index);
  }

  m_vertex_change_history.push_back(VertexUpdateEvent(
      VertexUpdateEvent::VERTEX_REMOVE, vertex_index, Vec2st(0, 0)));
}

// ---------------------------------------------------------
///
/// Remove deleted vertices and triangles from the mesh data structures
///
// ---------------------------------------------------------

void SurfTrack::defrag_mesh() {
  assert(!"depcrated; use defrag_mesh_from_scratch() instead.");

  //
  // First clear deleted vertices from the data structures
  //

  double start_time = get_time_in_seconds();

  // do a quick pass through to see if any vertices have been deleted
  bool any_deleted = false;
  for (size_t i = 0; i < get_num_vertices(); ++i) {
    if (m_mesh.vertex_is_deleted(i)) {
      any_deleted = true;
      break;
    }
  }

  // resize/allocate up front rather than via push_backs
  m_defragged_vertex_map.resize(get_num_vertices());

  if (!any_deleted) {
    for (size_t i = 0; i < get_num_vertices(); ++i) {
      m_defragged_vertex_map[i] = Vec2st(i, i);
    }

    double end_time = get_time_in_seconds();
    g_stats.add_to_double("total_clear_deleted_vertices_time",
                          end_time - start_time);

  } else {
    // Note: We could rebuild the mesh from scratch, rather than adding/removing
    // triangles, however this function is not a major computational bottleneck.

    size_t j = 0;

    std::vector<Vec3st> new_tris = m_mesh.get_triangles();

    for (size_t i = 0; i < get_num_vertices(); ++i) {
      if (!m_mesh.vertex_is_deleted(i)) {
        pm_positions[j] = pm_positions[i];
        pm_newpositions[j] = pm_newpositions[i];
        m_masses[j] = m_masses[i];

        m_defragged_vertex_map[i] = Vec2st(i, j);

        // Now rewire the triangles containing vertex i

        // copy this, since we'll be changing the original as we go
        std::vector<size_t> inc_tris = m_mesh.m_vertex_to_triangle_map[i];

        for (size_t t = 0; t < inc_tris.size(); ++t) {
          Vec3st triangle = m_mesh.get_triangle(inc_tris[t]);
          Vec2i tri_label = m_mesh.get_triangle_label(inc_tris[t]);

          assert(triangle[0] == i || triangle[1] == i || triangle[2] == i);
          if (triangle[0] == i) {
            triangle[0] = j;
          }
          if (triangle[1] == i) {
            triangle[1] = j;
          }
          if (triangle[2] == i) {
            triangle[2] = j;
          }

          // remove_triangle(inc_tris[t]);       // mark the triangle deleted
          // add_triangle(triangle, tri_label);  // add the updated triangle
          renumber_triangle(inc_tris[t], triangle);
        }

        ++j;
      }
    }

    pm_positions.resize(j);
    pm_newpositions.resize(j);
    m_masses.resize(j);
  }

  double end_time = get_time_in_seconds();

  g_stats.add_to_double("total_clear_deleted_vertices_time",
                        end_time - start_time);

  //
  // Now clear deleted triangles from the mesh
  //

  m_mesh.set_num_vertices(get_num_vertices());
  m_mesh.clear_deleted_triangles(&m_defragged_triangle_map);

  if (m_collision_safety) {
    rebuild_continuous_broad_phase();
  }
}

void SurfTrack::defrag_mesh_from_scratch(
    std::vector<size_t>& vertices_to_be_mapped) {
  defrag_mesh_from_scratch_manual(vertices_to_be_mapped);
}

void SurfTrack::defrag_mesh_from_scratch_manual(
    std::vector<size_t>& vertices_to_be_mapped) {
  double start_time = get_time_in_seconds();

  // defragment vertices
  std::vector<int> vm(get_num_vertices(), -1);
  size_t j = 0;
  for (size_t i = 0; i < get_num_vertices(); i++)
    if (!m_mesh.vertex_is_deleted(i)) vm[i] = j++;

  for (size_t i = 0; i < vertices_to_be_mapped.size(); i++)
    vertices_to_be_mapped[i] = vm[vertices_to_be_mapped[i]];

  for (size_t i = 0; i < m_mesh.m_vds.size(); i++) {
    m_mesh.m_vds[i]->compress(vm);
    m_mesh.m_vds[i]->resize(j);
  }

  for (size_t i = 0; i < vm.size(); i++) {
    if (vm[i] >= 0) {
      pm_positions[vm[i]] = pm_positions[i];
      pm_newpositions[vm[i]] = pm_newpositions[i];
      m_masses[vm[i]] = m_masses[i];
      m_mesh.m_is_boundary_vertex[vm[i]] = m_mesh.m_is_boundary_vertex[i];
      m_mesh.m_vertex_to_edge_map[vm[i]] = m_mesh.m_vertex_to_edge_map[i];
      m_mesh.m_vertex_to_triangle_map[vm[i]] =
          m_mesh.m_vertex_to_triangle_map[i];
    }
  }
  pm_positions.resize(j);
  pm_newpositions.resize(j);
  m_masses.resize(j);
  m_mesh.m_is_boundary_vertex.resize(j);
  m_mesh.m_vertex_to_edge_map.resize(j);
  m_mesh.m_vertex_to_triangle_map.resize(j);

  for (size_t i = 0; i < m_mesh.m_edges.size(); i++) {
    if (m_mesh.edge_is_deleted(i)) continue;

    assert(vm[m_mesh.m_edges[i][0]] >=
           0);  // an existing edge cannot reference a deleted vertex
    m_mesh.m_edges[i][0] = vm[m_mesh.m_edges[i][0]];
    assert(vm[m_mesh.m_edges[i][1]] >=
           0);  // an existing edge cannot reference a deleted vertex
    m_mesh.m_edges[i][1] = vm[m_mesh.m_edges[i][1]];
  }

  for (size_t i = 0; i < m_mesh.m_tris.size(); i++) {
    if (m_mesh.triangle_is_deleted(i)) continue;

    for (size_t j = 0; j < 3; j++) {
      assert(vm[m_mesh.m_tris[i][j]] >=
             0);  // an existing triangle cannot reference a deleted vertex
      m_mesh.m_tris[i][j] = vm[m_mesh.m_tris[i][j]];
    }
  }

  // defragment edges
  std::vector<int> em(m_mesh.ne(), -1);
  j = 0;
  for (size_t i = 0; i < m_mesh.ne(); i++)
    if (m_mesh.m_edge_to_triangle_map[i].size() != 0) em[i] = j++;

  for (size_t i = 0; i < m_mesh.m_eds.size(); i++) {
    m_mesh.m_eds[i]->compress(em);
    m_mesh.m_eds[i]->resize(j);
  }

  for (size_t i = 0; i < em.size(); i++) {
    if (em[i] >= 0) {
      m_mesh.m_is_boundary_edge[em[i]] = m_mesh.m_is_boundary_edge[i];
      m_mesh.m_edges[em[i]] = m_mesh.m_edges[i];
      m_mesh.m_edge_to_triangle_map[em[i]] = m_mesh.m_edge_to_triangle_map[i];
    }
  }
  m_mesh.m_is_boundary_edge.resize(j);
  m_mesh.m_edges.resize(j);
  m_mesh.m_edge_to_triangle_map.resize(j);

  for (size_t i = 0; i < m_mesh.m_vertex_to_edge_map.size(); i++) {
    if (m_mesh.vertex_is_deleted(i)) continue;

    for (size_t j = 0; j < m_mesh.m_vertex_to_edge_map[i].size(); j++) {
      assert(em[m_mesh.m_vertex_to_edge_map[i][j]] >=
             0);  // an existing vertex cannot be referenced by a deleted edge
      m_mesh.m_vertex_to_edge_map[i][j] = em[m_mesh.m_vertex_to_edge_map[i][j]];
    }
  }

  for (size_t i = 0; i < m_mesh.m_triangle_to_edge_map.size(); i++) {
    if (m_mesh.triangle_is_deleted(i)) continue;

    for (size_t j = 0; j < 3; j++) {
      assert(em[m_mesh.m_triangle_to_edge_map[i][j]] >=
             0);  // an existing triangle cannot reference a delete edge
      m_mesh.m_triangle_to_edge_map[i][j] =
          em[m_mesh.m_triangle_to_edge_map[i][j]];
    }
  }

  // defragment triangles
  std::vector<int> fm(m_mesh.nt(), -1);
  j = 0;
  for (size_t i = 0; i < m_mesh.nt(); i++)
    if (!m_mesh.triangle_is_deleted(i)) fm[i] = j++;

  for (size_t i = 0; i < m_mesh.m_fds.size(); i++) {
    m_mesh.m_fds[i]->compress(fm);
    m_mesh.m_fds[i]->resize(j);
  }

  for (size_t i = 0; i < fm.size(); i++) {
    if (fm[i] >= 0) {
      m_mesh.m_tris[fm[i]] = m_mesh.m_tris[i];
      m_mesh.m_triangle_labels[fm[i]] = m_mesh.m_triangle_labels[i];
      m_mesh.m_triangle_to_edge_map[fm[i]] = m_mesh.m_triangle_to_edge_map[i];
    }
  }
  m_mesh.m_tris.resize(j);
  m_mesh.m_triangle_labels.resize(j);
  m_mesh.m_triangle_to_edge_map.resize(j);

  for (size_t i = 0; i < m_mesh.m_vertex_to_triangle_map.size(); i++) {
    if (m_mesh.vertex_is_deleted(i)) continue;

    for (size_t j = 0; j < m_mesh.m_vertex_to_triangle_map[i].size(); j++) {
      assert(
          fm[m_mesh.m_vertex_to_triangle_map[i][j]] >=
          0);  // an existing vertex cannot be referenced by a deleted triangle
      m_mesh.m_vertex_to_triangle_map[i][j] =
          fm[m_mesh.m_vertex_to_triangle_map[i][j]];
    }
  }

  for (size_t i = 0; i < m_mesh.m_edge_to_triangle_map.size(); i++) {
    if (m_mesh.edge_is_deleted(i)) continue;

    for (size_t j = 0; j < m_mesh.m_edge_to_triangle_map[i].size(); j++) {
      assert(fm[m_mesh.m_edge_to_triangle_map[i][j]] >=
             0);  // an existing edge cannot be referenced by a deleted triangle
      m_mesh.m_edge_to_triangle_map[i][j] =
          fm[m_mesh.m_edge_to_triangle_map[i][j]];
    }
  }

  double end_time = get_time_in_seconds();
  g_stats.add_to_double("total_defrag_time", end_time - start_time);

  m_mesh.test_connectivity();

  if (m_collision_safety) {
    rebuild_continuous_broad_phase();
  }
}

void SurfTrack::defrag_mesh_from_scratch_copy(
    std::vector<size_t>& vertices_to_be_mapped) {
  double start_time = get_time_in_seconds();

  // defragment vertices
  std::vector<int> vm(get_num_vertices(), -1);
  size_t j = 0;
  for (size_t i = 0; i < get_num_vertices(); i++)
    if (!m_mesh.vertex_is_deleted(i)) vm[i] = j++;

  for (size_t i = 0; i < vertices_to_be_mapped.size(); i++)
    vertices_to_be_mapped[i] = vm[vertices_to_be_mapped[i]];

  for (size_t i = 0; i < m_mesh.m_vds.size(); i++) {
    m_mesh.m_vds[i]->compress(vm);
    m_mesh.m_vds[i]->resize(j);
  }

  for (size_t i = 0; i < vm.size(); i++) {
    if (vm[i] >= 0) {
      pm_positions[vm[i]] = pm_positions[i];
      pm_newpositions[vm[i]] = pm_newpositions[i];
      m_masses[vm[i]] = m_masses[i];
      m_mesh.m_is_boundary_vertex[vm[i]] = m_mesh.m_is_boundary_vertex[i];
      m_mesh.m_vertex_to_edge_map[vm[i]] = m_mesh.m_vertex_to_edge_map[i];
      m_mesh.m_vertex_to_triangle_map[vm[i]] =
          m_mesh.m_vertex_to_triangle_map[i];
    }
  }
  pm_positions.resize(j);
  pm_newpositions.resize(j);
  m_masses.resize(j);
  m_mesh.m_is_boundary_vertex.resize(j);
  m_mesh.m_vertex_to_edge_map.resize(j);
  m_mesh.m_vertex_to_triangle_map.resize(j);

  for (size_t i = 0; i < m_mesh.m_edges.size(); i++) {
    if (m_mesh.edge_is_deleted(i)) continue;

    assert(vm[m_mesh.m_edges[i][0]] >=
           0);  // an existing edge cannot reference a deleted vertex
    m_mesh.m_edges[i][0] = vm[m_mesh.m_edges[i][0]];
    assert(vm[m_mesh.m_edges[i][1]] >=
           0);  // an existing edge cannot reference a deleted vertex
    m_mesh.m_edges[i][1] = vm[m_mesh.m_edges[i][1]];
  }

  for (size_t i = 0; i < m_mesh.m_tris.size(); i++) {
    if (m_mesh.triangle_is_deleted(i)) continue;

    for (size_t j = 0; j < 3; j++) {
      assert(vm[m_mesh.m_tris[i][j]] >=
             0);  // an existing triangle cannot reference a deleted vertex
      m_mesh.m_tris[i][j] = vm[m_mesh.m_tris[i][j]];
    }
  }

  // defragment triangles and edges together by rebuilding the connectivity
  // the code below follows NonDestructiveTriMesh::clear_deleted_triangles();
  // the difference is that the data in m_eds and m_fds are salvaged and copied
  // to the right place.
  std::vector<Vec3st> new_tris;
  std::vector<Vec2i> new_labels;
  new_tris.reserve(m_mesh.nt());
  new_labels.reserve(m_mesh.nt());

  std::vector<int> fm(m_mesh.nt(), -1);
  j = 0;
  for (size_t i = 0; i < m_mesh.nt(); i++) {
    if (!m_mesh.triangle_is_deleted(i)) {
      fm[i] = j++;
      new_tris.push_back(m_mesh.m_tris[i]);
      new_labels.push_back(m_mesh.m_triangle_labels[i]);
    }
  }
  size_t new_nf = j;

  m_mesh.m_tris = new_tris;
  m_mesh.m_triangle_labels = new_labels;

  std::vector<Vec2st> old_edges = m_mesh.m_edges;

  m_mesh.update_connectivity();  // rebuild triangles and edges

  // infer the edge defrag map
  std::vector<int> em(old_edges.size(), -1);
  for (size_t i = 0; i < old_edges.size(); i++) {
    if (old_edges[i][0] == old_edges[i][1])  // deleted
      continue;

    size_t e = m_mesh.get_edge_index(vm[old_edges[i][0]], vm[old_edges[i][1]]);
    assert(e != m_mesh.ne());  // should be found
    em[i] = e;
  }
  size_t new_ne = m_mesh.ne();

  for (size_t i = 0; i < m_mesh.m_eds.size(); i++) {
    m_mesh.m_eds[i]->compress(em);
    m_mesh.m_eds[i]->resize(new_ne);
  }

  for (size_t i = 0; i < m_mesh.m_fds.size(); i++) {
    m_mesh.m_fds[i]->compress(fm);
    m_mesh.m_fds[i]->resize(new_nf);
  }

  double end_time = get_time_in_seconds();
  g_stats.add_to_double("total_defrag_time", end_time - start_time);

  m_mesh.test_connectivity();

  if (m_collision_safety) {
    rebuild_continuous_broad_phase();
  }
}

// --------------------------------------------------------
///
/// Fire an assert if any triangle has repeated vertices or if any zero-volume
/// tets are found.
///
// --------------------------------------------------------

void SurfTrack::assert_no_degenerate_triangles() {
  // for each triangle on the surface
  for (size_t i = 0; i < m_mesh.num_triangles(); ++i) {
    const Vec3st& current_triangle = m_mesh.get_triangle(i);

    if ((current_triangle[0] == 0) && (current_triangle[1] == 0) &&
        (current_triangle[2] == 0)) {
      // deleted triangle
      continue;
    }

    //
    // check if triangle has repeated vertices
    //

    assert(!((current_triangle[0] == current_triangle[1]) ||
             (current_triangle[1] == current_triangle[2]) ||
             (current_triangle[2] == current_triangle[0])));

    //
    // look for flaps
    //
    const Vec3st& tri_edges = m_mesh.m_triangle_to_edge_map[i];

    bool flap_found = false;

    for (unsigned int e = 0; e < 3 && flap_found == false; ++e) {
      const std::vector<size_t>& edge_tris =
          m_mesh.m_edge_to_triangle_map[tri_edges[e]];

      for (size_t t = 0; t < edge_tris.size(); ++t) {
        if (edge_tris[t] == i) {
          continue;
        }

        size_t other_triangle_index = edge_tris[t];
        const Vec3st& other_triangle =
            m_mesh.get_triangle(other_triangle_index);

        if ((other_triangle[0] == other_triangle[1]) ||
            (other_triangle[1] == other_triangle[2]) ||
            (other_triangle[2] == other_triangle[0])) {
          assert(!"repeated vertices");
        }

        if (((current_triangle[0] == other_triangle[0]) ||
             (current_triangle[0] == other_triangle[1]) ||
             (current_triangle[0] == other_triangle[2])) &&
            ((current_triangle[1] == other_triangle[0]) ||
             (current_triangle[1] == other_triangle[1]) ||
             (current_triangle[1] == other_triangle[2])) &&
            ((current_triangle[2] == other_triangle[0]) ||
             (current_triangle[2] == other_triangle[1]) ||
             (current_triangle[2] == other_triangle[2]))) {
          size_t common_edge = tri_edges[e];
          if (m_mesh.oriented(m_mesh.m_edges[common_edge][0],
                              m_mesh.m_edges[common_edge][1],
                              current_triangle) ==
              m_mesh.oriented(m_mesh.m_edges[common_edge][0],
                              m_mesh.m_edges[common_edge][1], other_triangle)) {
            assert(false);
            continue;
          }

          assert(false);
        }
      }
    }
  }
}

// --------------------------------------------------------
///
/// Fire an assert if any triangle has repeated vertices or if any zero-volume
/// tets are found.
///
// --------------------------------------------------------

bool SurfTrack::any_triangles_with_bad_angles() {
  // for each triangle on the surface
  for (size_t i = 0; i < m_mesh.num_triangles(); ++i)
    if (triangle_with_bad_angle(i)) return true;

  return false;
}

bool SurfTrack::triangle_with_bad_angle(size_t i) {
  Vec3st tri = m_mesh.m_tris[i];
  if (m_mesh.triangle_is_deleted(i)) return false;

  assert(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);

  Vec3d v0 = get_position(tri[0]);
  Vec3d v1 = get_position(tri[1]);
  Vec3d v2 = get_position(tri[2]);

  double min_angle = min_triangle_angle(v0, v1, v2);
  double max_angle = max_triangle_angle(v0, v1, v2);

  // these simply must be true at all times
  assert(min_angle >= 0);
  assert(max_angle < 1.000001 * M_PI);
  assert(min_angle == min_angle);
  assert(max_angle == max_angle);

  // if any triangles are outside our bounds, we have to keep going.

  if (rad2deg(min_angle) < m_min_triangle_angle ||
      rad2deg(max_angle) >= m_max_triangle_angle)
    return true;
  return false;
}

// --------------------------------------------------------
///
/// Delete flaps and zero-area triangles.
///
// --------------------------------------------------------

void SurfTrack::trim_degeneracies(std::vector<size_t>& triangle_indices) {
  // If we're not allowing non-manifold, assert we don't have any

  if (false == m_allow_non_manifold) {
    // check for edges incident on more than 2 triangles
    for (size_t i = 0; i < m_mesh.m_edge_to_triangle_map.size(); ++i) {
      if (m_mesh.edge_is_deleted(i)) {
        continue;
      }
      assert(m_mesh.m_edge_to_triangle_map[i].size() == 1 ||
             m_mesh.m_edge_to_triangle_map[i].size() == 2);
    }

    triangle_indices.clear();
    return;
  }

  for (size_t j = 0; j < triangle_indices.size(); ++j) {
    size_t i = triangle_indices[j];

    const Vec3st& current_triangle = m_mesh.get_triangle(i);

    if ((current_triangle[0] == 0) && (current_triangle[1] == 0) &&
        (current_triangle[2] == 0)) {
      continue;
    }

    //
    // look for triangles with repeated vertices
    //
    if ((current_triangle[0] == current_triangle[1]) ||
        (current_triangle[1] == current_triangle[2]) ||
        (current_triangle[2] == current_triangle[0])) {
      if (m_verbose) {
        std::cout << "deleting degenerate triangle " << i << ": "
                  << current_triangle << std::endl;
      }

      // delete it
      remove_triangle(i);

      continue;
    }

    //
    // look for flaps
    //
    const Vec3st& tri_edges = m_mesh.m_triangle_to_edge_map[i];

    bool flap_found = false;

    for (unsigned int e = 0; e < 3 && flap_found == false; ++e) {
      const std::vector<size_t>& edge_tris =
          m_mesh.m_edge_to_triangle_map[tri_edges[e]];

      for (size_t t = 0; t < edge_tris.size(); ++t) {
        if (edge_tris[t] == i) {
          continue;
        }

        size_t other_triangle_index = edge_tris[t];
        const Vec3st& other_triangle =
            m_mesh.get_triangle(other_triangle_index);

        //(just tests if the tri has repeated vertices)
        if (m_mesh.triangle_is_deleted(other_triangle_index)) continue;

        if (((current_triangle[0] == other_triangle[0]) ||
             (current_triangle[0] == other_triangle[1]) ||
             (current_triangle[0] == other_triangle[2])) &&
            ((current_triangle[1] == other_triangle[0]) ||
             (current_triangle[1] == other_triangle[1]) ||
             (current_triangle[1] == other_triangle[2])) &&
            ((current_triangle[2] == other_triangle[0]) ||
             (current_triangle[2] == other_triangle[1]) ||
             (current_triangle[2] == other_triangle[2]))) {
          if (false == m_allow_topology_changes) {
            std::cout << "flap found while topology changes disallowed"
                      << std::endl;
            std::cout << current_triangle << std::endl;
            std::cout << other_triangle << std::endl;
            assert(0);
          }

          Vec2i current_label = m_mesh.get_triangle_label(i);
          Vec2i other_label = m_mesh.get_triangle_label(other_triangle_index);

          size_t common_edge = tri_edges[e];
          bool orientation =
              (m_mesh.oriented(m_mesh.m_edges[common_edge][0],
                               m_mesh.m_edges[common_edge][1],
                               current_triangle) ==
               m_mesh.oriented(m_mesh.m_edges[common_edge][0],
                               m_mesh.m_edges[common_edge][1], other_triangle));

          int region_0;  // region behind surface a
          int region_1;  // region behind surface b
          int region_2;  // region between the two surfaces

          if (current_label[0] == other_label[0]) {
            region_0 = current_label[1];
            region_1 = other_label[1];
            region_2 = current_label[0];
            assert(!orientation);  // two triangles should have opposite
                                   // orientation
          } else if (current_label[1] == other_label[0]) {
            region_0 = current_label[0];
            region_1 = other_label[1];
            region_2 = current_label[1];
            assert(orientation);  // two triangles should have same orientation
          } else if (current_label[0] == other_label[1]) {
            region_0 = current_label[1];
            region_1 = other_label[0];
            region_2 = current_label[0];
            assert(orientation);  // two triangles should have same orientation
          } else if (current_label[1] == other_label[1]) {
            region_0 = current_label[0];
            region_1 = other_label[0];
            region_2 = current_label[1];
            assert(!orientation);  // two triangles should have opposite
                                   // orientation
          } else {
            // shouldn't happen
            assert(!"Face label inconsistency detected.");
          }

          std::vector<size_t> to_delete;
          std::vector<std::pair<size_t, Vec2i> > dirty;
          if (region_0 == region_1) {
            // Two parts of the same region meeting upon removal of the region
            // in between. This entire flap pair should be removed to open a
            // window connecting both sides.

            to_delete.push_back(i);
            to_delete.push_back(other_triangle_index);

          } else {
            // Three different regions. After the region in between is removed,
            // there should be an interface remaining to separate the two
            // regions on two sides

            to_delete.push_back(other_triangle_index);

            if (current_label[0] == region_0)
              m_mesh.set_triangle_label(i, Vec2i(region_0, region_1));
            else
              m_mesh.set_triangle_label(i, Vec2i(region_1, region_0));

            dirty.push_back(
                std::pair<size_t, Vec2i>(i, m_mesh.get_triangle_label(i)));
          }

          // the dangling vertex will be safely removed by the vertex cleanup
          // function

          // delete the triangles

          if (m_verbose) {
            std::cout << "flap: triangles << " << i << " [" << current_triangle
                      << "] and " << edge_tris[t] << " [" << other_triangle
                      << "]" << std::endl;
          }

          for (size_t k = 0; k < to_delete.size(); k++)
            remove_triangle(to_delete[k]);

          MeshUpdateEvent flap_delete(MeshUpdateEvent::FLAP_DELETE);
          flap_delete.m_deleted_tris = to_delete;
          flap_delete.m_dirty_tris = dirty;
          m_mesh_change_history.push_back(flap_delete);

          ////////////////////////////////////////////////////////////

          flap_found = true;
          break;
        }
      }
    }
  }

  triangle_indices.clear();
}

// --------------------------------------------------------
///
/// One pass: split long edges, flip non-delaunay edges, collapse short edges,
/// null-space smoothing
///
// --------------------------------------------------------

void SurfTrack::improve_mesh() {
  if (m_mesheventcallback)
    m_mesheventcallback->log() << "Improve mesh began" << std::endl;

  if (m_perform_improvement) {
    ////////////////////////////////////////////////////////////

    // standard mesh improvement pass is gentle, seeks to preserve volume,
    // prevent flips, preserve features, avoid popping.
    m_aggressive_mode = false;

    int i = 0;

    // edge splitting
    std::cout << "Splits\n";
    while (m_splitter.split_pass()) {
      if (m_mesheventcallback)
        m_mesheventcallback->log()
            << "Split pass " << i << " finished" << std::endl;
      i++;
      std::cout << "Splits\n";
    }

    //        // edge flipping
    //        std::cout << "Flips\n";
    //        m_flipper.flip_pass();
    //        if (m_mesheventcallback)
    //            m_mesheventcallback->log() << "Flip pass finished" <<
    //            std::endl;

    // edge collapsing
    i = 0;
    std::cout << "Collapses\n";
    while (m_collapser.collapse_pass()) {
      if (m_mesheventcallback)
        m_mesheventcallback->log()
            << "Collapse pass " << i << " finished" << std::endl;
      i++;
      std::cout << "Collapses\n";
    }

    // process t1 transitions (vertex separation)
    i = 0;

    while (m_t1_transition_enabled && m_t1transition.t1_pass()) {
      std::cout << "T1's\n";
      if (m_mesheventcallback)
        m_mesheventcallback->log()
            << "T1 pass " << i << " finished" << std::endl;
      i++;
    }

    //        // smoothing
    //        if ( m_perform_smoothing)
    //        {
    //            std::cout << "Smoothing\n";
    //            m_smoother.null_space_smoothing_pass( 1.0 );
    //            if (m_mesheventcallback)
    //                m_mesheventcallback->log() << "Smoothing pass finished" <<
    //                std::endl;
    //        }

    //
    //        ////////////////////////////////////////////////////////////
    //        //enter aggressive improvement mode to improve remaining bad
    //        triangles up to minimum bounds,
    //        //potentially at the expense of some shape deterioration. (aka.
    //        BEAST MODE!!1!1!!)
    //
    //        m_aggressive_mode = true;
    //        //m_verbose = true;
    //        i = 0;
    //        while(any_triangles_with_bad_angles() && i < 100) {
    //            //enter aggressive mode
    //
    //            std::cout << "Aggressive mesh improvement iteration #" << i <<
    //            "." << std::endl;
    //
    //            m_splitter.split_pass();
    //            if (m_mesheventcallback)
    //                m_mesheventcallback->log() << "Aggressive split pass " <<
    //                i << " finished" << std::endl;
    //
    //            //switch to delaunay criterion for this, since it is purported
    //            to produce better angles for a given vertex set.
    //            m_flipper.m_use_Delaunay_criterion = true;
    //            m_flipper.flip_pass();
    //            if (m_mesheventcallback)
    //                m_mesheventcallback->log() << "Aggressive flip pass " << i
    //                << " finished" << std::endl;
    //            m_flipper.m_use_Delaunay_criterion = false; //switch back to
    //            valence-based mode
    //
    //            //try to cut out early if things have already gotten better.
    //            if(!any_triangles_with_bad_angles())
    //                break;
    //
    //            m_collapser.collapse_pass();
    //            if (m_mesheventcallback)
    //                m_mesheventcallback->log() << "Aggressive collapse pass "
    //                << i << " finished" << std::endl;
    //
    //            //try to cut out early if things have already gotten better.
    //            if(!any_triangles_with_bad_angles())
    //                break;
    //
    //            if (m_perform_smoothing)
    //            {
    //                m_smoother.null_space_smoothing_pass( 1.0 );
    //                if (m_mesheventcallback)
    //                    m_mesheventcallback->log() << "Aggressive smoothing
    //                    pass " << i << " finished" << std::endl;
    //            }
    //
    //            i++;
    //            //start dumping warning messages if we're doing a lot of
    //            iterations.
    //
    //        }
    //
    //        m_aggressive_mode = false;
    //        //m_verbose = false;
    //
    //        assert_no_bad_labels();

    std::cout << "Done improvement\n" << std::endl;
    if (m_collision_safety) {
      assert_mesh_is_intersection_free(false);
    }
  }

  if (m_mesheventcallback)
    m_mesheventcallback->log() << "Improve mesh finished" << std::endl;
}

void SurfTrack::cut_mesh(const std::vector<std::pair<size_t, size_t> >& edges) {
  // edge cutting
  m_cutter.separate_edges_new(edges);

  if (m_collision_safety) {
    // std::cout << "Checking collisions after cutting.\n";
    assert_mesh_is_intersection_free(false);
  }
}

// --------------------------------------------------------
///
/// Perform a pass of merge attempts
///
// --------------------------------------------------------

void SurfTrack::topology_changes() {
  if (m_mesheventcallback)
    m_mesheventcallback->log() << "Topology changes began" << std::endl;

  if (false == m_allow_topology_changes) {
    return;
  }

  // bool merge_occurred = merge_occurred = m_merger.merge_pass(); //OLD MERGING
  // CODE
  bool merge_occurred = m_snapper.snap_pass();  // NEW MERGING CODE

  if (m_mesheventcallback)
    m_mesheventcallback->log() << "Snap pass finished" << std::endl;

  if (m_collision_safety) {
    assert_mesh_is_intersection_free(false);
  }

  if (m_mesheventcallback)
    m_mesheventcallback->log() << "Topology changes finished" << std::endl;
}

void SurfTrack::assert_no_bad_labels() {
  for (size_t i = 0; i < m_mesh.m_triangle_labels.size(); ++i) {
    if (m_mesh.triangle_is_deleted(i))  // skip dead tris.
      continue;

    Vec2i label = m_mesh.get_triangle_label(i);
    assert(label[0] != -1 && "Uninitialized label");
    assert(label[1] != -1 && "Uninitialized label");

    // TODO The next might have to change for open regions.
    assert(label[0] != label[1] && "Same front and back labels");
  }
}

int SurfTrack::vertex_feature_edge_count(size_t vertex) const {
  int count = 0;
  for (size_t i = 0; i < m_mesh.m_vertex_to_edge_map[vertex].size(); ++i) {
    count += (edge_is_feature(m_mesh.m_vertex_to_edge_map[vertex][i]) ? 1 : 0);
  }
  return count;
}
int SurfTrack::vertex_feature_edge_count(
    size_t vertex, const std::vector<Vec3d>& cached_normals) const {
  int count = 0;
  for (size_t i = 0; i < m_mesh.m_vertex_to_edge_map[vertex].size(); ++i) {
    count +=
        (edge_is_feature(m_mesh.m_vertex_to_edge_map[vertex][i], cached_normals)
             ? 1
             : 0);
  }
  return count;
}

// Return whether the given edge is a feature.
bool SurfTrack::edge_is_feature(size_t edge) const {
  if (edge_is_any_solid(edge)) {
    assert(m_solid_vertices_callback);
    if (m_solid_vertices_callback->solid_edge_is_feature(*this, edge))
      return true;
  }
  return get_largest_dihedral(edge) > m_feature_edge_angle_threshold;
}

bool SurfTrack::edge_is_feature(
    size_t edge, const std::vector<Vec3d>& cached_normals) const {
  if (edge_is_any_solid(edge)) {
    assert(m_solid_vertices_callback);
    if (m_solid_vertices_callback->solid_edge_is_feature(*this, edge))
      return true;
  }
  return get_largest_dihedral(edge, cached_normals) >
         m_feature_edge_angle_threshold;
}

bool SurfTrack::vertex_feature_is_smooth_ridge(size_t vertex) const {
  std::vector<size_t> feature_edges;
  for (size_t i = 0; i < m_mesh.m_vertex_to_edge_map[vertex].size(); i++)
    if (edge_is_feature(m_mesh.m_vertex_to_edge_map[vertex][i]))
      feature_edges.push_back(m_mesh.m_vertex_to_edge_map[vertex][i]);

  if (feature_edges.size() != 2)  // the feature cannot be a smooth ridge unless
                                  // the number of feature edges is 2
    return false;

  const Vec2st& e0 = m_mesh.m_edges[feature_edges[0]];
  const Vec2st& e1 = m_mesh.m_edges[feature_edges[1]];

  Vec3d t0 = get_position(e0[0]) - get_position(e0[1]);
  t0 /= mag(t0);
  Vec3d t1 = get_position(e1[0]) - get_position(e1[1]);
  t1 /= mag(t1);

  // adjust for orientation
  if ((e0[0] == vertex && e1[0] == vertex) ||
      (e0[1] == vertex && e1[1] == vertex))
    t1 = -t1;

  double angle = acos(dot(t0, t1));
  if (angle > m_feature_edge_angle_threshold) return false;

  return true;
}

}  // namespace LosTopos
