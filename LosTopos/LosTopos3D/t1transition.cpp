// ---------------------------------------------------------
//
//  t1transition.cpp
//  Christopher Batty, Fang Da 2014
//
//  Functions handling T1 transitions (edge popping and vertex popping).
//
// ---------------------------------------------------------

#include <queue>
#include <set>

#include <t1transition.h>
#include <broadphase.h>
#include <collisionqueries.h>
#include <runstats.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>


// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace LosTopos {
    
extern RunStats g_stats;

struct T1Transition::InteriorStencil
{
    Vec3st vertices;
    Vec3st vertex_indices;
    Vec2st edges;
    Vec2st edge_indices;
};



// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor.  Active SurfTrack object must be supplied.
///
// ---------------------------------------------------------

T1Transition::T1Transition(SurfTrack & surf, VelocityFieldCallback * vfc, bool remesh_boundaries) :
m_remesh_boundaries(remesh_boundaries),
m_pull_apart_distance(0.002),
m_pull_apart_tendency_threshold(0),
m_surf(surf),
m_velocity_field_callback(vfc)
{
    
}

template <class S, class T>
struct less_pair_first
{
    bool operator() (const std::pair<S, T> & x, const std::pair<S, T> & y) const { return x.first < y.first; }
};

template <class S, class T>
struct less_pair_second
{
    bool operator() (const std::pair<S, T> & x, const std::pair<S, T> & y) const { return x.second < y.second; }
};

struct SortableDirectionCandidate
{
    size_t vertex;
    Vec2i regions;
    Vec3d direction;
    double tendency;
    
    bool operator < (const SortableDirectionCandidate & sdc) const { return tendency < sdc.tendency; }
};

// --------------------------------------------------------
///
/// Perform vertex popping
///
// --------------------------------------------------------
bool T1Transition::t1_pass()
{
    if (m_surf.m_verbose)
        std::cout << "---------------------- T1 Transition: Vertex popping ----------------------" << std::endl;
    
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    bool pop_occurred = false;
    
    // find the region count
    int max_region = -1;
    for (size_t i = 0; i < mesh.nt(); i++)
    {
        if (mesh.get_triangle(i)[0] == mesh.get_triangle(i)[1] && mesh.get_triangle(i)[0] == mesh.get_triangle(i)[2])
            continue;
        Vec2i label = mesh.get_triangle_label(i);
        assert(label[0] >= 0);
        assert(label[1] >= 0);
        if (label[0] > max_region) max_region = label[0];
        if (label[1] > max_region) max_region = label[1];
    }
    int nregion = max_region + 1;
    
    // this is the incident matrix for a region graph. not every col/row need to be filled for a particular vertex becuase the vertex may not be incident to all regions.
    std::vector<std::vector<int> > region_graph(nregion, std::vector<int>(nregion, 0));
    
    // a list of candidate directions
    std::vector<SortableDirectionCandidate> candidates;
    
    // loop through all the vertices
    for (size_t i = 0; i < mesh.nv(); i++)
    {
        size_t xj = i;
        
        std::set<int> vertex_regions_set;
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[xj][i]);
            vertex_regions_set.insert(label[0]);
            vertex_regions_set.insert(label[1]);
        }
        
        std::vector<int> vertex_regions;
        vertex_regions.assign(vertex_regions_set.begin(), vertex_regions_set.end());
        
        // cull away the interior vertices of manifold surface patches
        if (vertex_regions_set.size() < 3)
            continue;
        
        //
        // Pull apart strategy: in order to cope with various configurations (such as a region-valence-5 vertex being a triple junction, or
        //  junctions on the BB walls), the following strategy is adopted:
        //
        // The regions incident on the vertex form a graph, with each region being a node and an edge exists between two regions iff the two
        //  regions share a triangle face in the mesh. The objective at the end of this processing, is to make sure this graph is complete for
        //  every vertex. For example the original region-valence-4 X-junction vertex (the "hourglass" junction) forms a graph like this:
        //
        //  A o--o B
        //    | /|
        //    |/ |
        //  C o--o D
        //
        //  where the edge between node A and D (the two bulbs of the hourglass) is missing because the two regions only meet at one vertex
        //  (the hourglass neck).
        //
        // The processing here makes an incomplete graph complete by pulling apart vertices. If two nodes (e.g. A and D in the figure above)
        //  do not share an edge, it means the center vertex is the only interface between the two regions, and it needs to be pulled apart.
        //  Region A and D each keep one of the two duplicates of the vertex. If a region has edges with both regions A and D, such as region
        //  C and B, will contain both duplicates (the region's shape is turned from a cone into a flat screwdriver's tip). Now look at the
        //  resulting graphs for the two duplicate vertices. For the duplicate that goes with region A, region D is no longer incident and
        //  thus removed from the graph:
        //
        //  A o--o B
        //    | /
        //    |/
        //  C o
        //
        //  and this is now a complete graph. Similarly the graph for the duplicate that goes with region D will not have node A, and it is
        //  complete too. Of course in more complex scenarios the graphs may not directly become complete immediately after we process one
        //  missing edge. We will visit the two resulting vertices again as if they are a regular vertex that may need to be pulled apart too.
        //
        // The algorithm pseudocode is as following:
        //
        //  Pop a vertex from the stack of vertices to be processed, construct a graph of region adjacency
        //  If the graph is already complete, skip this vertex
        //  Find the pair of disconnected nodes A and B from the graph that has the strongest tensile force after pulling apart (there may be more than one pair)
        //  Pull apart the vertex into two (a and b, corresponding to region A and B respectively), initialized to have the same coordinates and then pulled apart
        //  For each face incident to the center vertex in the original mesh
        //    If it is incident to region A, update it to use vertex a
        //    If it is incident to region B, update it to use vertex b
        //    Otherwise, update it to use vertex b
        //  For each edge incident to the center vertex in the original mesh, and incident to any face with A label
        //    Create a new triangle with both a and b
        //  Push vertex a and b on the stack to be visited next.
        //
        // The first loop essentially pulls region A from everything else around the vertex, leaving a ring of blank around the vertex. Then
        //  the next loop fills the blank. Note that this algorithm works even in presence of unresolved X-junction edges, if any (due to
        //  collision etc). Another interpretation is: region A and region B remain intact; all the other triangles (incident to neither A
        //  nor B) form a number of triangle fans, which may share edges/triangles between them. There are a number of "head" edges, which
        //  are located at the head of the fans from region A to region B. These edges are the edges incident to region A. The second loop
        //  basically sweeps these head edges into a triangles (from a to b).
        //
        // Collision test is performed before the operation to ensure collision safety. First the center vertex is moved to the position of
        //  new vertex b without changing connectivity. Collision is checked for this motion. Then vertex a is separated from vertex b,
        //  bringing all region A faces along with it. Collision is checked for this motion, using only region A faces as incident faces and
        //  region A edges as incident edges. This CCD only detects the collision of this motion against the rest of the mesh, excluding all
        //  faces/edges connected to both a and b (b faces/edges are excluded too because initially a and b coincide). Then an instantaneous
        //  intersection test is perform that finds intersections in the final configuration.
        //
        //
        
        // construct the region graph
        //region_graph.assign(nregion, std::vector<int>(nregion, 0));
        region_graph.resize(nregion); // Fang's big speedup
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[xj][i]);
            region_graph[label[0]][label[1]] = 1;
            region_graph[label[1]][label[0]] = 1;
        }
        
        // find missing edges, as candidates of pull-apart
        std::vector<std::pair<double, std::pair<Vec2i, Vec3d> > > candidate_pairs;  // components: <tensile_force, <(A, B), pull_apart_direction> >
        for (size_t i = 0; i < vertex_regions.size(); i++)
        {
            for (size_t j = i + 1; j < vertex_regions.size(); j++)
            {
                int A = vertex_regions[i];
                int B = vertex_regions[j];
                if (region_graph[A][B] == 0)
                {
                    Vec3d pull_apart_direction;
                    double pull_apart_tendency = 0;
                    if (m_velocity_field_callback)
                        pull_apart_tendency = try_pull_vertex_apart_using_velocity_field(xj, A, B, pull_apart_direction);
                    else
                        pull_apart_tendency = try_pull_vertex_apart_using_surface_tension(xj, A, B, pull_apart_direction);
                    
                    if (pull_apart_tendency > 0)
                    {
                        SortableDirectionCandidate candidate;
                        candidate.vertex = xj;
                        candidate.regions = Vec2i(A, B);
                        candidate.direction = pull_apart_direction;
                        candidate.tendency = pull_apart_tendency;
                        candidates.push_back(candidate);
                    }
                }
            }
        }
    }
    
    // sort the candidate pairs according to the strength of the tensile force
    std::sort(candidates.begin(), candidates.end());
    
    // process the candidates from top
    for ( ; candidates.size() > 0; candidates.pop_back())
    {
        size_t xj = candidates.back().vertex;
        int A = candidates.back().regions[0];
        int B = candidates.back().regions[1];
        Vec3d pull_apart_direction = candidates.back().direction;
        double pull_apart_tendency = candidates.back().tendency;
        
        // check if this candidate is still valid
        // due to processing other candidates prior to this, the vertex xj may have been deleted, region A and B may have
        //  already come into contact, or the tendency may have dropped to negative.
        if (mesh.m_vertex_to_triangle_map[xj].size() == 0)
            continue;
        
        bool contact = false;
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[xj][i]);
            if ((label[0] == A && label[1] == B) || (label[0] == B && label[1] == A))
            {
                contact = true;
                break;
            }
        }
        if (contact)
            continue;
        
        if (m_velocity_field_callback)
            pull_apart_tendency = try_pull_vertex_apart_using_velocity_field(xj, A, B, pull_apart_direction);
        else
            pull_apart_tendency = try_pull_vertex_apart_using_surface_tension(xj, A, B, pull_apart_direction);
        if (pull_apart_tendency < 0)
            continue;
        
        if (m_surf.m_verbose)
        {
            std::cout << "Attempting to pop vertex " << xj << " region " << A << " from region " << B << std::endl;
        }
        
        // compute the desired destination positions, enforcing constraints
        Vec3c original_solid = m_surf.vertex_is_solid_3(xj);
        Vec3d original_position = m_surf.get_position(xj);
        
        double mean_edge_length = 0;
        int edge_count = 0;
        for (size_t i = 0; i < mesh.m_vertex_to_edge_map[xj].size(); i++)
        {
            size_t v0 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][0];
            size_t v1 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][1];
            mean_edge_length += mag(m_surf.get_position(v1) - m_surf.get_position(v0));
            edge_count++;
        }
        assert(edge_count > 0);
        mean_edge_length /= edge_count;
        
        //        Vec3d pull_apart_offset = pull_apart_direction * mean_edge_length;
        Vec3d pull_apart_offset = pull_apart_direction;
        
        Vec3d a_desired_position = original_position + pull_apart_offset * m_pull_apart_distance;
        Vec3d b_desired_position = original_position - pull_apart_offset * m_pull_apart_distance;
        size_t a = static_cast<size_t>(~0);
        size_t b = static_cast<size_t>(~0);
        
        if (m_surf.vertex_is_any_solid(xj))
        {
            assert(m_surf.m_solid_vertices_callback);
            m_surf.m_solid_vertices_callback->generate_vertex_popped_positions(m_surf, xj, A, B, a_desired_position, b_desired_position);
        }
        
        // collision test
        // sort the incident faces and edges into those that go with nv0, and those that go with nv1 (the two groups are not necessarily disjoint)
        std::vector<size_t> A_faces;
        std::vector<size_t> A_edges;
        
        for (size_t j = 0; j < mesh.m_vertex_to_triangle_map[xj].size(); j++)
        {
            size_t triangle = mesh.m_vertex_to_triangle_map[xj][j];
            
            Vec2i label = mesh.get_triangle_label(triangle);
            if (label[0] == A || label[1] == A)
                A_faces.push_back(triangle);
        }
        
        for (size_t j = 0; j < mesh.m_vertex_to_edge_map[xj].size(); j++)
        {
            size_t edge = mesh.m_vertex_to_edge_map[xj][j];
            
            bool adjA = false;
            for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
            {
                Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][k]);
                if (label[0] == A || label[1] == A)
                    adjA = true;
            }
            
            if (adjA)
                A_edges.push_back(edge);
        }
        
        if (vertex_pseudo_motion_introduces_collision(xj, original_position, b_desired_position))
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: pulling vertex " << xj << " apart introduces collision." << std::endl;
            continue;
        }
        
        m_surf.set_position(xj, b_desired_position);
        if (vertex_pseudo_motion_introduces_collision(xj, b_desired_position, a_desired_position, A_faces, A_edges))
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: pulling vertex " << xj << " apart introduces collision." << std::endl;
            m_surf.set_position(xj, original_position);
            continue;
        }
        
        // check intersection in the final configuration
        const std::vector<Vec3d> & x = m_surf.get_positions();
        bool collision = false;
        
        // point-tet
        for (size_t j = 0; j < A_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(A_faces[j]);
            
            Vec3d low, high;
            minmax(x[t[0]], x[t[1]], x[t[2]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_vertices;
            m_surf.m_broad_phase->get_potential_vertex_collisions(low, high, true, true, overlapping_vertices);
            
            for (size_t k = 0; k < overlapping_vertices.size(); k++)
            {
                size_t ov = overlapping_vertices[k];
                if (mesh.m_vertex_to_triangle_map[ov].size() == 0)
                    continue;
                if (ov == t[0] || ov == t[2] || ov == t[1])
                    continue;
                
                if (point_tetrahedron_intersection(x[ov], ov, x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], a_desired_position, mesh.nv()))
                    collision = true;
            }
        }
        
        // edge-triangle
        for (size_t j = 0; j < A_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(A_faces[j]);
            if (t[1] == xj) std::swap(t[0], t[1]);
            if (t[2] == xj) std::swap(t[0], t[2]);
            
            Vec3d low, high;
            minmax(x[t[1]], x[t[2]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_edges;
            m_surf.m_broad_phase->get_potential_edge_collisions(low, high, true, true, overlapping_edges);
            
            for (size_t k = 0; k < overlapping_edges.size(); k++)
            {
                const Vec2st & e = mesh.m_edges[overlapping_edges[k]];
                if (e[0] == e[1])
                    continue;
                if (e[0] == t[1] || e[1] == t[1] || e[0] == t[2] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < A_edges.size(); l++)
                    if (A_edges[l] == overlapping_edges[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[0]], e[0], x[e[1]], e[1], x[t[1]], t[1], x[t[2]], t[2], a_desired_position, mesh.nv(), true))
                    collision = true;
            }
        }
        
        // triangle-edge
        for (size_t j = 0; j < A_edges.size(); j++)
        {
            Vec2st e = mesh.m_edges[A_edges[j]];
            if (e[1] == xj) std::swap(e[0], e[1]);
            
            Vec3d low, high;
            minmax(x[e[1]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_triangles;
            m_surf.m_broad_phase->get_potential_triangle_collisions(low, high, true, true, overlapping_triangles);
            
            for (size_t k = 0; k < overlapping_triangles.size(); k++)
            {
                const Vec3st & t = mesh.get_triangle(overlapping_triangles[k]);
                if (t[0] == t[1] || t[0] == t[2] || t[1] == t[2])
                    continue;
                if (e[1] == t[0] || e[1] == t[1] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < A_faces.size(); l++)
                    if (A_faces[l] == overlapping_triangles[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[1]], e[1], a_desired_position, mesh.nv(), x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], true))
                    collision = true;
            }
        }
        
        if (collision)
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: collision introduced." << std::endl;
            m_surf.set_position(xj, original_position);
            continue;
        }
        
        // all clear, now perform the pull-apart operation
        
        void * data = NULL;
        if (m_surf.m_mesheventcallback)
            m_surf.m_mesheventcallback->pre_t1(m_surf, xj, &data);

        // pull apart
        std::vector<size_t> verts_to_delete;
        std::vector<Vec3d> verts_to_create;
        std::vector<size_t> verts_created;
        
        a = m_surf.add_vertex(a_desired_position, m_surf.m_masses[xj]);
        b = m_surf.add_vertex(b_desired_position, m_surf.m_masses[xj]);
        
        m_surf.set_remesh_velocity(a, m_surf.get_remesh_velocity(xj));
        m_surf.set_remesh_velocity(b, m_surf.get_remesh_velocity(xj));
        for (int i = 0; i < 3; i++)
        {
            if (original_solid[i])
            {
                m_surf.m_masses[a][i] = std::numeric_limits<double>::infinity();
                m_surf.m_masses[b][i] = std::numeric_limits<double>::infinity();
            }
        }
        
        verts_to_delete.push_back(xj);
        verts_to_create.push_back(a_desired_position);
        verts_to_create.push_back(b_desired_position);
        verts_created.push_back(a);
        verts_created.push_back(b);
        
        // update the face connectivities
        std::vector<size_t> faces_to_delete;
        std::vector<Vec3st> faces_to_create;
        std::vector<Vec2i> face_labels_to_create;
        std::vector<size_t> faces_created;
        
        triangulate_popped_vertex(xj, A, B, a, b, faces_to_delete, faces_to_create, face_labels_to_create);
        
        if (m_surf.m_verbose)
        {
            std::cout << "Vertex " << xj << " (" << m_surf.get_position(xj) << " is splitted into vertex " << a << " (" << a_desired_position << ") and vertex " << b << " (" << b_desired_position << ")" << std::endl;
        }
        
        // apply the deletion/addition
        assert(faces_to_create.size() == face_labels_to_create.size());
        for (size_t i = 0; i < faces_to_create.size(); i++)
        {
            size_t nf = m_surf.add_triangle(faces_to_create[i], face_labels_to_create[i]);
            faces_created.push_back(nf);
        }
        
        for (size_t i = 0; i < faces_to_delete.size(); i++)
            m_surf.remove_triangle(faces_to_delete[i]);
        
        m_surf.remove_vertex(xj);
        
        // Add to new history log
        MeshUpdateEvent vertpop(MeshUpdateEvent::VERTEX_POP);
        vertpop.m_deleted_tris = faces_to_delete;
        vertpop.m_created_tris = faces_created;
        vertpop.m_created_tri_data = faces_to_create;
        vertpop.m_created_tri_labels = face_labels_to_create;
        vertpop.m_deleted_verts = verts_to_delete;
        vertpop.m_created_verts = verts_created;
        vertpop.m_created_vert_data = verts_to_create;
        m_surf.m_mesh_change_history.push_back(vertpop);
        
        pop_occurred = true;
        
        if (m_surf.m_mesheventcallback)
            m_surf.m_mesheventcallback->post_t1(m_surf, xj, a, b, data);
        
    }
    
    return pop_occurred;
}


// --------------------------------------------------------
///
/// Decide whether to cut an X-junction vertex between two given regions
///
// --------------------------------------------------------

double T1Transition::try_pull_vertex_apart_using_surface_tension(size_t xj, int A, int B, Vec3d & pull_apart_direction)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    // compute surface tension pulling force to see if this vertex pair needs to be pulled open.
    // first find the 1-ring neighbors in the cone of region A and those in the cone of region B
    std::vector<Vec3d> vertsA;
    std::vector<Vec3d> vertsB;
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t other_edge = static_cast<size_t>(~0);
        for (l = 0; l < 3; l++)
        {
            if (mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][0] != xj &&
                mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][1] != xj)
                other_edge = mesh.m_triangle_to_edge_map[triangle][l];
        }
        assert(other_edge < mesh.ne());
        size_t v0 = mesh.m_edges[other_edge][0];
        size_t v1 = mesh.m_edges[other_edge][1];
        
        Vec3d x0 = m_surf.get_position(v0);
        Vec3d x1 = m_surf.get_position(v1);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            vertsA.push_back(x0);
            vertsA.push_back(x1);
        } else if (label[0] == B || label[1] == B)
        {
            vertsB.push_back(x0);
            vertsB.push_back(x1);
        }
    }
    assert(vertsA.size() > 0);
    assert(vertsB.size() > 0);
    
    // compute the centroids of the two cones
    Vec3d centroidA(0, 0, 0);
    Vec3d centroidB(0, 0, 0);
    for (size_t i = 0; i < vertsA.size(); i++)
        centroidA += vertsA[i];
    centroidA /= vertsA.size();
    for (size_t i = 0; i < vertsB.size(); i++)
        centroidB += vertsB[i];
    centroidB /= vertsB.size();
    
    // the pull apart direction is along the line between the two centroids
    pull_apart_direction = (centroidA - centroidB);
    pull_apart_direction /= mag(pull_apart_direction);
    
    // compute the mean edge length around vertex xj
    double mean_edge_length = 0;
    int edge_count = 0;
    for (size_t i = 0; i < mesh.m_vertex_to_edge_map[xj].size(); i++)
    {
        size_t v0 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][0];
        size_t v1 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][1];
        mean_edge_length += mag(m_surf.get_position(v1) - m_surf.get_position(v0));
        edge_count++;
    }
    assert(edge_count > 0);
    mean_edge_length /= edge_count;
    
    // do a trial pull-apart
    Vec3d xxj = m_surf.get_position(xj);
    Vec3d force_a(0);
    Vec3d force_b(0);
    
    std::vector<size_t> faces_to_delete;
    std::vector<Vec3st> faces_to_create;
    std::vector<Vec2i> face_labels_to_create;
    
    size_t a = mesh.nv() + 1;   // placeholder vertices
    size_t b = mesh.nv() + 2;
    triangulate_popped_vertex(xj, A, B, a, b, faces_to_delete, faces_to_create, face_labels_to_create);
    
    // compute pre-pull-apart surface area
    double pre_area = 0;
    for (size_t i = 0; i < faces_to_delete.size(); i++)
        pre_area += m_surf.get_triangle_area(faces_to_delete[i]);
    
    // compute post-pull-apart surface area
    double post_area = 0;
    for (size_t i = 0; i < faces_to_create.size(); i++)
    {
        Vec3st & t = faces_to_create[i];
        Vec3d pos[3];
        for (int j = 0; j < 3; j++)
        {
            if (t[j] == a)
                //                pos[j] = xxj + pull_apart_direction * mean_edge_length * m_pull_apart_distance;
                pos[j] = xxj + pull_apart_direction * m_pull_apart_distance;
            else if (t[j] == b)
                //                pos[j] = xxj - pull_apart_direction * mean_edge_length * m_pull_apart_distance;
                pos[j] = xxj - pull_apart_direction * m_pull_apart_distance;
            else
                pos[j] = m_surf.get_position(t[j]);
        }
        post_area += 0.5 * mag(cross(pos[1] - pos[0], pos[2] - pos[0]));
    }
    
    double area_diff = (pre_area - post_area);
    
    return area_diff;
}

// --------------------------------------------------------
///
/// Decide whether to cut an X-junction vertex between two given regions
///
// --------------------------------------------------------

double T1Transition::try_pull_vertex_apart_using_velocity_field(size_t xj, int A, int B, Vec3d & pull_apart_direction)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    // compute surface tension pulling force to see if this vertex pair needs to be pulled open.
    // first find the 1-ring neighbors in the cone of region A and those in the cone of region B
    std::vector<Vec3d> vertsA;
    std::vector<Vec3d> vertsB;
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t other_edge = static_cast<size_t>(~0);
        for (l = 0; l < 3; l++)
        {
            if (mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][0] != xj &&
                mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][1] != xj)
                other_edge = mesh.m_triangle_to_edge_map[triangle][l];
        }
        assert(other_edge < mesh.ne());
        size_t v0 = mesh.m_edges[other_edge][0];
        size_t v1 = mesh.m_edges[other_edge][1];
        
        Vec3d x0 = m_surf.get_position(v0);
        Vec3d x1 = m_surf.get_position(v1);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            vertsA.push_back(x0);
            vertsA.push_back(x1);
        } else if (label[0] == B || label[1] == B)
        {
            vertsB.push_back(x0);
            vertsB.push_back(x1);
        }
    }
    assert(vertsA.size() > 0);
    assert(vertsB.size() > 0);
    
    // compute the centroids of the two cones
    Vec3d centroidA(0, 0, 0);
    Vec3d centroidB(0, 0, 0);
    for (size_t i = 0; i < vertsA.size(); i++)
        centroidA += vertsA[i];
    centroidA /= vertsA.size();
    for (size_t i = 0; i < vertsB.size(); i++)
        centroidB += vertsB[i];
    centroidB /= vertsB.size();
    
    // the pull apart direction is along the line between the two centroids
    pull_apart_direction = (centroidA - centroidB);
    pull_apart_direction /= mag(pull_apart_direction);
    
    // compute the mean edge length around vertex xj
    double mean_edge_length = 0;
    int edge_count = 0;
    for (size_t i = 0; i < mesh.m_vertex_to_edge_map[xj].size(); i++)
    {
        size_t v0 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][0];
        size_t v1 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][1];
        mean_edge_length += mag(m_surf.get_position(v1) - m_surf.get_position(v0));
        edge_count++;
    }
    assert(edge_count > 0);
    mean_edge_length /= edge_count;
    
    assert(m_velocity_field_callback);
    
    // decide the final positions
    Vec3d xxj = m_surf.get_position(xj);
    //    Vec3d x_a = xxj + pull_apart_direction * mean_edge_length * m_pull_apart_distance;
    //    Vec3d x_b = xxj - pull_apart_direction * mean_edge_length * m_pull_apart_distance;
    Vec3d x_a = xxj + pull_apart_direction * m_pull_apart_distance;
    Vec3d x_b = xxj - pull_apart_direction * m_pull_apart_distance;
    
    //    Vec3d vxj = m_velocity_field_callback->sampleVelocity(xxj);
    Vec3d v_a = m_velocity_field_callback->sampleVelocity(x_a);
    Vec3d v_b = m_velocity_field_callback->sampleVelocity(x_b);
    
    //    double divergence = dot(v_a - v_b, pull_apart_direction) / (mean_edge_length * m_pull_apart_distance * 2);
    double divergence = dot(v_a - v_b, pull_apart_direction) / (m_pull_apart_distance * 2);
    
    return divergence;
}

void T1Transition::triangulate_popped_vertex(size_t xj, int A, int B, size_t a, size_t b, std::vector<size_t> & faces_to_delete, std::vector<Vec3st> & faces_to_create, std::vector<Vec2i> & face_labels_to_create)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        faces_to_delete.push_back(triangle);
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t edge2 = static_cast<size_t>(~0);
        for (l = 0; l < 3; l++)
        {
            size_t e = mesh.m_triangle_to_edge_map[triangle][l];
            if (mesh.m_edges[e][0] != xj && mesh.m_edges[e][1] != xj)
                edge2 = e;
        }
        assert(edge2 < mesh.ne());
        size_t v0 = mesh.m_edges[edge2][0];
        size_t v1 = mesh.m_edges[edge2][1];
        
        if (!mesh.oriented(v0, v1, mesh.get_triangle(triangle)))
            std::swap(v0, v1);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            faces_to_create.push_back(Vec3st(a, v0, v1));
            face_labels_to_create.push_back(label);
        } else
        {
            faces_to_create.push_back(Vec3st(b, v0, v1));
            face_labels_to_create.push_back(label);
        }
    }
    
    std::vector<size_t> A_edges;
    for (size_t j = 0; j < mesh.m_vertex_to_edge_map[xj].size(); j++)
    {
        size_t edge = mesh.m_vertex_to_edge_map[xj][j];
        for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][k]);
            if (label[0] == A || label[1] == A)
            {
                A_edges.push_back(edge);
                break;
            }
        }
    }
    
    // sweep A region edges
    for (size_t i = 0; i < A_edges.size(); i++)
    {
        size_t edge = A_edges[i];
        size_t v2 = (mesh.m_edges[edge][0] == xj ? mesh.m_edges[edge][1] : mesh.m_edges[edge][0]);
        
        int upper_region = -1;  // the region on the top when looking down the edge from xj to v2, with region B on the right
        int lower_region = -1;
        for (size_t j = 0; j < mesh.m_edge_to_triangle_map[edge].size(); j++)
        {
            size_t triangle = mesh.m_edge_to_triangle_map[edge][j];
            bool oriented = mesh.oriented(xj, v2, mesh.get_triangle(triangle));
            
            Vec2i label = mesh.get_triangle_label(triangle);
            
            if ((label[0] == A &&  oriented) ||
                (label[1] == A && !oriented))
            {
                upper_region = (label[0] == A ? label[1] : label[0]);
            }
            if ((label[0] == A && !oriented) ||
                (label[1] == A &&  oriented))
            {
                lower_region = (label[0] == A ? label[1] : label[0]);
            }
        }
        
        if (upper_region >= 0 && lower_region >= 0) // if this is not true, then the neighborhood around this edge is not complete, which can oly happen on the boundary.
        {
            if (upper_region == lower_region)
            {
                // this means either this edge is just a manifold edge between region A faces (thus pulling apart xj doesn't affect this edge), or it is an X junction edge with the same region above and below (pulling this edge apart creates a tunnel connecting them)
                // in either case, no triangle should be created.
            } else
            {
                faces_to_create.push_back(Vec3st(a, b, v2));
                face_labels_to_create.push_back(Vec2i(lower_region, upper_region));
            }
        }
    }
    
}

bool T1Transition::pulling_vertex_apart_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos0, const Vec3d & newpos1)
{
    bool collision0 = vertex_pseudo_motion_introduces_collision(v, oldpos, newpos0);
    bool collision1 = vertex_pseudo_motion_introduces_collision(v, oldpos, newpos1);
    
    return collision0 || collision1;
}

bool T1Transition::vertex_pseudo_motion_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos)
{
    // code adapted from EdgeSplitter::split_edge_pseudo_motion_introduces_intersection()
    
    if (!m_surf.m_collision_safety)
        return false;
    
    NonDestructiveTriMesh & m_mesh = m_surf.m_mesh;
    
    const std::vector<Vec3d> & x = m_surf.get_positions();
    
    std::vector<size_t> & tris = m_surf.m_mesh.m_vertex_to_triangle_map[v];
    std::vector<size_t> & edges = m_surf.m_mesh.m_vertex_to_edge_map[v];
    std::vector<size_t> edge_other_endpoints(edges.size());
    
    for (size_t i = 0; i < edges.size(); i++)
        edge_other_endpoints[i] = (m_mesh.m_edges[edges[i]][0] == v ? m_mesh.m_edges[edges[i]][1] : m_mesh.m_edges[edges[i]][0]);
    
    // new point vs all triangles
    {
        
        Vec3d aabb_low, aabb_high;
        minmax(oldpos, newpos, aabb_low, aabb_high);
        
        aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_triangles;
        m_surf.m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles);
        
        for (size_t i = 0; i < overlapping_triangles.size(); i++)
        {
            // exclude incident triangles
            if (m_mesh.get_triangle(overlapping_triangles[i])[0] == v ||
                m_mesh.get_triangle(overlapping_triangles[i])[1] == v ||
                m_mesh.get_triangle(overlapping_triangles[i])[2] == v)
                continue;
            
            Vec3st sorted_triangle = sort_triangle(m_mesh.get_triangle(overlapping_triangles[i]));
            size_t a = sorted_triangle[0];
            size_t b = sorted_triangle[1];
            size_t c = sorted_triangle[2];
            
            double t_zero_distance;
            check_point_triangle_proximity(oldpos, x[a], x[b], x[c], t_zero_distance);
            if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                return true;
            
            if (point_triangle_collision(oldpos, newpos, v, x[a], x[a], a, x[b], x[b], b, x[c], x[c], c))
            {
                if (m_surf.m_verbose)
                    std::cout << "Popping collision: point triangle: with triangle " << overlapping_triangles[i] << std::endl;
                return true;
            }
        }
        
    }
    
    // new edges vs all edges
    {
        Vec3d edge_aabb_low, edge_aabb_high;
        
        // do one big query into the broad phase for all new edges
        minmax(oldpos, newpos, edge_aabb_low, edge_aabb_high);
        for (size_t i = 0; i < edge_other_endpoints.size(); ++i)
            update_minmax(m_surf.get_position(edge_other_endpoints[i]), edge_aabb_low, edge_aabb_high);
        
        edge_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        edge_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_edges;
        m_surf.m_broad_phase->get_potential_edge_collisions(edge_aabb_low, edge_aabb_high, true, true, overlapping_edges);
        
        for (size_t i = 0; i < overlapping_edges.size(); i++)
        {
            if (m_mesh.m_edges[overlapping_edges[i]][0] == m_mesh.m_edges[overlapping_edges[i]][1])
                continue;
            
            for (size_t j = 0; j < edges.size(); j++)
            {
                // exclude adjacent edges
                if (m_mesh.get_common_vertex(edges[j], overlapping_edges[i]) < m_mesh.nv())
                    continue;
                
                size_t n = edge_other_endpoints[j];
                size_t e0 = m_mesh.m_edges[overlapping_edges[i]][0];
                size_t e1 = m_mesh.m_edges[overlapping_edges[i]][1];
                if (e0 > e1)
                    std::swap(e0, e1);
                
                double t_zero_distance;
                check_edge_edge_proximity(oldpos, x[n], x[e0], x[e1], t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                bool collision = (n < v ?
                                  segment_segment_collision(x[n], x[n], n, oldpos, newpos, v, x[e0], x[e0], e0, x[e1], x[e1], e1) :
                                  segment_segment_collision(oldpos, newpos, v, x[n], x[n], n, x[e0], x[e0], e0, x[e1], x[e1], e1));
                
                if (collision)
                {
                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: edge edge: edge other vertex = " << edge_other_endpoints[j] << " edge = " << overlapping_edges[i] << std::endl;
                    return true;
                }
            }
        }
    }
    
    // new triangles vs all points
    {
        Vec3d triangle_aabb_low, triangle_aabb_high;
        
        // do one big query into the broad phase for all new triangles
        minmax(oldpos, newpos, triangle_aabb_low, triangle_aabb_high);
        for (size_t i = 0; i < edge_other_endpoints.size(); ++i)
            update_minmax(m_surf.get_position(edge_other_endpoints[i]), triangle_aabb_low, triangle_aabb_high);
        
        triangle_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        triangle_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_vertices;
        m_surf.m_broad_phase->get_potential_vertex_collisions(triangle_aabb_low, triangle_aabb_high, true, true, overlapping_vertices);
        
        for (size_t i = 0; i < overlapping_vertices.size(); i++)
        {
            if (m_mesh.m_vertex_to_triangle_map[overlapping_vertices[i]].empty())
                continue;
            
            const Vec3d & vert = m_surf.get_position(overlapping_vertices[i]);
            
            for (size_t j = 0; j < tris.size(); j++)
            {
                // exclude incident triangles
                if (m_mesh.get_triangle(tris[j])[0] == overlapping_vertices[i] ||
                    m_mesh.get_triangle(tris[j])[1] == overlapping_vertices[i] ||
                    m_mesh.get_triangle(tris[j])[2] == overlapping_vertices[i])
                    continue;
                
                Vec3st sorted_triangle = sort_triangle(m_mesh.get_triangle(tris[j]));
                size_t a = sorted_triangle[0];
                size_t b = sorted_triangle[1];
                size_t c = sorted_triangle[2];
                
                Vec3d oldxa = (a == v ? oldpos : x[a]);
                Vec3d newxa = (a == v ? newpos : x[a]);
                Vec3d oldxb = (b == v ? oldpos : x[b]);
                Vec3d newxb = (b == v ? newpos : x[b]);
                Vec3d oldxc = (c == v ? oldpos : x[c]);
                Vec3d newxc = (c == v ? newpos : x[c]);
                
                double t_zero_distance;
                check_point_triangle_proximity(vert, oldxa, oldxb, oldxc, t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                if (point_triangle_collision(vert, vert, overlapping_vertices[i], oldxa, newxa, a, oldxb, newxb, b, oldxc, newxc, c))
                {
                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: triangle point: with triangle " << tris[j] << " with vertex " << overlapping_vertices[i] << std::endl;
                    return true;
                }
            }
        }
    }
    
    return false;
}

bool T1Transition::vertex_pseudo_motion_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos, const std::vector<size_t> & tris, const std::vector<size_t> & edges)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    const std::vector<Vec3d> & x = m_surf.get_positions();
    
    // sanity check: all tris triangles and all edges edges are actually incident to v
    for (size_t i = 0; i < tris.size(); i++)
    {
        const Vec3st & t = mesh.get_triangle(tris[i]);
        assert(t[0] == v || t[1] == v || t[2] == v);
    }
    
    for (size_t i = 0; i < edges.size(); i++)
    {
        const Vec2st & e = mesh.m_edges[edges[i]];
        assert(e[0] == v || e[1] == v);
    }
    
    // new point vs all triangles
    {
        Vec3d aabb_low, aabb_high;
        minmax(oldpos, newpos, aabb_low, aabb_high);
        
        aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_triangles;
        m_surf.m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles);
        
        for (size_t i = 0; i < overlapping_triangles.size(); i++)
        {
            const Vec3st & t = mesh.get_triangle(overlapping_triangles[i]);
            
            // exclude incident triangles
            if (t[0] == v || t[1] == v || t[2] == v)
                continue;
            
            Vec3st sorted_triangle = sort_triangle(t);
            size_t a = sorted_triangle[0];
            size_t b = sorted_triangle[1];
            size_t c = sorted_triangle[2];
            
            double t_zero_distance;
            check_point_triangle_proximity(oldpos, x[a], x[b], x[c], t_zero_distance);
            if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                return true;
            
            if (point_triangle_collision(oldpos, newpos, v, x[a], x[a], a, x[b], x[b], b, x[c], x[c], c))
            {
                if (m_surf.m_verbose)
                    std::cout << "Popping collision: point triangle: with triangle " << overlapping_triangles[i] << std::endl;
                return true;
            }
        }
        
    }
    
    // new edges vs all edges
    {
        Vec3d edge_aabb_low, edge_aabb_high;
        
        // do one big query into the broad phase for all new edges
        minmax(oldpos, newpos, edge_aabb_low, edge_aabb_high);
        for (size_t i = 0; i < edges.size(); ++i)
        {
            size_t other_endpoint = (mesh.m_edges[edges[i]][0] == v ? mesh.m_edges[edges[i]][1] : mesh.m_edges[edges[i]][0]);
            update_minmax(m_surf.get_position(other_endpoint), edge_aabb_low, edge_aabb_high);
        }
        
        edge_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        edge_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_edges;
        m_surf.m_broad_phase->get_potential_edge_collisions(edge_aabb_low, edge_aabb_high, true, true, overlapping_edges);
        
        for (size_t i = 0; i < overlapping_edges.size(); i++)
        {
            const Vec2st & e = mesh.m_edges[overlapping_edges[i]];
            
            if (e[0] == e[1])
                continue;
            
            for (size_t j = 0; j < edges.size(); j++)
            {
                // exclude adjacent edges
                if (mesh.get_common_vertex(edges[j], overlapping_edges[i]) < mesh.nv())
                    continue;
                
                size_t n = (mesh.m_edges[edges[j]][0] == v ? mesh.m_edges[edges[j]][1] : mesh.m_edges[edges[j]][0]);
                size_t e0 = e[0];
                size_t e1 = e[1];
                if (e0 > e1)
                    std::swap(e0, e1);
                
                double t_zero_distance;
                check_edge_edge_proximity(oldpos, x[n], x[e0], x[e1], t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                bool collision = (n < v ?
                                  segment_segment_collision(x[n], x[n], n, oldpos, newpos, v, x[e0], x[e0], e0, x[e1], x[e1], e1) :
                                  segment_segment_collision(oldpos, newpos, v, x[n], x[n], n, x[e0], x[e0], e0, x[e1], x[e1], e1));
                
                if (collision)
                {
                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: edge edge: edge other vertex = " << n << " edge = " << overlapping_edges[i] << std::endl;
                    return true;
                }
            }
        }      
    }
    
    // new triangles vs all points
    {
        Vec3d triangle_aabb_low, triangle_aabb_high;
        
        // do one big query into the broad phase for all new triangles
        minmax(oldpos, newpos, triangle_aabb_low, triangle_aabb_high);
        for (size_t i = 0; i < tris.size(); ++i)
        {
            const Vec3st & t = mesh.get_triangle(tris[i]);
            Vec2st other_vertices;
            if (t[0] == v)
                other_vertices = Vec2st(t[1], t[2]);
            else if (t[1] == v)
                other_vertices = Vec2st(t[2], t[0]);
            else if (t[2] == v)
                other_vertices = Vec2st(t[0], t[1]);
            else
                assert(!"triangle in tris does not contain vertex v");
            update_minmax(m_surf.get_position(other_vertices[0]), triangle_aabb_low, triangle_aabb_high);
            update_minmax(m_surf.get_position(other_vertices[1]), triangle_aabb_low, triangle_aabb_high);
        }
        
        triangle_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        triangle_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_vertices;
        m_surf.m_broad_phase->get_potential_vertex_collisions(triangle_aabb_low, triangle_aabb_high, true, true, overlapping_vertices);
        
        for (size_t i = 0; i < overlapping_vertices.size(); i++)
        {
            size_t & ov = overlapping_vertices[i];
            
            if (mesh.m_vertex_to_triangle_map[ov].empty()) 
                continue; 
            
            const Vec3d & vert = m_surf.get_position(ov);
            
            for (size_t j = 0; j < tris.size(); j++)
            {
                const Vec3st & t = mesh.get_triangle(tris[j]);
                
                // exclude incident triangles
                if (t[0] == ov || t[1] == ov || t[2] == ov)
                    continue;
                
                Vec3st sorted_triangle = sort_triangle(t);
                size_t a = sorted_triangle[0];
                size_t b = sorted_triangle[1];
                size_t c = sorted_triangle[2];
                
                Vec3d oldxa = (a == v ? oldpos : x[a]);
                Vec3d newxa = (a == v ? newpos : x[a]);
                Vec3d oldxb = (b == v ? oldpos : x[b]);
                Vec3d newxb = (b == v ? newpos : x[b]);
                Vec3d oldxc = (c == v ? oldpos : x[c]);
                Vec3d newxc = (c == v ? newpos : x[c]);
                
                double t_zero_distance;
                check_point_triangle_proximity(vert, oldxa, oldxb, oldxc, t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                if (point_triangle_collision(vert, vert, overlapping_vertices[i], oldxa, newxa, a, oldxb, newxb, b, oldxc, newxc, c))
                {
                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: triangle point: with triangle " << tris[j] << " with vertex " << ov << std::endl;
                    return true;
                }
            }
        }
    }
    
    return false;
    
    
    
}

    
    
}
