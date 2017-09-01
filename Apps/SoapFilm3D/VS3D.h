//
//  VS3D.h
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#ifndef __MultiTracker__VS3D__
#define __MultiTracker__VS3D__

#include <iostream>
#include <unordered_set>
#include "eigenheaders.h"
#include "surftrack.h"
#include "MathDefs.h"
#include "Force.h"
#include "SceneStepper.h"

class Sim;
class Scenes;

class VS3D : public LosTopos::SurfTrack::SolidVerticesCallback, public LosTopos::T1Transition::VelocityFieldCallback, public LosTopos::SurfTrack::MeshEventCallback
{
    friend class Sim;
    friend class Scenes;
    friend VecXd BiotSavart(VS3D & vs, const VecXd & dx);
    friend VecXd BiotSavart_naive(VS3D & vs, const VecXd & dx);
    friend VecXd BiotSavart_fmmtl(VS3D & vs, const VecXd & dx);
    
public:
    VS3D(const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<size_t> & constrained_vertices = std::vector<size_t>(),  const std::vector<Vec3d> & constrained_positions = std::vector<Vec3d>(), const std::vector<Vec3d> & constrained_velocities = std::vector<Vec3d>(), const std::vector<unsigned char> & constrained_fixed = std::vector<unsigned char>());
    ~VS3D();
    
    class SimOptions
    {
    public:
        bool implicit;
        bool pbd;
        bool looped;
        double smoothing_coef;
        double damping_coef;
        double sigma;
        double gravity;
        double radius;
        double density;
        double stretching;
        double bending;

        SimOptions() : implicit(false), pbd(false), smoothing_coef(0), damping_coef(1), sigma(1), gravity(0)
        { }
    };
    
    SimOptions & simOptions() { return m_sim_options; }
    
public:
    const LosTopos::SurfTrack * surfTrack() const { return m_st; }
          LosTopos::SurfTrack * surfTrack()       { return m_st; }
    const LosTopos::NonDestructiveTriMesh & mesh() const { return m_st->m_mesh; }
          LosTopos::NonDestructiveTriMesh & mesh()       { return m_st->m_mesh; }
    
    Vec3d pos(size_t v) const { return vc(m_st->pm_positions[v]); }
    
    double step(double dt);
    
    void update_dbg_quantities();
    
    const std::vector<size_t> & constrainedVertices() const { return m_constrained_vertices; }
          std::vector<size_t> & constrainedVertices()       { return m_constrained_vertices; }
    const std::vector<Vec3d> & constrainedPositions() const { return m_constrained_positions; }
          std::vector<Vec3d> & constrainedPositions()       { return m_constrained_positions; }
    const std::vector<Vec3d> & constrainedVelocities() const { return m_constrained_velocities; }
          std::vector<Vec3d> & constrainedVelocities()       { return m_constrained_velocities; }
    const std::vector<double> & constrainedMass() const { return m_constrained_mass; }
          std::vector<double> & constrainedMass()       { return m_constrained_mass; }
    const std::vector<unsigned char> & constrainedFixed() const { return m_constrained_fixed; }
          std::vector<unsigned char> & constrainedFixed()       { return m_constrained_fixed; }
    
    void accumulateGradU( VectorXs& F, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs() );
    
    void accumulateddUdxdx( TripletXs& A, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs() );
    
    // Kind of a misnomer.
    void accumulateddUdxdv( TripletXs& A, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs() );
    
    void preCompute( const VectorXs& dx, const VectorXs& dv, const scalar& dt );
    
    void stepConstrainted( const scalar& dt );
    
    double delta() const { return m_delta; }
    
    int nregion() const { return m_nregion; }
    
    bool isVertexConstrained( size_t vert );
public:
    class GammaType
    {
    public:
        GammaType() { }
        GammaType(int nregion) { values = MatXd::Zero(nregion, nregion); }
        void setZero() { values.setZero(); }
        
        void   set(const Vec2i & l,           double v) { assert_valid(l[0], l[1]); (l[0] < l[1] ? values(l[0], l[1]) : values(l[1], l[0])) = v * (l[0] < l[1] ? 1 : -1); }
        void   set(const LosTopos::Vec2i & l, double v) { assert_valid(l[0], l[1]); (l[0] < l[1] ? values(l[0], l[1]) : values(l[1], l[0])) = v * (l[0] < l[1] ? 1 : -1); }
        void   set(int l0, int l1,            double v) { assert_valid(l0,   l1);   (l0   < l1   ? values(l0,   l1)   : values(l1,   l0))   = v * (l0   < l1   ? 1 : -1); }
        
        double get(const Vec2i & l) const               { assert_valid(l[0], l[1]); return (l[0] < l[1] ? values(l[0], l[1]) : -values(l[1], l[0])); }
        double get(const LosTopos::Vec2i & l) const     { assert_valid(l[0], l[1]); return (l[0] < l[1] ? values(l[0], l[1]) : -values(l[1], l[0])); }
        double get(int l0, int l1)  const               { assert_valid(l0,   l1);   return (l0   < l1   ? values(l0,   l1)   : -values(l1,   l0)); }

        void assert_valid(int l0, int l1) const { assert(values.cols() == values.rows()); assert(l0 >= 0); assert(l1 >= 0); assert(l0 < values.cols()); assert(l1 < values.rows()); }
        
    public:
        MatXd values;   // only uses the upper triangular half (values(i,j), i < j). values(i,j) satisfies that the tangential velocity jump delta u = u(i) - u(j) = grad values(i,j)
    };

    const GammaType & Gamma(size_t v) const { return (*m_Gamma)[v]; }
          GammaType & Gamma(size_t v)       { return (*m_Gamma)[v]; }
    
protected:
    void step_explicit(double dt);
    void step_implicit(double dt);
    void step_PBD_implicit(double dt);
    
protected:
    // SurfTrack::SolidVerticesCallback method
    bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
    bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
    bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
    
    // T1Transition::VelocityFieldCallback methods
    LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
    bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);

    // SurfTrack::MeshEventCallback
    void pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data);
    
    void pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data);
    
    void pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_flip(const LosTopos::SurfTrack & st, size_t e, void * data);
    
    void pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data);
    void post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data);
    
    void pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data);
    void post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data);
    
    void pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data);
    void post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data);
    
    void pre_smoothing(const LosTopos::SurfTrack & st, void ** data);
    void post_smoothing(const LosTopos::SurfTrack & st, void * data);

    std::ostream & log() { static std::stringstream ss; return ss; }

protected:
    LosTopos::SurfTrack * m_st;

    SimOptions m_sim_options;
    double m_delta;     // Biot-Savart regularization parameter
    
    // sheet internal dynamics
    int m_nregion;
    LosTopos::NonDestructiveTriMesh::VertexData<GammaType> * m_Gamma;     // average circulation of a vertex \Gamma (one scalar value for each region pair incident to the vertex)

    std::vector<Vec3d> m_dbg_t1;
    std::vector<Vec3d> m_dbg_t2;
    std::vector<std::vector<double> > m_dbg_e1;
    std::vector<std::vector<double> > m_dbg_v1;
    std::vector<Vec3d> m_dbg_v2;
    
    // constrained vertices
    std::vector<size_t> m_constrained_vertices;
    std::vector<Vec3d> m_constrained_positions;
    std::vector<Vec3d> m_constrained_velocities;
    std::vector<double> m_constrained_mass;
    std::vector<unsigned char> m_constrained_fixed;
    std::unordered_set<size_t> m_constrained_mapping;
    
    std::vector< Force* > m_forces;
    
    std::vector<Vec3d> m_obefc;     // open boundary extra face centers
    std::vector<Vec3d> m_obefe;     // open boundary extra face vorticity direction (i.e. edge tangent)
    std::vector<Vec3d> m_obefn;     // open boundary extra face normals
    std::vector<double> m_obefv;    // open boundary extra face vorticity magnitudes
    
    SceneStepper* m_constraint_stepper;
};

#endif /* defined(__MultiTracker__VS3D__) */
