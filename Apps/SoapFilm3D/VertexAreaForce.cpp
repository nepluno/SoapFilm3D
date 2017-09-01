#include "VertexAreaForce.h"
#include "VS3D.h"

using namespace LosTopos;

VertexAreaForce::VertexAreaForce( VS3D* parent, const scalar& k )
: Force()
, m_parent(parent)
, m_stiffness(k)
{
}

VertexAreaForce::~VertexAreaForce()
{}

void VertexAreaForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%3 == 0 );
}

inline void compute_gradE(const Vector3s& x0, const Vector3s& x1, const Vector3s& x2, Vector3s& gradE)
{
    const scalar t2 = x0(0)-x1(0);
    const scalar t3 = x0(0)-x2(0);
    const scalar t4 = x0(1)-x2(1);
    const scalar t5 = t2*t4;
    const scalar t6 = x0(1)-x1(1);
    const scalar t13 = t3*t6;
    const scalar t7 = t5-t13;
    const scalar t8 = x0(2)-x2(2);
    const scalar t9 = t2*t8;
    const scalar t10 = x0(2)-x1(2);
    const scalar t18 = t3*t10;
    const scalar t11 = t9-t18;
    const scalar t15 = t6*t8;
    const scalar t16 = t4*t10;
    const scalar t12 = -t15+t16;
    const scalar t14 = x1(2)-x2(2);
    const scalar t17 = t7*t7;
    const scalar t19 = t11*t11;
    const scalar t20 = t15-t16;
    const scalar t21 = x1(0)-x2(0);
    const scalar t22 = x1(1)-x2(1);
    const scalar t23 = t20*t20;
    const scalar t24 = t17+t19+t23;
    if(t24 <= 1e-24) return;
    const scalar t25 = 1.0/sqrt(t24);
    gradE(0) += (t11*t14*2.0+t7*t22*2.0)*1.0/sqrt(t17+t19+t12*t12)*(1.0/2.0);
    gradE(1) += t25*(t7*t21*2.0-t14*t20*2.0)*(-1.0/2.0);
    gradE(2) += t25*(t11*t21*2.0+t20*t22*2.0)*(-1.0/2.0);
}

inline void compute_hessE(const Vector3s& x0, const Vector3s& x1, const Vector3s& x2, Matrix3s& hessE)
{
    const scalar t2 = x1(1)-x2(1);
    const scalar t3 = x1(2)-x2(2);
    const scalar t5 = x0(0)-x1(0);
    const scalar t6 = x0(0)-x2(0);
    const scalar t8 = x0(1)-x1(1);
    const scalar t11 = x0(1)-x2(1);
    const scalar t13 = t5*t11;
    const scalar t14 = t6*t8;
    const scalar t4 = t13-t14;
    const scalar t9 = x0(2)-x2(2);
    const scalar t10 = x0(2)-x1(2);
    const scalar t15 = t5*t9;
    const scalar t16 = t6*t10;
    const scalar t7 = -t15+t16;
    const scalar t20 = t8*t9;
    const scalar t21 = t10*t11;
    const scalar t12 = t20-t21;
    const scalar t19 = t15-t16;
    const scalar t25 = t2*t4*2.0;
    const scalar t26 = t3*t19*2.0;
    const scalar t17 = t25+t26;
    const scalar t18 = t4*t4;
    const scalar t22 = t12*t12;
    const scalar t23 = t19*t19;
    const scalar t24 = t18+t22+t23;
    if(t24 <= 1e-24) return;
    const scalar t27 = x1(0)-x2(0);
    const scalar t28 = 1.0/pow(t24,3.0/2.0);
    const scalar t29 = 1.0/sqrt(t24);
    const scalar t30 = t3*t12*2.0;
    const scalar t35 = t4*t27*2.0;
    const scalar t31 = t30-t35;
    const scalar t32 = -t2*t27*t29-t17*t28*t31*(1.0/4.0);
    const scalar t33 = t3*t3;
    const scalar t34 = t33*2.0;
    const scalar t36 = t30-t35;
    const scalar t37 = t19*t27*2.0;
    const scalar t38 = t2*t12*2.0;
    const scalar t39 = t37+t38;
    const scalar t40 = t17*t28*t39*(1.0/4.0);
    const scalar t41 = t40-t3*t27*t29;
    const scalar t42 = t28*t39*(t30-t35)*(1.0/4.0);
    const scalar t43 = t42-t2*t3*t29;
    const scalar t44 = t27*t27;
    const scalar t45 = t44*2.0;
    const scalar t46 = t2*t2;
    const scalar t47 = t46*2.0;
    if(t18+t22+t7*t7 <= 1e-24) return;
    hessE(0, 0) += (t34+t47)*1.0/sqrt(t18+t22+t7*t7)*(1.0/2.0)-(t17*t17)*t28*(1.0/4.0);
    hessE(0, 1) += t32;
    hessE(0, 2) += t41;
    hessE(1, 0) += t32;
    hessE(1, 1) += t29*(t34+t45)*(1.0/2.0)-t28*(t36*t36)*(1.0/4.0);
    hessE(1, 2) += t43;
    hessE(2, 0) += t41;
    hessE(2, 1) += t43;
    hessE(2, 2) += t29*(t45+t47)*(1.0/2.0)-t28*(t39*t39)*(1.0/4.0);
}


void VertexAreaForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%3 == 0 );
    
    const std::vector<size_t>& verts = m_parent->constrainedVertices();
    for(size_t i = 0; i < verts.size(); ++i)
    {
        int pidx = verts[i];
        const std::vector<size_t>& tris = m_parent->mesh().m_vertex_to_triangle_map[pidx];
        for(size_t tri_idx : tris)
        {
            const Vec3st& tri = m_parent->mesh().m_tris[tri_idx];
            int i1, i2;
            if(tri[0] == pidx) i1 = tri[1], i2 = tri[2];
            else if(tri[1] == pidx) i1 = tri[0], i2 = tri[2];
            else i1 = tri[0], i2 = tri[1];
            
            const Vector3s& x0 = x.segment<3>(i * 3);
            const Vector3s& x1 = m_parent->pos(i1);
            const Vector3s& x2 = m_parent->pos(i2);
            
            Vector3s ge = Vector3s::Zero();
            compute_gradE(x0, x1, x2, ge);
            gradE.segment<3>(i * 3) += m_stiffness * ge;
        }
    }
}

void VertexAreaForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%3 == 0 );
    const std::vector<size_t>& verts = m_parent->constrainedVertices();
    for(size_t i = 0; i < verts.size(); ++i)
    {
        int pidx = verts[i];
        const std::vector<size_t>& tris = m_parent->mesh().m_vertex_to_triangle_map[pidx];
        for(size_t tri_idx : tris)
        {
            const Vec3st& tri = m_parent->mesh().m_tris[tri_idx];
            int i1, i2;
            if(tri[0] == pidx) i1 = tri[1], i2 = tri[2];
            else if(tri[1] == pidx) i1 = tri[0], i2 = tri[2];
            else i1 = tri[0], i2 = tri[1];
            
            const Vector3s& x0 = x.segment<3>(i * 3);
            const Vector3s& x1 = m_parent->pos(i1);
            const Vector3s& x2 = m_parent->pos(i2);
            
            Matrix3s he = Matrix3s::Zero();
            compute_hessE(x0, x1, x2, he);
            for(int r = 0; r < 3; ++r) {
                for(int s = 0; s < 3; ++s) {
                    if(he(r, s) != 0.0)
                        hessE.push_back(Triplets(i * 3 + r, i * 3 + s, m_stiffness * he(r, s)));
                }
            }
        }
    }
}

void VertexAreaForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%3 == 0 );
    // Nothing to do.
}

Force* VertexAreaForce::createNewCopy()
{
    return new VertexAreaForce(*this);
}
