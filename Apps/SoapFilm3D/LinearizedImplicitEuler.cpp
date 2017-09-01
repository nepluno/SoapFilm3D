#include "LinearizedImplicitEuler.h"
#include "VS3D.h"

// TODO: Save space by creating dx, dv, rhs, A only once.

LinearizedImplicitEuler::LinearizedImplicitEuler()
: SceneStepper()
{}

LinearizedImplicitEuler::~LinearizedImplicitEuler()
{}

bool LinearizedImplicitEuler::stepScene( VS3D& scene, scalar dt )
{
    const int ndof = scene.constrainedPositions().size() * 3;
    const VecXd& x = Eigen::Map<VecXd>((double*) &scene.constrainedPositions()[0], ndof);
    const VecXd& v = Eigen::Map<VecXd>((double*) &scene.constrainedVelocities()[0], ndof);
    
    const VecXd& m = Eigen::Map<VecXd>((double*) &scene.constrainedMass()[0], ndof);
    
    assert(x.size() == v.size());
    assert( ndof % 3 == 0 );
    int nprts = scene.constrainedPositions().size();
    assert( 3 * nprts == ndof );
    
    // We specify the current iterate as a change from last timestep's solution
    // Linearizing about the previous timestep's solution implies dx = dt*v0
    VectorXs dx = dt*v;
    // Linearizing about the previous timestep's solution implies dv = 0
    VectorXs dv = VectorXs::Zero(ndof);
    
    scene.preCompute(dx, dv, dt);
    // RHS of linear system is force. rhs == f

    SparseXs A(ndof,ndof);
    TripletXs tri_A;
    
    SparseXs Av(ndof, ndof);
    TripletXs tri_Av;
    
    SparseXs M(ndof, ndof);
    TripletXs mm;
    
    SparseXs C(ndof, ndof);
    TripletXs cc;
    
    for(int i = 0; i < ndof; ++i) {
        if(scene.constrainedFixed()[i / 3]) {
            mm.push_back(Triplets(i, i, 1.0));
        } else {
            mm.push_back(Triplets(i, i, m(i)));
            cc.push_back(Triplets(i, i, 1.0));
        }
    }
    M.setFromTriplets(mm.begin(), mm.end());
    C.setFromTriplets(cc.begin(), cc.end());
    
    tri_A.clear();
    scene.accumulateddUdxdx(tri_A);
    
    tri_Av.clear();
    scene.accumulateddUdxdv(tri_Av);
    
    // lhs == -df/dx
    A.setFromTriplets(tri_A.begin(), tri_A.end());
    
    A *= dt;
    
    SparseXs B = A;
    // lhs == -h*df/dx
    
    Av.setFromTriplets(tri_Av.begin(), tri_Av.end());
    // lhs == -h*df/dv -h^2*df/dx
    // lhs == -df/dv -h*df/dx
    A += Av;
    
    // lhs == -h*df/dv -h^2*df/dx
    A *= dt;
    // For scripted DoFs, zero out the rows/cols of fixed degrees of freedom
    A = C * A;
    
    // lhs == M -h*df/dv -h^2*df/dx
    A += M;
    
    Eigen::SimplicialLDLT<SparseXs> solver(A);
    
    VectorXs rhs = VectorXs::Zero(ndof);
    scene.accumulateGradU(rhs,dx,dv);
    rhs *= -1.0;
    
    
    // rhs == f - hessE * v
    rhs -= C * (B * v);
    
    // rhs == h*f
    rhs *= dt;
    // For scripted DoFs, zero the elements of fixed degrees of freedom
    zeroFixedDoFs( scene, rhs );
    
    // Change in velocity returned by linear solve
    VectorXs dqdot = solver.solve(rhs);
    
    Eigen::Map<VecXd>((double*) &scene.constrainedVelocities()[0], ndof) += dqdot;
    Eigen::Map<VecXd>((double*) &scene.constrainedPositions()[0], ndof) += dt * Eigen::Map<VecXd>((double*) &scene.constrainedVelocities()[0], ndof);
    
    return true;
}

std::string LinearizedImplicitEuler::getName() const
{
    return "Linearized Implicit Euler";
}

void LinearizedImplicitEuler::zeroFixedDoFs( const VS3D& scene, VectorXs& vec )
{
    int nprts = scene.constrainedPositions().size();
    for( int i = 0; i < nprts; ++i ) if( scene.constrainedFixed()[i] ) vec.segment<3>(3 * i).setZero();
}

