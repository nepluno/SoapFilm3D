//
//  VS3DImplicit.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/27.
//
//

#include "VS3D.h"

namespace
{
    Mat3d skewSymmetric(const Vec3d & v)
    {
        Mat3d ss = Mat3d::Zero();
        ss(0, 1) = -v(2);
        ss(1, 0) =  v(2);
        ss(0, 2) =  v(1);
        ss(2, 0) = -v(1);
        ss(1, 2) = -v(0);
        ss(2, 1) =  v(0);
        return ss;
    }

    double angleAroundAxis(const Vec3d & v0, const Vec3d & v1, const Vec3d & a)   // angle from v0 to v1 around axis a
    {
        double asq = a.squaredNorm();
        assert(asq != 0);
        
        Vec3d u = v0 - v0.dot(a) * a / asq;
        assert(u.squaredNorm() != 0);
        u.normalize();
        
        Vec3d v = a.cross(u).normalized();
        
        return atan2(v1.dot(v), v1.dot(u));
    }
    
}

namespace
{
    void compute_dvdGamma(MatXd & dvdGamma, const VecXd & xs, VS3D & vs, const std::vector<std::pair<size_t, Vec2i> > & Gamma_map, std::vector<int> & Gamma_map_inv, std::vector<int> & rp_count)
    {
        size_t nv = vs.mesh().nv();
        size_t nt = vs.mesh().nt();
        
        size_t ndof = Gamma_map.size();
        
        dvdGamma = MatXd::Zero(nv * 3, ndof);
        for (size_t i = 0; i < nv; i++)
        {
            Vec3d x = xs.segment<3>(i * 3);
            
            for (size_t j = 0; j < nt; j++)
            {
                LosTopos::Vec2i l = vs.mesh().get_triangle_label(j);
                Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
                
                LosTopos::Vec3st t = vs.mesh().get_triangle(j);
                Vec3d xp = (xs.segment<3>(t[0] * 3) + xs.segment<3>(t[1] * 3) + xs.segment<3>(t[2] * 3)) / 3;
                
                Vec3d e01 = xs.segment<3>(t[1] * 3) - xs.segment<3>(t[0] * 3);
                Vec3d e12 = xs.segment<3>(t[2] * 3) - xs.segment<3>(t[1] * 3);
                Vec3d e20 = xs.segment<3>(t[0] * 3) - xs.segment<3>(t[2] * 3);
                
//                Vec3d gamma = e01 * (*m_Gamma)[t[2]].get(l) +
//                              e12 * (*m_Gamma)[t[0]].get(l) +
//                              e20 * (*m_Gamma)[t[1]].get(l);
                
                Vec3d dx = x - xp;
                
//                double dxn = dx.norm();
                double dxn = sqrt(dx.squaredNorm() + vs.delta() * vs.delta());
                
//                v += gamma.cross(dx) / (dxn * dxn * dxn);
////                v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 - exp(-dxn / m_delta));
                
                Mat3d dx_ss = skewSymmetric(dx); // skew symmetric matrix corresponding to dx
                
                Mat3d dvdgamma = -dx_ss / (dxn * dxn * dxn);
//                Mat3d dvdgamma = -dx_ss / (dxn * dxn * dxn) * (1 - exp(-dxn / m_delta));
                
                int ti0 = -1;
                for (int k = 0; k < rp_count[t[0]]; k++)
                    if (Gamma_map[Gamma_map_inv[t[0]] + k].second == rp)
                        ti0 = k;
                assert(ti0 >= 0);
                int ti1 = -1;
                for (int k = 0; k < rp_count[t[1]]; k++)
                    if (Gamma_map[Gamma_map_inv[t[1]] + k].second == rp)
                        ti1 = k;
                assert(ti1 >= 0);
                int ti2 = -1;
                for (int k = 0; k < rp_count[t[2]]; k++)
                    if (Gamma_map[Gamma_map_inv[t[2]] + k].second == rp)
                        ti2 = k;
                assert(ti2 >= 0);
                
                dvdGamma.block<3, 1>(i * 3, Gamma_map_inv[t[0]] + ti0) += -dvdgamma * e12 * (l[0] < l[1] ? 1 : -1);
                dvdGamma.block<3, 1>(i * 3, Gamma_map_inv[t[1]] + ti1) += -dvdgamma * e20 * (l[0] < l[1] ? 1 : -1);
                dvdGamma.block<3, 1>(i * 3, Gamma_map_inv[t[2]] + ti2) += -dvdgamma * e01 * (l[0] < l[1] ? 1 : -1);
            }
        }
        dvdGamma /= (4 * M_PI);
    }
    
    void compute_dHdx(VecXd & H, MatXd & dHdx, const VecXd & xs, VS3D & vs, const std::vector<std::pair<size_t, Vec2i> > & Gamma_map, std::vector<int> & Gamma_map_inv, std::vector<int> & rp_count)
    {
        size_t nv = vs.mesh().nv();
        size_t ne = vs.mesh().ne();
        size_t nt = vs.mesh().nt();
        
        size_t ndof = Gamma_map.size();
        int nregion = vs.nregion();

        H = VecXd::Zero(ndof);
        dHdx = MatXd::Zero(ndof, nv * 3);

        for (int region = 0; region < nregion; region++)
        {
            // compute the vertex areas and their derivatives first (for this region only)
            VecXd vertexAreas = VecXd::Zero(nv);
            MatXd dadx = MatXd::Zero(nv, nv * 3);
            for (size_t i = 0; i < nt; i++)
            {
                LosTopos::Vec3st t = vs.mesh().get_triangle(i);
                LosTopos::Vec2i l = vs.mesh().get_triangle_label(i);
                if (l[0] == region || l[1] == region)
                {
                    Vec3d x0 = xs.segment<3>(t[0] * 3);
                    Vec3d x1 = xs.segment<3>(t[1] * 3);
                    Vec3d x2 = xs.segment<3>(t[2] * 3);

                    Vec3d n = (x1 - x0).cross(x2 - x0);
                    double a = n.norm() / 2;
                    n.normalize();

                    //&&&& TODO: correct for the different area on triple junctions
                    for (int j = 0; j < 3; j++)
                    {
                        vertexAreas[t[j]] += a / 3;
                        dadx.block<1, 3>(t[j], t[0] * 3) += n.cross(x2 - x1) / 2 / 3;
                        dadx.block<1, 3>(t[j], t[1] * 3) += n.cross(x0 - x2) / 2 / 3;
                        dadx.block<1, 3>(t[j], t[2] * 3) += n.cross(x1 - x0) / 2 / 3;
                    }
                }
            }

            // compute edge curvatures and accumulate them into dHdx
            for (size_t i = 0; i < ne; i++)
            {
                std::vector<size_t> incident_faces; // faces incident to edge i that have the label of interest (assume there are only two of them for now; this can be false only when complex collision prevents immediate T1 resolution, which is not expected to happen for bubble complexes.)
                std::set<int> incident_regions;
                for (size_t j = 0; j < vs.mesh().m_edge_to_triangle_map[i].size(); j++)
                {
                    const LosTopos::Vec2i & l = vs.mesh().get_triangle_label(vs.mesh().m_edge_to_triangle_map[i][j]);
                    if (l[0] == region || l[1] == region)
                        incident_faces.push_back(j);
                    incident_regions.insert(l[0]);
                    incident_regions.insert(l[1]);
                }
                if (incident_faces.size() == 0)
                    continue;

                if (incident_faces.size() != 2)
                {
//                    std::cout << "Warning: incident_faces.size() != 2" << std::endl;
                    continue;
                }

                int v0 = vs.mesh().m_edges[i][0];
                int v1 = vs.mesh().m_edges[i][1];
                Vec3d x0 = xs.segment<3>(v0 * 3);
                Vec3d x1 = xs.segment<3>(v1 * 3);
                Vec3d et = (x1 - x0);

                int ti0 = vs.mesh().m_edge_to_triangle_map[i][0];
                int ti1 = vs.mesh().m_edge_to_triangle_map[i][1];
                LosTopos::Vec3st t0 = vs.mesh().get_triangle(ti0);
                LosTopos::Vec3st t1 = vs.mesh().get_triangle(ti1);

                if (vs.mesh().get_triangle_label(ti0)[vs.mesh().oriented(v0, v1, t0) ? 1 : 0] == region)
                    std::swap(ti0, ti1),
                    std::swap(t0, t1);  // the region of interest should be to the CCW direction of t0 when looking along the direciton of edge i

                int vo0 = vs.mesh().get_third_vertex(v0, v1, t0);
                int vo1 = vs.mesh().get_third_vertex(v0, v1, t1);
                Vec3d xo0 = xs.segment<3>(vo0 * 3);
                Vec3d xo1 = xs.segment<3>(vo1 * 3);

                assert(vs.mesh().get_triangle_label(ti0)[vs.mesh().oriented(v0, v1, t0) ? 0 : 1] == region);
                assert(vs.mesh().get_triangle_label(ti1)[vs.mesh().oriented(v0, v1, t1) ? 1 : 0] == region);

                Vec3d n0 = (xs.segment<3>(t0[1] * 3) - xs.segment<3>(t0[0] * 3)).cross(xs.segment<3>(t0[2] * 3) - xs.segment<3>(t0[0] * 3));
                if (vs.mesh().get_triangle_label(ti0)[1] == region)
                    n0 = -n0;   // n0 should point away from the region of interest

                Vec3d n1 = (xs.segment<3>(t1[1] * 3) - xs.segment<3>(t1[0] * 3)).cross(xs.segment<3>(t1[2] * 3) - xs.segment<3>(t1[0] * 3));
                if (vs.mesh().get_triangle_label(ti1)[1] == region)
                    n1 = -n1;   // n1 should point away from the region of interest

                double a0 = n0.norm() / 2;
                double a1 = n1.norm() / 2;
                n0.normalize();
                n1.normalize();

                // integral curvature on edge
                double kappa_i = n0.cross(n1).dot(et);  // the infinitesimal curvature
//                double kappa_i = angleAroundAxis(n0, n1, et) * et.norm(); // the finite curvature; it's what step_explicit() used but not used here

                Mat3d dn0dxo0 =  et.cross(n0) * n0.transpose() / (a0 * 2);
                Mat3d dn1dxo1 = -et.cross(n1) * n1.transpose() / (a1 * 2);

                Mat3d dn0dx0 = (xo0 - x1).cross(n0) * n0.transpose() / (a0 * 2);
                Mat3d dn0dx1 = (x0 - xo0).cross(n0) * n0.transpose() / (a0 * 2);
                Mat3d dn1dx0 = (x1 - xo1).cross(n1) * n1.transpose() / (a1 * 2);
                Mat3d dn1dx1 = (xo1 - x0).cross(n1) * n1.transpose() / (a1 * 2);

                Mat3d n0ss = skewSymmetric(n0);
                Mat3d n1ss = skewSymmetric(n1);

                // derivatives of edge (integral) curvatures
                Eigen::Matrix<double, 1, 3> dcurvaturedx0 =  (-n0.cross(n1).transpose() + et.transpose() * n0ss * dn1dx0  - et.transpose() * n1ss * dn0dx0 );
                Eigen::Matrix<double, 1, 3> dcurvaturedx1 =  ( n0.cross(n1).transpose() + et.transpose() * n0ss * dn1dx1  - et.transpose() * n1ss * dn0dx1 );
                Eigen::Matrix<double, 1, 3> dcurvaturedxo0 = (                                                            - et.transpose() * n1ss * dn0dxo0);
                Eigen::Matrix<double, 1, 3> dcurvaturedxo1 = (                            et.transpose() * n0ss * dn1dxo1                                  );

                // derivatives of vertex (pointwise) curvatures
                Eigen::Matrix<double, 1, 3> dcurvature0dx0 =  dcurvaturedx0  / vertexAreas[v0] - kappa_i / (vertexAreas[v0] * vertexAreas[v0]) * dadx.block<1, 3>(v0, v0  * 3);
                Eigen::Matrix<double, 1, 3> dcurvature0dx1 =  dcurvaturedx1  / vertexAreas[v0] - kappa_i / (vertexAreas[v0] * vertexAreas[v0]) * dadx.block<1, 3>(v0, v1  * 3);
                Eigen::Matrix<double, 1, 3> dcurvature0dxo0 = dcurvaturedxo0 / vertexAreas[v0] - kappa_i / (vertexAreas[v0] * vertexAreas[v0]) * dadx.block<1, 3>(v0, vo0 * 3);
                Eigen::Matrix<double, 1, 3> dcurvature0dxo1 = dcurvaturedxo1 / vertexAreas[v0] - kappa_i / (vertexAreas[v0] * vertexAreas[v0]) * dadx.block<1, 3>(v0, vo1 * 3);
                Eigen::Matrix<double, 1, 3> dcurvature1dx0 =  dcurvaturedx0  / vertexAreas[v1] - kappa_i / (vertexAreas[v1] * vertexAreas[v1]) * dadx.block<1, 3>(v1, v0  * 3);
                Eigen::Matrix<double, 1, 3> dcurvature1dx1 =  dcurvaturedx1  / vertexAreas[v1] - kappa_i / (vertexAreas[v1] * vertexAreas[v1]) * dadx.block<1, 3>(v1, v1  * 3);
                Eigen::Matrix<double, 1, 3> dcurvature1dxo0 = dcurvaturedxo0 / vertexAreas[v1] - kappa_i / (vertexAreas[v1] * vertexAreas[v1]) * dadx.block<1, 3>(v1, vo0 * 3);
                Eigen::Matrix<double, 1, 3> dcurvature1dxo1 = dcurvaturedxo1 / vertexAreas[v1] - kappa_i / (vertexAreas[v1] * vertexAreas[v1]) * dadx.block<1, 3>(v1, vo1 * 3);

                // assemble into dHdx, according to the region pairs
                for (int k = 0; k < rp_count[v0]; k++)
                {
                    if (Gamma_map[Gamma_map_inv[v0] + k].second[0] == region)
                    {
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, v0  * 3) += dcurvature0dx0;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, v1  * 3) += dcurvature0dx1;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, vo0 * 3) += dcurvature0dxo0;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, vo1 * 3) += dcurvature0dxo1;
                        H[Gamma_map_inv[v0] + k] += kappa_i / vertexAreas[v0];
                    }
                    if (Gamma_map[Gamma_map_inv[v0] + k].second[1] == region)
                    {
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, v0  * 3) -= dcurvature0dx0;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, v1  * 3) -= dcurvature0dx1;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, vo0 * 3) -= dcurvature0dxo0;
                        dHdx.block<1, 3>(Gamma_map_inv[v0] + k, vo1 * 3) -= dcurvature0dxo1;
                        H[Gamma_map_inv[v0] + k] -= kappa_i / vertexAreas[v0];
                    }
                }

                for (int k = 0; k < rp_count[v1]; k++)
                {
                    if (Gamma_map[Gamma_map_inv[v1] + k].second[0] == region)
                    {
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, v0  * 3) += dcurvature1dx0;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, v1  * 3) += dcurvature1dx1;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, vo0 * 3) += dcurvature1dxo0;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, vo1 * 3) += dcurvature1dxo1;
                        H[Gamma_map_inv[v1] + k] += kappa_i / vertexAreas[v1];
                    }
                    if (Gamma_map[Gamma_map_inv[v1] + k].second[1] == region)
                    {
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, v0  * 3) -= dcurvature1dx0;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, v1  * 3) -= dcurvature1dx1;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, vo0 * 3) -= dcurvature1dxo0;
                        dHdx.block<1, 3>(Gamma_map_inv[v1] + k, vo1 * 3) -= dcurvature1dxo1;
                        H[Gamma_map_inv[v1] + k] -= kappa_i / vertexAreas[v1];
                    }
                }
            }
        }
    }
    
    void indexGammaDofs(VS3D & vs, std::vector<std::pair<size_t, Vec2i> > & Gamma_map, std::vector<int> & Gamma_map_inv, std::vector<int> & rp_count)
    {
        size_t nv = vs.mesh().nv();
        
        for (size_t i = 0; i < nv; i++)
        {
            std::set<Vec2i, Vec2iComp> rps;
            for (size_t j = 0; j < vs.mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                LosTopos::Vec2i l = vs.mesh().get_triangle_label(vs.mesh().m_vertex_to_triangle_map[i][j]);
                Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
                rps.insert(rp);
            }
            rp_count.push_back(rps.size());
            Gamma_map_inv.push_back(Gamma_map.size());
            for (std::set<Vec2i, Vec2iComp>::iterator j = rps.begin(); j != rps.end(); j++)
                Gamma_map.push_back(std::pair<size_t, Vec2i>(i, *j));
        }
    }
}

void VS3D::step_implicit(double dt)
{
    size_t nv = mesh().nv();
    
    // Gamma dof indexing
    std::vector<std::pair<size_t, Vec2i> > Gamma_map;   // which vertex and region pair each Gamma dof corresponds to
    std::vector<int> Gamma_map_inv;     // the starting index in Gamma_map of the Gamma dofs on a given vertex
    std::vector<int> rp_count;          // the number of Gamma dofs on a given vertex
    indexGammaDofs(*this, Gamma_map, Gamma_map_inv, rp_count);
    
    
    
    // save the old state for convenience during the computation below
    VecXd oldpos = VecXd::Zero(nv * 3);
    for (int i = 0; i < nv; i++)
        oldpos.segment<3>(i * 3) = pos(i);
    VecXd newpos = oldpos;
    
    size_t ndof = Gamma_map.size();
    VecXd oldGamma = VecXd::Zero(ndof);
    for (size_t i = 0; i < ndof; i++)
        oldGamma[i] = (*m_Gamma)[Gamma_map[i].first].get(Gamma_map[i].second);
    VecXd newGamma = oldGamma;
    
    VecXd dGamma = VecXd::Zero(ndof);
    
    
    
    // iterate until convergence
    bool converged = false;
    int iter = 0;
    while (!converged && iter < 100)
    {
        std::cout << "Newton iteration " << iter << ": " << std::endl;
        iter++;
        
        MatXd dvdGamma = MatXd::Zero(nv * 3, ndof);
        compute_dvdGamma(dvdGamma, newpos, *this, Gamma_map, Gamma_map_inv, rp_count);
        
        VecXd H = VecXd::Zero(ndof);
        MatXd dHdx = MatXd::Zero(ndof, nv * 3);
        compute_dHdx(H, dHdx, newpos, *this, Gamma_map, Gamma_map_inv, rp_count);
        
        // implicit solve
        MatXd A = -dt * simOptions().sigma * (dHdx * dt * dvdGamma);
        for (size_t i = 0; i < ndof; i++)
            A(i, i) += 1;

        VecXd rhs = -dGamma + dt * simOptions().sigma * H;
        
        VecXd ddGamma = A.partialPivLu().solve(rhs);
        dGamma += ddGamma;
        
        if (ddGamma.norm() < 1e-8)
            converged = true;
        
        // Biot-Savart to update the mesh
        newGamma = oldGamma + dGamma;
        newpos = oldpos + dvdGamma * newGamma * dt;
    }
    
    
    
    // accept the result of the implicit solve
    for (size_t i = 0; i < ndof; i++)
        (*m_Gamma)[Gamma_map[i].first].set(Gamma_map[i].second, newGamma[i]);
    
    for (size_t i = 0; i < nv; i++)
        m_st->pm_newpositions[i] = vc(newpos.segment<3>(i * 3));
    
    
    
    // damping
//    for (size_t i = 0; i < ndof; i++)
//        (*m_Gamma)[Gamma_map[i].first].set(Gamma_map[i].second, (*m_Gamma)[Gamma_map[i].first].get(Gamma_map[i].second) * pow(simOptions().damping_coef, dt));

}

void VS3D::step_PBD_implicit(double dt)
{
    size_t nv = mesh().nv();

    // Gamma dof indexing
    std::vector<std::pair<size_t, Vec2i> > Gamma_map;   // which vertex and region pair each Gamma dof corresponds to
    std::vector<int> Gamma_map_inv;     // the starting index in Gamma_map of the Gamma dofs on a given vertex
    std::vector<int> rp_count;          // the number of Gamma dofs on a given vertex
    indexGammaDofs(*this, Gamma_map, Gamma_map_inv, rp_count);

    
    
    // save the old state for convenience during the computation below
    VecXd oldpos = VecXd::Zero(nv * 3);
    for (int i = 0; i < nv; i++)
        oldpos.segment<3>(i * 3) = pos(i);
    VecXd newpos = oldpos;
    
    size_t ndof = Gamma_map.size();
    VecXd oldGamma = VecXd::Zero(ndof);
    for (size_t i = 0; i < ndof; i++)
        oldGamma[i] = (*m_Gamma)[Gamma_map[i].first].get(Gamma_map[i].second);
    VecXd newGamma = oldGamma;
    
    VecXd dGamma = VecXd::Zero(ndof);

    
    
    // iterate until convergence
    bool converged = false;
    int iter = 0;
    while (!converged && iter < 100)
    {
        std::cout << "PBS-style iteration " << iter << ": " << std::endl;
        iter++;

        MatXd dvdGamma = MatXd::Zero(nv * 3, ndof);
        compute_dvdGamma(dvdGamma, newpos, *this, Gamma_map, Gamma_map_inv, rp_count);
        
        VecXd H = VecXd::Zero(ndof);
        MatXd dHdx = MatXd::Zero(ndof, nv * 3);
        compute_dHdx(H, dHdx, newpos, *this, Gamma_map, Gamma_map_inv, rp_count);

        // PBD-style solve
        MatXd gradC = -dt * simOptions().sigma * (dHdx * dt * dvdGamma);
        for (size_t i = 0; i < ndof; i++)
            gradC(i, i) += 1;
        
        VecXd C = dGamma - dt * simOptions().sigma * H;
        
        // solve one row at a time
        VecXd ddGamma = VecXd::Zero(ndof);
        for (size_t i = 0; i < ndof; i++)
        {
            VecXd gradCi = gradC.row(i);
            double Ci = C[i];
            double lambda = -Ci / gradCi.squaredNorm();
            ddGamma += lambda * gradCi;
        }
        dGamma += ddGamma;
        
        if (ddGamma.norm() < 1e-8)
            converged = true;
        
        // Biot-Savart to update the mesh
        newGamma = oldGamma + dGamma;
        newpos = oldpos + dvdGamma * newGamma * dt;
    }
    
    
    
    // accept the result of the implicit solve
    for (size_t i = 0; i < ndof; i++)
        (*m_Gamma)[Gamma_map[i].first].set(Gamma_map[i].second, newGamma[i]);
    
    for (size_t i = 0; i < nv; i++)
        m_st->pm_newpositions[i] = vc(newpos.segment<3>(i * 3));
    
    
    
    // damping
//    for (size_t i = 0; i < ndof; i++)
//        (*m_Gamma)[Gamma_map[i].first].set(Gamma_map[i].second, (*m_Gamma)[Gamma_map[i].first].get(Gamma_map[i].second) * pow(simOptions().damping_coef, dt));
    
}

