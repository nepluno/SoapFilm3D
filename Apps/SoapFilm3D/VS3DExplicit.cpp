//
//  VS3DExplicit.cpp
//  MultiTracker
//
//  Created by Fang Da on 15/1/27.
//
//

#include "VS3D.h"
#include "SimOptions.h"

#ifndef WIN32
#include "fmmtl/fmmtl/KernelMatrix.hpp"
#include "fmmtl/kernel/BiotSpherical.hpp"
#include "fmmtl/kernel/RMSpherical.hpp"
#include "fmmtl/fmmtl/Direct.hpp"
#include "fmmtl/fmmtl/util/Clock.hpp"
#endif

namespace
{
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

    VecXd BiotSavart_naive(VS3D & vs, const VecXd & dx)
    {
        VecXd vel = VecXd::Zero(vs.mesh().nv() * 3);
        for (size_t i = 0; i < vs.mesh().nv(); i++)
        {
            Vec3d v(0, 0, 0);
            Vec3d x = vs.pos(i);
            
            for (size_t j = 0; j < vs.mesh().nt(); j++)
            {
                LosTopos::Vec3st t = vs.mesh().get_triangle(j);
                if (vs.surfTrack()->vertex_is_any_solid(t[0]) && vs.surfTrack()->vertex_is_any_solid(t[1]) && vs.surfTrack()->vertex_is_any_solid(t[2]))
                    continue;   // all-solid faces don't contribute vorticity.
                
                LosTopos::Vec2i l = vs.mesh().get_triangle_label(j);
                Vec3d x0 = vs.pos(t[0]) + dx.segment<3>(t[0] * 3);
                Vec3d x1 = vs.pos(t[1]) + dx.segment<3>(t[1] * 3);
                Vec3d x2 = vs.pos(t[2]) + dx.segment<3>(t[2] * 3);
                
                Vec3d xp = (x0 + x1 + x2) / 3;
                
                Vec3d e01 = x1 - x0;
                Vec3d e12 = x2 - x1;
                Vec3d e20 = x0 - x2;
                
                Vec3d gamma =  -(e01 * vs.Gamma(t[2]).get(l) +
                                 e12 * vs.Gamma(t[0]).get(l) +
                                 e20 * vs.Gamma(t[1]).get(l));
                
                Vec3d dx = x - xp;
//                double dxn = dx.norm();
                double dxn = sqrt(dx.squaredNorm() + vs.delta() * vs.delta());
                
                v += gamma.cross(dx) / (dxn * dxn * dxn);
//                v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 - exp(-dxn / m_delta));
            }
            
            // open boundary extra face contributions
            if (vs.m_obefv.size() == vs.m_obefe.size() && vs.m_obefv.size() == vs.m_obefc.size())
            {
                for (size_t j = 0; j < vs.m_obefv.size(); j++)
                {
                    Vec3d gamma = vs.m_obefe[j] * vs.m_obefv[j];
                    Vec3d xp = vs.m_obefc[j];
                    
                    Vec3d dx = x - xp;
//                    double dxn = dx.norm();
                    double dxn = sqrt(dx.squaredNorm() + vs.delta() * vs.delta());
                    
                    v += gamma.cross(dx) / (dxn * dxn * dxn);
//                    v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 - exp(-dxn / m_delta));
                }
            }
        
            v /= (4 * M_PI);
            vel.segment<3>(i * 3) = v;
        }
        
        return vel;
    }
#ifndef WIN32
    VecXd BiotSavart_fmmtl(VS3D & vs, const VecXd & dx)
    {
        // code adapted from FMMTL example test "error_biot.cpp"
        
        // Init the FMM Kernel and options
        FMMOptions opts = get_options(0, NULL);
//        typedef BiotSpherical kernel_type;
        typedef RMSpherical kernel_type;
        
        // Init kernel
        kernel_type K(vs.delta());
        
        typedef kernel_type::point_type point_type;
        typedef kernel_type::source_type source_type;
        typedef kernel_type::target_type target_type;
        typedef kernel_type::charge_type charge_type;
        typedef kernel_type::result_type result_type;
        
        // Init points and charges
        std::vector<source_type> sources;
        std::vector<target_type> targets;
        std::vector<charge_type> charges;
        
        for (size_t i = 0; i < vs.mesh().nv(); i++)
        {
            Vec3d x = vs.pos(i);
            targets.push_back(Vec<3, double>(x[0], x[1], x[2]));
        }
        
        for (size_t j = 0; j < vs.mesh().nt(); j++)
        {
            LosTopos::Vec3st t = vs.mesh().get_triangle(j);
            if (vs.surfTrack()->vertex_is_any_solid(t[0]) && vs.surfTrack()->vertex_is_any_solid(t[1]) && vs.surfTrack()->vertex_is_any_solid(t[2]))
                continue;   // all-solid faces don't contribute vorticity.
            
            LosTopos::Vec2i l = vs.mesh().get_triangle_label(j);
            Vec3d x0 = vs.pos(t[0]) + dx.segment<3>(t[0] * 3);
            Vec3d x1 = vs.pos(t[1]) + dx.segment<3>(t[1] * 3);
            Vec3d x2 = vs.pos(t[2]) + dx.segment<3>(t[2] * 3);
            
            Vec3d xp = (x0 + x1 + x2) / 3;
            
            Vec3d e01 = x1 - x0;
            Vec3d e12 = x2 - x1;
            Vec3d e20 = x0 - x2;
            
            Vec3d gamma =  -(e01 * vs.Gamma(t[2]).get(l) +
                             e12 * vs.Gamma(t[0]).get(l) +
                             e20 * vs.Gamma(t[1]).get(l));
            
//            Vec3d dx = x - xp;
////            double dxn = dx.norm();
//            double dxn = sqrt(dx.squaredNorm() + vs.delta() * vs.delta());
//            
//            v += gamma.cross(dx) / (dxn * dxn * dxn);
////            v += gamma.cross(dx) / (dxn * dxn * dxn) * (1 - exp(-dxn / m_delta));
            
            sources.push_back(Vec<3, double>(xp[0], xp[1], xp[2]));
            charges.push_back(Vec<3, double>(gamma[0], gamma[1], gamma[2]));
        }
        
        
        
        // Build the FMM
        fmmtl::kernel_matrix<kernel_type> A = K(targets, sources);
        A.set_options(opts);
        
        // Execute the FMM
        Clock t1;
        std::vector<result_type> result = A * charges;
        double time1 = t1.seconds();
        
        VecXd vel = VecXd::Zero(vs.mesh().nv() * 3);
        for (size_t i = 0; i < vs.mesh().nv(); i++)
        {
            vel[i * 3 + 0] = result[i][0];
            vel[i * 3 + 1] = result[i][1];
            vel[i * 3 + 2] = result[i][2];
        }
        
        vel /= (4 * M_PI);
        
        return vel;
    }
#endif
    VecXd BiotSavart(VS3D & vs, const VecXd & dx)
    {
#ifndef WIN32
        if (Options::boolValue("fmmtl"))
            return BiotSavart_fmmtl(vs, dx);
        else
#endif
            return BiotSavart_naive(vs, dx);
    }

//#define FANGS_VERSION
//#define FANGS_PATCHED

void VS3D::step_explicit(double dt)
{
    std::cout << "Explicit time stepping" << std::endl;
    m_dbg_t1.clear();
    m_dbg_t2.clear();
    m_dbg_e1.clear();
    m_dbg_v1.clear();
    
    m_dbg_t1.resize(mesh().nt());
    m_dbg_t2.resize(mesh().nt());
    m_dbg_e1.resize(mesh().ne(), std::vector<double>(m_nregion, 0));
    m_dbg_v1.resize(mesh().nv(), std::vector<double>(m_nregion, 0));
    
    // integrate surface tension force
    
    std::vector<std::vector<double> > curvature(mesh().ne(), std::vector<double>(m_nregion, 0));  // edge-aligned curvature (signed scalar)
    std::vector<std::vector<double> > mean_curvatures(mesh().nv(), std::vector<double>(m_nregion, 0));   // vertex-aligned mean curvature (signed scalar)
    std::vector< double > avg_vertex_areas(mesh().nv(), 0.);
  
  for (size_t i = 0; i < mesh().nv(); i++)
  {
    std::set<int> incident_region;
    double vertex_area = 0.0;
    for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
    {
      const LosTopos::Vec3st & t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
      const LosTopos::Vec2i & l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
      
      for(int r = 0; r < 2; ++r) {
        incident_region.insert(l[r]);
      }
      
      Vec3d x0 = pos(t[0]);
      Vec3d x1 = pos(t[1]);
      Vec3d x2 = pos(t[2]);
      double area = (x1 - x0).cross(x2 - x0).norm() / 2;
      
      vertex_area += area / 3 * 2.0;
    }
    
    if(incident_region.size() > 0) {
      avg_vertex_areas[i] = vertex_area / (double) incident_region.size();
    }
  }
  
    for (size_t region = 0; region < m_nregion; region++)
    {
        for (size_t i = 0; i < mesh().ne(); i++)
        {
            std::vector<size_t> incident_faces; // faces incident to edge i that have the label of interest (assume there are only two of them for now; this can be false only when complex collision prevents immediate T1 resolution, which is not expected to happen for bubble complexes.)
            std::set<int> incident_regions;
            for (size_t j = 0; j < mesh().m_edge_to_triangle_map[i].size(); j++)
            {
                const LosTopos::Vec2i & l = mesh().get_triangle_label(mesh().m_edge_to_triangle_map[i][j]);
                if (l[0] == region || l[1] == region)
                    incident_faces.push_back(j);
                incident_regions.insert(l[0]);
                incident_regions.insert(l[1]);
            }
            if (incident_faces.size() == 0)
                continue;
//            assert(incident_faces.size() == 2);
            if (incident_faces.size() != 2)
            {
//                std::cout << "Warning: incident_faces.size() != 2" << std::endl;
                curvature[i][region] = 0;
                continue;
            }
            bool nonmanifold = (incident_regions.size() > 2);
            
//            assert(mesh().m_edge_to_triangle_map[i].size() == 2);
            int v0 = mesh().m_edges[i][0];
            int v1 = mesh().m_edges[i][1];
            Vec3d x0 = pos(v0);
            Vec3d x1 = pos(v1);
            Vec3d et = (x1 - x0);
            
            int ti0 = mesh().m_edge_to_triangle_map[i][incident_faces[0]];
            int ti1 = mesh().m_edge_to_triangle_map[i][incident_faces[1]];
            LosTopos::Vec3st t0 = mesh().get_triangle(ti0);
            LosTopos::Vec3st t1 = mesh().get_triangle(ti1);
            
            if (mesh().get_triangle_label(ti0)[mesh().oriented(v0, v1, t0) ? 1 : 0] == region)
                std::swap(ti0, ti1),
                std::swap(t0, t1);  // the region of interest should be to the CCW direction of t0 when looking along the direciton of edge i
            
            assert(mesh().get_triangle_label(ti0)[mesh().oriented(v0, v1, t0) ? 0 : 1] == region);
            assert(mesh().get_triangle_label(ti1)[mesh().oriented(v0, v1, t1) ? 1 : 0] == region);
            
            Vec3d n0 = (pos(t0[1]) - pos(t0[0])).cross(pos(t0[2]) - pos(t0[0]));
            if (mesh().get_triangle_label(ti0)[1] == region)
                n0 = -n0;   // n0 should point away from the region of interest
            
            Vec3d n1 = (pos(t1[1]) - pos(t1[0])).cross(pos(t1[2]) - pos(t1[0]));
            if (mesh().get_triangle_label(ti1)[1] == region)
                n1 = -n1;   // n1 should point away from the region of interest
            
            double edgeArea = (n0.norm() + n1.norm()) / 2 / 3;
            n0.normalize();
            n1.normalize();
            
//            double curvature_i = 0;
//            if (nonmanifold)
//            {
//                curvature_i = n0.cross(n1).dot(et) / (et.norm() * m_delta); // integral curvature along the curve orghotonal to the edge in plane
//            } else
//            {
//                curvature_i = n0.cross(n1).dot(et) / edgeArea;
////                curvature_i = angleAroundAxis(n0, n1, et) * et.norm() / edgeArea;
//            }

//            double curvature_i = n0.cross(n1).dot(et);
//            double curvature_i = (n0 - n1).dot((n0 + n1).normalized().cross(et));
            double curvature_i = angleAroundAxis(n0, n1, et) * et.norm();
            
            
            curvature[i][region] = curvature_i;
//            if (region == 1)
//                m_dbg_e1[i] = curvature_i;
            m_dbg_e1[i][region] = curvature_i;
        }
        
//        for (size_t i = 0; i < mesh().nv(); i++)
//        {
//            double mean_curvature = 0;
//            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
//                mean_curvature += curvature[mesh().m_vertex_to_edge_map[i][j]][region];
//            mean_curvature /= mesh().m_vertex_to_edge_map[i].size();
//
//            (*m_Gamma)[i][region] += simOptions().sigma * mean_curvature * dt;
//            m_dbg_v1[i][region] = mean_curvature;
//        }
        
        for (size_t i = 0; i < mesh().nv(); i++)
        {
            double mean_curvature = 0;
            
            Mat3d second_fundamental_form = Mat3d::Zero();
            int counter = 0;

#ifdef FANGS_VERSION
            double vertex_area = 0;
            double triple_junction_length_sum = 0;
#else
            double vertex_area = avg_vertex_areas[i];
#endif
            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                size_t e = mesh().m_vertex_to_edge_map[i][j];
                bool incident_to_region = false;
                for (size_t k = 0; k < mesh().m_edge_to_triangle_map[e].size(); k++)
                {
                    const LosTopos::Vec2i & l = mesh().get_triangle_label(mesh().m_edge_to_triangle_map[e][k]);
                    if (l[0] == region || l[1] == region)
                    {
                        incident_to_region = true;
                        break;
                    }
                }
                
                if (incident_to_region)
                {
                    Vec3d et = (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).normalized();
                    second_fundamental_form += et * et.transpose() * curvature[e][region];
                    mean_curvature += curvature[e][region];
                    counter++;
#ifdef FANGS_VERSION
                    if (mesh().m_edge_to_triangle_map[e].size() > 2)
                        triple_junction_length_sum += (pos(mesh().m_edges[e][1]) - pos(mesh().m_edges[e][0])).norm();
#endif
                }
            }
#ifdef FANGS_VERSION
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                const LosTopos::Vec3st & t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
                const LosTopos::Vec2i & l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
                if (l[0] == region || l[1] == region)
                {
                    Vec3d x0 = pos(t[0]);
                    Vec3d x1 = pos(t[1]);
                    Vec3d x2 = pos(t[2]);
                    double area = (x1 - x0).cross(x2 - x0).norm() / 2;
                    vertex_area += area / 3;
                }
            }
#endif
//            if (counter == 0)
//                mean_curvature = 0;
//            else
////                mean_curvature = second_fundamental_form.trace() / counter / 3;
//                mean_curvature /= counter;
            
#ifdef FANGS_VERSION
#ifdef FANGS_PATCHED
            if (triple_junction_length_sum > 0) // this means the vertex is a triple junction vertex; the vertex domain should then be smaller.
                vertex_area = triple_junction_length_sum * m_delta * 1.0;
#endif
#endif
            
            if (vertex_area == 0)
                mean_curvature = 0;
            else
                mean_curvature = mean_curvature / (vertex_area * 2);
            
//            (*m_Gamma)[i][region] += simOptions().sigma * mean_curvature * dt;
            m_dbg_v1[i][region] = mean_curvature;
            
            mean_curvatures[i][region] = mean_curvature;
            
        }
    }
    
    
    // Integrate vertex mean curvatures into vertex Gammas (skipping constrained vertices)
    std::vector<bool> constrained(mesh().nv(), false);
    for (size_t i = 0; i < m_constrained_vertices.size(); i++)
        constrained[m_constrained_vertices[i]] = true;
    
    for (size_t i = 0; i < mesh().nv(); i++)
    {
        if (constrained[i])
            continue;
        
        std::set<Vec2i, Vec2iComp> incident_region_pairs;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            LosTopos::Vec2i l = mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[i][j]);
            incident_region_pairs.insert(l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
        }
        
        for (std::set<Vec2i, Vec2iComp>::iterator j = incident_region_pairs.begin(); j != incident_region_pairs.end(); j++)
        {
            Vec2i rp = *j;
            double mean_curvature = mean_curvatures[i][rp[0]] - mean_curvatures[i][rp[1]];
            (*m_Gamma)[i].set(rp, (*m_Gamma)[i].get(rp) + simOptions().sigma * mean_curvature * dt);
        }
    }
    
    
    // Biot-Savart integral based on vortex sheet strength gamma: Stock 2006, Eq. 2.26
    VecXd v;
    if (Options::boolValue("RK4-velocity-integration"))
    {
        // RK4 integration
        VecXd v1 = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
        VecXd v2 = BiotSavart(*this, v1 * dt * 0.5);
        VecXd v3 = BiotSavart(*this, v2 * dt * 0.5);
        VecXd v4 = BiotSavart(*this, v3 * dt);
        v = (v1 + 2 * v2 + 2 * v3 + v4) / 6;
    } else
    {
        // explicit Euler
        v = BiotSavart(*this, VecXd::Zero(mesh().nv() * 3));
    }
  
    for (size_t i = 0; i < mesh().nv(); i++)
        m_st->pm_newpositions[i] = m_st->pm_positions[i] + vc(v.segment<3>(i * 3)) * dt;


    // damping
    for (size_t i = 0; i < mesh().nv(); i++)
        (*m_Gamma)[i].values *= pow(simOptions().damping_coef, dt);
    
    std::cout << "Explicit time stepping finished" << std::endl;
}
