// ---------------------------------------------------------
//
//  meshrenderer.cpp
//  Tyson Brochu 2011
//  Christopher Batty, Fang Da 2014
//
//  OpenGL rendering for a triangle mesh.
//
// ---------------------------------------------------------

#include <meshrenderer.h>

#ifndef NO_GUI

#include <dynamicsurface.h>
#include <gluvi.h>
#include "trianglequality.h"

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// OpenGL render a mesh surface
///
// ---------------------------------------------------------

namespace LosTopos {

void MeshRenderer::render( const DynamicSurface& surface )
{

   glDisable(GL_LIGHTING);
   glDepthFunc(GL_LEQUAL);

   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(1.0f, 1.0f);      //  allow the wireframe to show through
   
   glLineWidth(5);
   glBegin(GL_LINES);

   for(size_t e = 0; e < surface.m_mesh.m_edges.size(); e++)
   {

      /*if(surface.edge_is_feature(e)) {
         const Vec2st& edge = surface.m_mesh.m_edges[e];
         const Vec3d& vtx0 = surface.get_position(edge[0]);
         const Vec3d& vtx1 = surface.get_position(edge[1]);
         glColor3d(0,1,0);
         glVertex3d(vtx0[0], vtx0[1], vtx0[2]);
         glVertex3d(vtx1[0], vtx1[1], vtx1[2]);
      }*/
   }

   glEnd(); 

   /*
    //
    // edges
    //
    
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    
    if ( render_edges )
    {
        glLineWidth(2);
        glBegin(GL_LINES);
        
        for(size_t e = 0; e < surface.m_mesh.m_edges.size(); e++)
        {
            if ( surface.m_mesh.m_is_boundary_edge[e] )
            {
                glColor3d(1,0,0);
            }
            else
            {
                glColor3d(0,0,0);
            }
            
            const Vec2st& edge = surface.m_mesh.m_edges[e];
            const Vec3d& vtx0 = surface.get_position(edge[0]);
            const Vec3d& vtx1 = surface.get_position(edge[1]);
            glVertex3d(vtx0[0], vtx0[1], vtx0[2]);
            glVertex3d(vtx1[0], vtx1[1], vtx1[2]);
        }
        
        glEnd(); 
    }
    
    //
    // vertices
    //
    
  

    if ( render_vertex_rank )
    {
        glPointSize(10);
        glBegin(GL_POINTS);
        
        for ( size_t v = 0; v < surface.get_num_vertices(); ++v )
        {
            if ( surface.m_mesh.m_vertex_to_triangle_map[v].empty() )
            {
                continue;
            }
            
            if ( surface.vertex_is_solid(v) )
            {
                glColor3f( 1.0f, 0.0f, 0.0f );
            }
            else
            {
                glColor3f( 0.0f, 1.0f, 0.0f );
            }
            
            glVertex3dv( surface.get_position(v).v );      
            
        }
        glEnd();
    }   
    
    //
    // triangles
    //
    
    if ( render_fill_triangles )
    {
        if ( two_sided )
        {
            glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );
        }
        else
        {
            glEnable(GL_CULL_FACE);
        }
        
        glEnable(GL_LIGHTING);
        glShadeModel(GL_SMOOTH);
        Gluvi::set_generic_lights();
        Gluvi::set_generic_material(1.0f, 1.0f, 1.0f, GL_FRONT);   // exterior surface colour
        Gluvi::set_generic_material(1.0f, 1.0f, 1.0f, GL_BACK);
        
        //if ( !smooth_shading )
        //{
        //glDisable(GL_LIGHTING);
        //glColor3d(1,1,1);
        //}
        
        if ( render_edges || render_nonmanifold_curve)
        {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0f, 1.0f);      //  allow the wireframe to show through
        }
        
        glEnable(GL_DEPTH_TEST);
//        glDepthMask(GL_TRUE);
        glDepthMask(GL_FALSE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        
        glBegin(GL_TRIANGLES);
        if(smooth_shading) {
           for(size_t i = 0; i < surface.m_mesh.num_triangles(); i++)
           {
               const Vec3st& tri = surface.m_mesh.get_triangle(i);
            
               const Vec3d& v0 = surface.get_position(tri[0]);
               const Vec3d& v1 = surface.get_position(tri[1]);
               const Vec3d& v2 = surface.get_position(tri[2]);
            
               glNormal3dv( surface.get_vertex_normal(tri[0]).v );
               glVertex3d(v0[0], v0[1], v0[2]);
            
               glNormal3dv( surface.get_vertex_normal(tri[1]).v );
               glVertex3d(v1[0], v1[1], v1[2]);
            
               glNormal3dv( surface.get_vertex_normal(tri[2]).v );
               glVertex3d(v2[0], v2[1], v2[2]);
            
           }
        }
        else {
           for(size_t i = 0; i < surface.m_mesh.num_triangles(); i++)
           {
              const Vec3st& tri = surface.m_mesh.get_triangle(i);

              const Vec3d& v0 = surface.get_position(tri[0]);
              const Vec3d& v1 = surface.get_position(tri[1]);
              const Vec3d& v2 = surface.get_position(tri[2]);
              
              const Vec3d& normal = surface.get_triangle_normal(i);
              glNormal3dv( normal.v );
              glVertex3d(v0[0], v0[1], v0[2]);

              glNormal3dv( normal.v );
              glVertex3d(v1[0], v1[1], v1[2]);

              glNormal3dv( normal.v );
              glVertex3d(v2[0], v2[1], v2[2]);

           }
        }
        
        glEnd();
        
       
        if ( render_edges || render_nonmanifold_curve)
        {
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
        
        glDisable(GL_LIGHTING);
        
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }*/
    

    
}


// ---------------------------------------------------------
///
/// OpenGL render a mesh surface
///
// ---------------------------------------------------------

void MeshRenderer::render(const std::vector<Vec3d>& xs,
                          const std::vector<Vec3d>& normals,
                          const std::vector<Vec3st>& triangles,
                          const std::vector<Vec3d>& tri_normals,
                          const std::vector<Vec2st>& edges,
                          const std::vector<int>& ranks,
                          const std::vector<bool>& edge_manifold,
                          const std::vector<Vec2i>& labels)
{
    
    //
    // edges
    //
    
   std::vector<int> vert_incidences(xs.size(),0);
   for(unsigned int i = 0; i < triangles.size(); ++i) {
      Vec3st tri = triangles[i];
      vert_incidences[tri[0]]++;
      vert_incidences[tri[1]]++;
      vert_incidences[tri[2]]++;
   }

    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    
    if ( render_edges )
    {
        glLineWidth(2);
        glColor3f( 0.0f, 0.0f, 0.0f );
        glBegin(GL_LINES);
        for(size_t e = 0; e < edges.size(); e++)
        {
            const Vec2st& edge = edges[e];
            const Vec3d& vtx0 = xs[edge[0]];
            const Vec3d& vtx1 = xs[edge[1]];

            if(render_nonmanifold_curve && !edge_manifold[e])
               glColor3d(1,0,0);
            else
               glColor3f(0,0,0);
            glVertex3dv( vtx0.v );
            glVertex3dv( vtx1.v );
        }
        glEnd(); 
    }

    
    if ( render_nonmanifold_curve )
    {
       glLineWidth(2);
       glColor3f( 0.0f, 0.0f, 0.0f );
       glBegin(GL_LINES);
       for(size_t e = 0; e < edges.size(); e++)
       {
          const Vec2st& edge = edges[e];
          const Vec3d& vtx0 = xs[edge[0]];
          const Vec3d& vtx1 = xs[edge[1]];
           if(!edge_manifold[e]) {
              glColor3d(1,0,0);
              glVertex3dv( vtx0.v );
              glVertex3dv( vtx1.v );
          }

         
       }
       glEnd(); 
    }

    //glLineWidth(3);
    //glColor3f(0,0,1);
    //glBegin(GL_LINES);
    //for(size_t i = 0; i < triangles.size(); i++)
    //{
    //   const Vec3st& tri = triangles[i];

    //   const Vec3d& v0 = xs[tri[0]];
    //   const Vec3d& v1 = xs[tri[1]];
    //   const Vec3d& v2 = xs[tri[2]];
    //   Vec3d avg = 0.3333*(v0+v1+v2);
    //   Vec3d normal = cross(v1-v0, v2-v0);
    //   normalize(normal);

    //   glVertex3dv(avg.v);
    //   glVertex3dv((avg - 0.01*normal).v);

    //}
    //glEnd();
    
    //
    // vertices
    //

    //Vertex valences, sort of.
    //glPointSize(10);
    //glBegin(GL_POINTS);
    //for ( size_t v = 0; v < xs.size(); ++v )
    //{
    //   if(vert_incidences[v] < 5) {
    //      glColor3d(1,0,0);
    //      glVertex3dv( xs[v].v );     
    //   }
    //   else if(vert_incidences[v] > 7){
    //      glColor3d(0,0,1);
    //      glVertex3dv( xs[v].v );         
    //   }
    //}
    //glEnd();

    if ( render_vertex_rank )
    {
        glPointSize(7);
        glBegin(GL_POINTS);
        for ( size_t v = 0; v < xs.size(); ++v )
        {
           if(ranks[v] == 0)
              glColor3d(0,0,0);
           else if(ranks[v] == 1)
              glColor3d(1,0,0);
           else if(ranks[v] == 2)
              glColor3d(0,1,0);
           else if(ranks[v] == 3)
              glColor3d(0,0,1);
           
           if(ranks[v] > 1)
               glVertex3dv( xs[v].v );               

        }
        glEnd();
    }   
    
    
    //
    // face labels
    //
    
    const static float colors[][3] = 
    {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 0, 1, 1 },
        { 1, 0, 1 },
        { 1, 1, 0 },
        { 0, 0.5, 1 },
        { 1, 0, 0.5 },
        { 0.5, 1, 0 },
        { 0, 1, 0.5 },
        { 0.5, 0, 1 },
        { 1, 0.5, 0 }
    };
    
    if (render_face_labels)
    {
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        for(size_t i = 0; i < triangles.size(); i++)
        {
            const Vec3st & tri = triangles[i];
            const Vec2i & label = labels[i];
            Vec3d x0 = xs[tri[0]];
            Vec3d x1 = xs[tri[1]];
            Vec3d x2 = xs[tri[2]];
            Vec3d n = cross(x1 - x0, x2 - x0);
            n /= mag(n);
            Vec3d c = (x0 + x1 + x2) / 3;
            double e = (mag(x1 - x0) + mag(x2 - x1) + mag(x0 - x2)) / 3;
            glColor3fv(colors[label[0]]);
            glVertex3dv(c.v);
            glVertex3dv((c - n * e * 0.1).v);
            glColor3fv(colors[label[1]]);
            glVertex3dv(c.v);
            glVertex3dv((c + n * e * 0.1).v);
        } 
        glEnd();
    }
    
    //
    // triangles
    //
    
    if ( render_fill_triangles )
    {
        if ( two_sided )
        {
            glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );
        }
        else
        {
            glEnable(GL_CULL_FACE);
        }
        
        glEnable(GL_LIGHTING);
        glShadeModel(GL_SMOOTH);
        Gluvi::set_generic_lights();
        Gluvi::set_generic_material(1.0f, 1.0f, 1.0f, GL_FRONT);   // exterior surface colour
        Gluvi::set_generic_material(1.0f, 1.0f, 1.0f, GL_BACK);
        
        /*if ( !smooth_shading )
        {
            glDisable(GL_LIGHTING);
            glColor3d(1,1,1);
        }*/
        
        if ( render_edges || render_nonmanifold_curve)
        {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0f, 1.0f);      //  allow the wireframe to show through
        }
        
         if(render_transparent_triangles) {
            glDisable(GL_DEPTH_TEST);
            //glDepthMask(GL_TRUE);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable(GL_BLEND);

         }
         else {
            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);
         }
        
        glBegin(GL_TRIANGLES);
        if(smooth_shading) {
           for(size_t i = 0; i < triangles.size(); i++)
           {
               const Vec3st& tri = triangles[i];
               glNormal3dv( normals[tri[0]].v );
               glVertex3dv( xs[tri[0]].v );
               glNormal3dv( normals[tri[1]].v );
               glVertex3dv( xs[tri[1]].v );
               glNormal3dv( normals[tri[2]].v );
               glVertex3dv( xs[tri[2]].v );
           }      
        }
        else {
           for(size_t i = 0; i < triangles.size(); i++)
           {
              const Vec3st& tri = triangles[i];
              const Vec3d& normal = tri_normals[i];
              glNormal3dv( normal.v );
              glVertex3dv( xs[tri[0]].v );
              glNormal3dv( normal.v );
              glVertex3dv( xs[tri[1]].v );
              glNormal3dv( normal.v );
              glVertex3dv( xs[tri[2]].v );
           }
        }
        glEnd();
        
        if ( render_edges || render_nonmanifold_curve)
        {
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
        
        glDisable(GL_LIGHTING);
        
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }
    
}

}

#endif // ndef NO_GUI

