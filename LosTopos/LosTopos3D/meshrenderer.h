// ---------------------------------------------------------
//
//  meshrenderer.h
//  Tyson Brochu 2011
//  Christopher Batty, Fang Da 2014
//
//  OpenGL rendering for a triangle mesh.
//
// ---------------------------------------------------------

#ifndef LOSTOPOS_MESHRENDERER_H
#define LOSTOPOS_MESHRENDERER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace LosTopos {

class DynamicSurface;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Mesh rendering object.  Contains current rendering options and functions for doing OpenGL render of a mesh.
///
// ---------------------------------------------------------

class MeshRenderer
{
    
public:
    
    /// Constructor
    ///
    MeshRenderer() :
    render_edges( true ),
    render_fill_triangles( true ),
    render_vertex_rank( false ),
    smooth_shading( false ),
    two_sided( true ),
    render_nonmanifold_curve(true),
    render_face_labels(false)
    {}
    
    /// Whether to show mesh edges (wireframe)
    ///
    bool render_edges;
    
    /// Whether to render filled triangles
    ///
    bool render_fill_triangles;
    
    /// Whether to render transparent triangles
    ///
    bool render_transparent_triangles;

    /// Whether to render the primary-space rank for each vertex
    ///
    bool render_vertex_rank;
    
    /// Whether to use smooth or flat shading
    ///    
    bool smooth_shading;

    /// Render both sides of the triangles
    ///    
    bool two_sided;

    /// Render the non-manifold curve highlighted red
    ///    
    bool render_nonmanifold_curve;
    
    /// Render the face labels
    ///    
    bool render_face_labels;
    
    /// Display the surface in OpenGL using the current options settings
    ///
    void render( const DynamicSurface& surface );

    /// Display the specified geometry in OpenGL using the current options settings
    ///
    void render(const std::vector<Vec3d>& xs,
                const std::vector<Vec3d>& normals,
                const std::vector<Vec3st>& triangles,
                const std::vector<Vec3d>& tri_normals,
                const std::vector<Vec2st>& edges,
                const std::vector<int>& ranks,
                const std::vector<bool>& edge_manifold,
                const std::vector<Vec2i>& labels);
    
};

}

#endif
