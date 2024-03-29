//
//  PRRenderer.h
//  MultiTracker
//
//  Created by Fang Da on 1/15/15.
//
//

#ifndef __MultiTracker__PRRenderer__
#define __MultiTracker__PRRenderer__

#include <fstream>
#include <iostream>
#include <sstream>

#include "VS3D.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#ifdef _MSC_VER
#define NOMINMAX
#include <Windows.h>
#include <glew.h>
#else
#include <GL/glew.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#endif
#include "Shader.h"

class PRRenderer {
 public:
  PRRenderer(VS3D* vs);

 public:
  void render();

 protected:
  LosTopos::NonDestructiveTriMesh& mesh() { return m_vs->mesh(); }

  bool create_cube_map(GLuint& tex);
  bool load_cube_map_side(GLuint tex, GLenum side, const unsigned char* content,
                          unsigned int content_size);

 protected:
  VS3D* m_vs;
  Shader m_shader_bubble;
  Shader m_shader_env;
  GLuint m_tex_env;
};

#endif /* defined(__MultiTracker__PRRenderer__) */
