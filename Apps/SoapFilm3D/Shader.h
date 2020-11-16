//
//  Shader.h
//  RenderingTank
//
//  Created by Fang Da on 3/28/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//
//  Version History
//      1.0     4/7/2013    first version
//      1.0.1   4/7/2013    added an integer verion of setUniform, useful for setting texture samplers
//

#ifndef RenderingTank_Shader_h
#define RenderingTank_Shader_h
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
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

#include <string>
#include <vector>
#include <map>
#include "eigenheaders.h"

typedef Eigen::Matrix<float, 3, 1> Vec3;
typedef Eigen::Matrix<float, 4, 1> Vec4;
typedef Eigen::Matrix<float, 3, 3> Mat3;
typedef Eigen::Matrix<float, 4, 4> Mat4;

class Shader
{
public:
    Shader();
    
public:
    void activate();
    static void deactivate();
    
    static Shader * currentlyActiveShader() { return s_current; }
    
    static void setUniform(const std::string & name, int i);
    static void setUniform(const std::string & name, float f);
    static void setUniform(const std::string & name, const Vec3 & v);
    static void setUniform(const std::string & name, const Vec4 & v);
    static void setUniform(const std::string & name, const Vec3 * v, int n);
    static void setUniform(const std::string & name, const Vec4 * v, int n);
    static void setUniform(const std::string & name, const Mat3 & m, bool transpose = false);
    static void setUniform(const std::string & name, const Mat4 & m, bool transpose = false);
    
public:
    // call these before loading the shader from file (if necessary)
    void setVertexAttribName(const std::string & name, GLuint index);
    
    // load sources (.vert and .frag), compile and link
    bool loadFromFile(const std::string & filename);
    
    // query uniform locations
    GLint getUniformLocation(const std::string & name);
    
    // query attribute locations
    GLuint getVertexAttribLocation(const std::string & name);
    
protected:
    void printProgramInfoLog(GLuint prog, const std::string & name, const std::string & stage);
    void printShaderInfoLog(GLuint shader, const std::string & name, const std::string & stage);
    
protected:
    std::map<std::string, GLuint> m_vertex_attributes;
    std::map<std::string, GLint> m_uniform_locations;

    GLuint m_program;
    GLuint m_vs;
    GLuint m_fs;

    static Shader * s_current;
    
};


#endif
