
//  Shader.cpp
//  RenderingTank
//
//  Created by Fang Da on 3/28/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "Shader.h"

#include <assert.h>

#include <fstream>
#include <iostream>

Shader* Shader::s_current = NULL;

Shader::Shader() : m_program(0), m_vs(0), m_fs(0) {}

void Shader::activate() {
  glUseProgram(m_program);
  s_current = this;
  assert(glGetError() == GL_NO_ERROR);
}

void Shader::deactivate() {
  glUseProgram(0);
  s_current = NULL;
}

void Shader::setUniform(const std::string& name, int i) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniform1i(current->getUniformLocation(name), i);
}

void Shader::setUniform(const std::string& name, float f) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniform1f(current->getUniformLocation(name), f);
}

void Shader::setUniform(const std::string& name, const Vec3& v) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniform3f(current->getUniformLocation(name), v.x(), v.y(), v.z());
}

void Shader::setUniform(const std::string& name, const Vec4& v) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniform4f(current->getUniformLocation(name), v.x(), v.y(), v.z(), v.w());
}

void Shader::setUniform(const std::string& name, const Vec3* v, int n) {
  Shader* current = currentlyActiveShader();
  assert(current);
  float* buffer = new float[n * 3];
  for (int i = 0; i < n; i++)
    memcpy(buffer + i * 3, (v + i)->data(), sizeof(float) * 3);
  glUniform3fv(current->getUniformLocation(name), n, buffer);
  delete[] buffer;
}

void Shader::setUniform(const std::string& name, const Vec4* v, int n) {
  Shader* current = currentlyActiveShader();
  assert(current);
  float* buffer = new float[n * 4];
  for (int i = 0; i < n; i++)
    memcpy(buffer + i * 4, (v + i)->data(), sizeof(float) * 4);
  glUniform4fv(current->getUniformLocation(name), n, buffer);
  delete[] buffer;
}

void Shader::setUniform(const std::string& name, const Mat3& m,
                        bool transpose) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniformMatrix3fv(current->getUniformLocation(name), 1, GL_FALSE, m.data());
}

void Shader::setUniform(const std::string& name, const Mat4& m,
                        bool transpose) {
  Shader* current = currentlyActiveShader();
  assert(current);
  glUniformMatrix4fv(current->getUniformLocation(name), 1, GL_FALSE, m.data());
}

bool Shader::loadFromFile(const std::string& name) {
  std::string vs_str;
  std::string fs_str;
  const char* vs_source = NULL;
  const char* fs_source = NULL;

  std::ifstream vs_file((name + ".vert").c_str());
  if (vs_file.is_open()) {
    vs_str = std::string((std::istreambuf_iterator<char>(vs_file)),
                         std::istreambuf_iterator<char>());
    vs_source = vs_str.c_str();
  }

  std::ifstream fs_file((name + ".frag").c_str());
  if (fs_file.is_open()) {
    fs_str = std::string((std::istreambuf_iterator<char>(fs_file)),
                         std::istreambuf_iterator<char>());
    fs_source = fs_str.c_str();
  }

  if (!vs_source && !fs_source) {
    std::cout << "Error when loading shader " << name << ": files not found."
              << std::endl;
    return false;
  }

  m_program = glCreateProgram();
  m_vs = glCreateShader(GL_VERTEX_SHADER);
  m_fs = glCreateShader(GL_FRAGMENT_SHADER);
  printProgramInfoLog(m_program, name, "creating shader");

  glShaderSource(m_vs, 1, &vs_source, NULL);
  printShaderInfoLog(m_vs, name, "loading vertex shader source");
  glShaderSource(m_fs, 1, &fs_source, NULL);
  printShaderInfoLog(m_fs, name, "loading fragment shader source");

  glCompileShader(m_vs);
  printShaderInfoLog(m_vs, name, "compiling vertex shader");
  glCompileShader(m_fs);
  printShaderInfoLog(m_fs, name, "compiling fragment shader");

  glAttachShader(m_program, m_vs);
  printShaderInfoLog(m_vs, name, "attaching vertex shader to program");
  glAttachShader(m_program, m_fs);
  printShaderInfoLog(m_fs, name, "attaching fragment shader to program");

  for (std::map<std::string, GLuint>::iterator i = m_vertex_attributes.begin();
       i != m_vertex_attributes.end(); i++)
    glBindAttribLocation(m_program, i->second, i->first.c_str());

  glLinkProgram(m_program);
  printProgramInfoLog(m_program, name, "linking");

  GLint ls = 0;
  glGetProgramiv(m_program, GL_LINK_STATUS, &ls);
  return (ls == GL_TRUE);
}

void Shader::setVertexAttribName(const std::string& name, GLuint index) {
  m_vertex_attributes[name] = index;
}

GLint Shader::getUniformLocation(const std::string& name) {
  std::map<std::string, GLint>::iterator entry = m_uniform_locations.find(name);
  if (entry != m_uniform_locations.end()) {
    // this entry is found in the cache
    return entry->second;
  } else {
    // perform the query in OpenGL, and cache the result
    GLint location = glGetUniformLocation(m_program, name.c_str());
    m_uniform_locations[name] = location;
    return location;
  }
}

GLuint Shader::getVertexAttribLocation(const std::string& name) {
  std::map<std::string, GLuint>::iterator entry =
      m_vertex_attributes.find(name);

  // this entry must have been set by setVertexAttribName()
  assert(entry != m_vertex_attributes.end());
  return entry->second;
}

void Shader::printProgramInfoLog(GLuint prog, const std::string& name,
                                 const std::string& stage) {
  GLint infologLength = 0, charsWritten = 0;
  glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infologLength);
  if (infologLength > 2) {
    GLchar* infoLog = new GLchar[infologLength];
    glGetProgramInfoLog(prog, infologLength, &charsWritten, infoLog);
    std::cerr << "Shader " << name << " at " << stage << ": " << std::endl
              << infoLog << std::endl;
    delete infoLog;
  }
}

void Shader::printShaderInfoLog(GLuint shader, const std::string& name,
                                const std::string& stage) {
  GLint infologLength = 0, charsWritten = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infologLength);
  if (infologLength > 2) {
    GLchar* infoLog = new GLchar[infologLength];
    glGetShaderInfoLog(shader, infologLength, &charsWritten, infoLog);
    std::cerr << "Shader " << name << " at " << stage << ": " << std::endl
              << infoLog << std::endl;
    delete infoLog;
  }
}
