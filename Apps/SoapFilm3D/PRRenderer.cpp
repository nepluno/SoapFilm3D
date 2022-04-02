//
//  PRRenderer.cpp
//  MultiTracker
//
//  Created by Fang Da on 1/15/15.
//
//

#include "PRRenderer.h"
#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "textures/beach.h"

PRRenderer::PRRenderer(VS3D *vs)
    : m_vs(vs), m_shader_bubble(), m_shader_env(), m_tex_env(0) {
#ifndef NO_SHADER
  m_shader_bubble.setVertexAttribName("a_position", 0);
  m_shader_bubble.setVertexAttribName("a_normal", 1);
  m_shader_bubble.loadFromFile("Apps/SoapFilm3D/PR");

  m_shader_env.setVertexAttribName("a_position", 0);
  m_shader_env.loadFromFile("Apps/SoapFilm3D/env");
#endif
  create_cube_map(m_tex_env);
}

namespace {
class FaceOrderComp {
 public:
  bool operator()(const std::pair<size_t, double> &f1,
                  const std::pair<size_t, double> &f2) const {
    return f1.second < f2.second;
  }
};
}  // namespace

void PRRenderer::render() {
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_CUBE_MAP, m_tex_env);

  glClearDepth(1.0);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  Mat4 MV;
  Mat4 PJ;
  {
    float mv[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mv);
    float pj[16];
    glGetFloatv(GL_PROJECTION_MATRIX, pj);
    MV << mv[0], mv[4], mv[8], mv[12], mv[1], mv[5], mv[9], mv[13], mv[2],
        mv[6], mv[10], mv[14], mv[3], mv[7], mv[11], mv[15];
    PJ << pj[0], pj[4], pj[8], pj[12], pj[1], pj[5], pj[9], pj[13], pj[2],
        pj[6], pj[10], pj[14], pj[3], pj[7], pj[11], pj[15];
  }
  Mat4 MVP = PJ * MV;

  Vec4 cam_pos_h = MV.inverse() * Vec4(0, 0, 1, 1);
  cam_pos_h /= cam_pos_h.w();
  Vec3 cam_pos(cam_pos_h.x(), cam_pos_h.y(), cam_pos_h.z());

  // render env map
  float env_vb[] = {
      -10.0f, 10.0f,  -10.0f, -10.0f, -10.0f, -10.0f, 10.0f,  -10.0f, -10.0f,
      10.0f,  -10.0f, -10.0f, 10.0f,  10.0f,  -10.0f, -10.0f, 10.0f,  -10.0f,

      -10.0f, -10.0f, 10.0f,  -10.0f, -10.0f, -10.0f, -10.0f, 10.0f,  -10.0f,
      -10.0f, 10.0f,  -10.0f, -10.0f, 10.0f,  10.0f,  -10.0f, -10.0f, 10.0f,

      10.0f,  -10.0f, -10.0f, 10.0f,  -10.0f, 10.0f,  10.0f,  10.0f,  10.0f,
      10.0f,  10.0f,  10.0f,  10.0f,  10.0f,  -10.0f, 10.0f,  -10.0f, -10.0f,

      -10.0f, -10.0f, 10.0f,  -10.0f, 10.0f,  10.0f,  10.0f,  10.0f,  10.0f,
      10.0f,  10.0f,  10.0f,  10.0f,  -10.0f, 10.0f,  -10.0f, -10.0f, 10.0f,

      -10.0f, 10.0f,  -10.0f, 10.0f,  10.0f,  -10.0f, 10.0f,  10.0f,  10.0f,
      10.0f,  10.0f,  10.0f,  -10.0f, 10.0f,  10.0f,  -10.0f, 10.0f,  -10.0f,

      -10.0f, -10.0f, -10.0f, -10.0f, -10.0f, 10.0f,  10.0f,  -10.0f, -10.0f,
      10.0f,  -10.0f, -10.0f, -10.0f, -10.0f, 10.0f,  10.0f,  -10.0f, 10.0f};

  for (int i = 0; i < 36; i++) {
    env_vb[i * 3 + 0] += cam_pos[0];
    env_vb[i * 3 + 1] += cam_pos[1];
    env_vb[i * 3 + 2] += cam_pos[2];
  }

  m_shader_env.activate();

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glCullFace(GL_FRONT_AND_BACK);

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, env_vb);

  glUniformMatrix4fv(m_shader_env.getUniformLocation("u_mat_mvp"), 1, GL_FALSE,
                     MVP.data());

  m_shader_env.setUniform("u_tex_env", 0);
  m_shader_env.setUniform("u_camera_pos", cam_pos);

  glDrawArrays(GL_TRIANGLES, 0, 36);

  m_shader_env.deactivate();

  // render the mesh, back-to-front
  std::vector<std::pair<size_t, double> > face_ordering;
  for (size_t i = 0; i < mesh().nt(); i++) {
    LosTopos::Vec3st t = mesh().get_triangle(i);
    if (m_vs->surfTrack()->triangle_is_all_solid(i)) continue;
    Vec3 x0(m_vs->pos(t[0]).x(), m_vs->pos(t[0]).y(), m_vs->pos(t[0]).z());
    Vec3 x1(m_vs->pos(t[1]).x(), m_vs->pos(t[1]).y(), m_vs->pos(t[1]).z());
    Vec3 x2(m_vs->pos(t[2]).x(), m_vs->pos(t[2]).y(), m_vs->pos(t[2]).z());
    Vec3 c = (x0 + x1 + x2) / 3;
    face_ordering.push_back(std::pair<size_t, double>(i, c.dot(cam_pos)));
  }
  size_t nt = face_ordering.size();

  FaceOrderComp comp;
  std::sort(face_ordering.begin(), face_ordering.end(), comp);

  float *vb = new float[nt * 9];
  float *nb = new float[nt * 9];

  for (size_t ii = 0; ii < nt; ii++) {
    size_t i = face_ordering[ii].first;
    LosTopos::Vec3st t = mesh().get_triangle(i);
    Vec3 x0(m_vs->pos(t[0]).x(), m_vs->pos(t[0]).y(), m_vs->pos(t[0]).z());
    Vec3 x1(m_vs->pos(t[1]).x(), m_vs->pos(t[1]).y(), m_vs->pos(t[1]).z());
    Vec3 x2(m_vs->pos(t[2]).x(), m_vs->pos(t[2]).y(), m_vs->pos(t[2]).z());
    Vec3 n[3];
    LosTopos::Vec2i l = mesh().get_triangle_label(i);
    Vec2i rp = (l[0] < l[1] ? Vec2i(l[0], l[1]) : Vec2i(l[1], l[0]));
    for (int j = 0; j < 3; j++) {
      n[j] = Vec3(0, 0, 0);
      for (size_t k = 0; k < mesh().m_vertex_to_triangle_map[t[j]].size();
           k++) {
        LosTopos::Vec2i ln =
            mesh().get_triangle_label(mesh().m_vertex_to_triangle_map[t[j]][k]);
        Vec2i rpn = (ln[0] < ln[1] ? Vec2i(ln[0], ln[1]) : Vec2i(ln[1], ln[0]));
        if (rp == rpn) {
          LosTopos::Vec3st tn =
              mesh().get_triangle(mesh().m_vertex_to_triangle_map[t[j]][k]);
          Vec3 xn0(m_vs->pos(tn[0]).x(), m_vs->pos(tn[0]).y(),
                   m_vs->pos(tn[0]).z());
          Vec3 xn1(m_vs->pos(tn[1]).x(), m_vs->pos(tn[1]).y(),
                   m_vs->pos(tn[1]).z());
          Vec3 xn2(m_vs->pos(tn[2]).x(), m_vs->pos(tn[2]).y(),
                   m_vs->pos(tn[2]).z());
          n[j] += (xn1 - xn0).cross(xn2 - xn0) * (ln[0] == l[0] ? 1 : -1);
        }
      }
      n[j].normalize();
      n[j] *= 3;
    }

    vb[ii * 9 + 0] = x0[0];
    vb[ii * 9 + 1] = x0[1];
    vb[ii * 9 + 2] = x0[2];
    vb[ii * 9 + 3] = x1[0];
    vb[ii * 9 + 4] = x1[1];
    vb[ii * 9 + 5] = x1[2];
    vb[ii * 9 + 6] = x2[0];
    vb[ii * 9 + 7] = x2[1];
    vb[ii * 9 + 8] = x2[2];
    nb[ii * 9 + 0] = n[0][0];
    nb[ii * 9 + 1] = n[0][1];
    nb[ii * 9 + 2] = n[0][2];
    nb[ii * 9 + 3] = n[1][0];
    nb[ii * 9 + 4] = n[1][1];
    nb[ii * 9 + 5] = n[1][2];
    nb[ii * 9 + 6] = n[2][0];
    nb[ii * 9 + 7] = n[2][1];
    nb[ii * 9 + 8] = n[2][2];
  }

  m_shader_bubble.activate();

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glCullFace(GL_FRONT_AND_BACK);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, vb);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nb);

  m_shader_bubble.setUniform("u_camera_pos", cam_pos);

  glUniformMatrix4fv(m_shader_bubble.getUniformLocation("u_mat_mvp"), 1,
                     GL_FALSE, MVP.data());

  glUniform3f(m_shader_bubble.getUniformLocation("u_light_direction"),
              0.26726124191, 0.53452248382, -0.80178372573);
  glUniform4f(m_shader_bubble.getUniformLocation("u_light_ambient"), 0.3, 0.3,
              0.3, 1);
  glUniform4f(m_shader_bubble.getUniformLocation("u_light_diffuse"), 1, 1, 1,
              1);
  glUniform4f(m_shader_bubble.getUniformLocation("u_light_specular"), 1, 1, 1,
              1);
  glUniform1f(m_shader_bubble.getUniformLocation("u_light_specularity"), 120);

  m_shader_bubble.setUniform("u_tex_env", 0);

  glDrawArrays(GL_TRIANGLES, 0, nt * 3);

  m_shader_bubble.deactivate();

  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);

  delete[] vb;
  delete[] nb;
  assert(glGetError() == GL_NO_ERROR);
}

bool PRRenderer::create_cube_map(GLuint &tex) {
  // generate a cube-map texture to hold all the sides
  glActiveTexture(GL_TEXTURE0);
  glGenTextures(1, &tex);

  // load each image and copy into a side of the cube-map texture
  bool success;
  success = load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
                               texture::beach::bottom,
                               sizeof(texture::beach::bottom));
  assert(success);
  success =
      load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
                         texture::beach::top, sizeof(texture::beach::top));
  assert(success);
  success =
      load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
                         texture::beach::front, sizeof(texture::beach::front));
  assert(success);
  success =
      load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
                         texture::beach::back, sizeof(texture::beach::back));
  assert(success);
  success =
      load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
                         texture::beach::left, sizeof(texture::beach::left));
  assert(success);
  success =
      load_cube_map_side(tex, GL_TEXTURE_CUBE_MAP_POSITIVE_X,
                         texture::beach::right, sizeof(texture::beach::right));
  assert(success);

  // format cube map texture
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  return true;
}

bool PRRenderer::load_cube_map_side(GLuint tex, GLenum side,
                                    const unsigned char *content,
                                    unsigned int content_size) {
  int w = 256;
  int h = 256;
  int channels = 3;
  unsigned char *image_data =
      stbi_load_from_memory(content, content_size, &w, &h, &channels, 0);

  // copy image data into 'target' side of cube map
  glBindTexture(GL_TEXTURE_CUBE_MAP, tex);
  glTexImage2D(side, 0, GL_RGBA, w, h, 0, channels == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE,
               image_data);
  stbi_image_free(image_data);

  return true;
}
