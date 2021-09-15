//
//  PR.vert

// precision mediump float;

attribute vec4 a_position;
attribute vec3 a_normal;

uniform mat4 u_mat_mvp;
uniform mat4 u_mat_mv;

varying vec4 v_position_world;
// varying vec4 v_position_clip;
varying vec3 v_normal;

void main() {
  v_position_world = a_position;
  v_normal = a_normal;

  gl_Position = u_mat_mvp * a_position;
  //    v_position_clip = gl_Position / gl_Position.w;
}
