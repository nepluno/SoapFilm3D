//
//  env.vert

// precision mediump float;

attribute vec4 a_position;

uniform mat4 u_mat_mvp;

varying vec4 v_position_world;

void main() {
  v_position_world = a_position;

  gl_Position = u_mat_mvp * a_position;
}
