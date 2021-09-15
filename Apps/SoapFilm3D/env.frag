//
//  env.frag

// precision mediump float;

uniform samplerCube u_tex_env;

uniform mat4 u_mat_mvp;
uniform vec3 u_camera_pos;

varying vec4 v_position_world;

void main() {
  vec3 viewvec = v_position_world.xyz / v_position_world.w - u_camera_pos;
  gl_FragColor = textureCube(u_tex_env, viewvec);
}
