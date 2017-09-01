//
//  PR.frag

//precision mediump float;

uniform sampler2D u_tex_depth;
uniform samplerCube u_tex_env;

uniform mat4 u_mat_mvp;
uniform mat4 u_mat_mv;

uniform vec3 u_light_direction;
uniform vec4 u_light_diffuse;
uniform vec4 u_light_ambient;
uniform vec4 u_light_specular;
uniform float u_light_specularity;
uniform vec3 u_camera_pos;
uniform float u_ssao_coef;

uniform vec3 u_ssao_samples[128];

varying vec4 v_position_world;
//varying vec4 v_position_clip;
varying vec3 v_normal;

void main()
{
    vec3 normal = normalize(v_normal);

    // phong shading
    vec4 shading_diffuse = u_light_diffuse * clamp(dot(normal, -u_light_direction), 0.0, 1.0);

    vec3 viewvec = normalize(u_camera_pos - v_position_world.xyz);
    vec3 lightvec = -u_light_direction;
    vec3 halfway = normalize(viewvec + lightvec);
    vec4 shading_specular = u_light_specular * pow(clamp(dot(halfway, normal), 0.0, 1.0), u_light_specularity);

    vec4 shading_ambient = vec4(0.3, 0.3, 0.3, 0.3);
    
    // env map reflection
    vec3 reflectiondir = reflect(-viewvec, normal);
    vec4 reflection = textureCube(u_tex_env, reflectiondir);
    
    // fresnel effect
    float angle = abs(dot(viewvec, normal));
    float alpha = 0.3 + 0.9 * pow(1.0 - angle, 2.0);
    
//    gl_FragColor = (shading_diffuse * 0.6 + shading_specular * 0.5 + shading_ambient) * 0.3 + reflection * 0.7;;
    gl_FragColor = vec4(reflection.xyz * 1.5, alpha);
//    gl_FragColor = vec4(alpha, 0.0, 0.0, 1.0);

}



