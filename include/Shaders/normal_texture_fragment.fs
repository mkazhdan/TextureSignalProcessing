R"(
#version 400

uniform vec3 light_direction;
uniform vec3 light_diffuse;
uniform vec3 light_specular;
uniform vec3 light_ambient;
uniform float specular_falloff;

uniform sampler2D normal_texture;
in vec2 v_texture;
in vec3 v_eye_to_position;
layout (location = 0) out vec4 FragColor;

void main(){
   //vec3 surface_normal = normalize(vec3(texture(normal_texture, v_texture))*2.f-vec3(1.f,1.f,1.f));
   vec3 surface_normal = normalize(vec3(texture(normal_texture, v_texture)));
   vec3 position_to_light = normalize(-light_direction);
   vec3 lighting_color = light_ambient;
   
   if(dot(position_to_light,surface_normal) > 0.f){
		//Diffuse attenuation
		float diffuse_attenuation = max(dot(surface_normal,position_to_light), 0.0);
		lighting_color += diffuse_attenuation*light_diffuse;

		//Specular attenuation
		//vec3 position_to_eye = normalize(-v_eye_to_position);
        //vec3 specular_direction = 2.f*surface_normal*dot(position_to_light,surface_normal) - position_to_light;
        //float specular_attenuation = pow(max(dot(specular_direction,position_to_eye), 0.0),specular_falloff);
		//lighting_color += specular_attenuation*light_specular;
    }

   FragColor = vec4(lighting_color,1.0);
}
)"