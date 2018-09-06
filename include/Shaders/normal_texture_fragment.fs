R"(
#version 400

uniform vec3 light_direction;

uniform vec3 light_ambient;
uniform vec3 light_diffuse;
uniform vec3 light_specular;

uniform vec3 shape_ambient;
uniform vec3 shape_diffuse;
uniform vec3 shape_specular;
uniform float shape_specular_shininess;

uniform sampler2D normal_texture;
in vec2 v_texture;
in vec3 v_eye_to_position;
layout (location = 0) out vec4 FragColor;

void main()
{
	vec3 surface_normal = normalize( vec3( texture( normal_texture , v_texture ) ) );
	vec3 position_to_light = normalize( -light_direction );
	vec3 max_color = vec3( 1. , 1. , 1. );

	// Ambient component
	vec3 lighting_color = light_ambient * shape_ambient;
   
   float normal_dot_light = dot( position_to_light , surface_normal );
	if( normal_dot_light>0.f )
	{
		// Diffuse component
		lighting_color += normal_dot_light * ( light_diffuse * shape_diffuse );

		// Specular component
		vec3 position_to_eye = normalize( -v_eye_to_position );
        vec3 specular_direction = 2.f * surface_normal * normal_dot_light - position_to_light;
        float specular_attenuation = pow( max( dot( specular_direction , position_to_eye ) , 0.0 ) , shape_specular_shininess );
		lighting_color += specular_attenuation * ( light_specular * shape_specular );
	}

	FragColor = vec4( min( lighting_color , max_color ) , 1.0 );
}
)"