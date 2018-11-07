R"(
#version 400

uniform mat4 eye_projection;
uniform mat4 world_to_eye;

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec2 vertex_texture;
//in vec3 vertex_position;
//in vec2 vertex_texture;

out vec2 v_texture;
out vec3 v_eye_to_position;

void main()
{
	v_texture = vertex_texture;
	v_eye_to_position = vec3( world_to_eye * vec4( vertex_position , 1.0 ) );
	gl_Position = eye_projection * world_to_eye * vec4( vertex_position , 1.0 );
}
)"