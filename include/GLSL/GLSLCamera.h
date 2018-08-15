#ifndef GLSL_CAMERA_INCLUDED
#define GLSL_CAMERA_INCLUDED

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using glm::vec2;
using glm::vec3;
using glm::ivec3;
using glm::vec4;
using glm::mat4;
using glm::mat3;

class Camera
{
public:
	float heightAngle;
	float aspectRatio;
	float nearest_plane;
	float farthest_plane;
	vec3 position;
	vec3 direction;
	vec3 up;
	vec3 right;
	mat4 projection_matrix;
	mat4 world_to_camera;
	void rotateUp(vec3 center, float angle);
	void rotateRight(vec3 center, float angle);
	void rotateDirection(vec3 center, float angle);
	void moveForward(float dist);
	void moveRight(float dist);
	void moveUp(float dist);
	void updateTransformations();
};

vec3 CentralizedRotation(const vec3 center, const vec3 direction, const vec3 position, const double angle)
{
	vec3 centralized_position = position - center;
	vec3 v1 = glm::cross(direction, centralized_position);
	vec3 v2 = glm::cross(direction, v1);
	vec3 rotated_centralized_position = centralized_position + v1*(float(sin(angle))) + v2*(float(1.f - cos(angle)));

	vec3 result = center + rotated_centralized_position;

	return result;
}

void Camera::rotateUp(vec3 center, float angle)
{
	right = CentralizedRotation(center, up, right, angle);
	right = glm::normalize(right);
	direction = CentralizedRotation(center, up, direction, angle);
	direction = glm::normalize(direction);
	position = CentralizedRotation(center, up, position, angle);
}

void Camera::rotateDirection(vec3 center, float angle)
{
	right = CentralizedRotation(center, direction, right, angle);
	right = glm::normalize(right);
	up = CentralizedRotation(center, direction, up, angle);
	up = glm::normalize(up);
	position = CentralizedRotation(center, direction, position, angle);
}

void Camera::rotateRight( vec3 center , float angle )
{
	up = CentralizedRotation( center , right , up , angle );
	up = glm::normalize(up);
	direction = CentralizedRotation(center, right, direction, angle);
	direction = glm::normalize(direction);
	position = CentralizedRotation(center, right, position, angle);
}
void Camera::moveForward(float dist) {
	position += (direction*dist);
}
void Camera::moveRight(float dist) {
	position += (right*dist);
}
void Camera::moveUp(float dist) {
	position += (up*dist);
}

void  Camera::updateTransformations( void )
{
	projection_matrix = glm::perspective( heightAngle , aspectRatio , nearest_plane , farthest_plane );
	world_to_camera = glm::lookAt( position , position + direction , up );
}

#endif //GLSL_CAMERA_INCLUDED