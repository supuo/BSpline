#include "Camera.h"

Camera::Camera(glm::vec3 _position):
	position(_position) {
	updateCameraVectors();
}

glm::mat4 Camera::getViewMatrix() const {
	return lookAt(position, position + front, up);
}

void Camera::processKeyboard(CameraMovement direction) {
	if (direction == CameraMovement::FORWARD) position += front * movementSpeed;
	if (direction == CameraMovement::BACKWARD) position -= front * movementSpeed;
	if (direction == CameraMovement::LEFT) position -= right * movementSpeed;
	if (direction == CameraMovement::RIGHT) position += right * movementSpeed;
}

void Camera::processMouseMovement(float xoffset, float yoffset) {
	xoffset *= mouseSensitivity;
	yoffset *= mouseSensitivity;

	yaw += xoffset;
	pitch += yoffset;

	if (pitch > 89.0f) pitch = 89.0f;
	if (pitch < -89.0f) pitch = -89.0f;

	updateCameraVectors();
}

void Camera::processMouseScroll(float yoffset) {
	zoom -= static_cast<float>(yoffset);
	if (zoom < 1.0f) zoom = 1.0f;
	if (zoom > 45.0f) zoom = 45.0f;
}

float Camera::getZoom() const {
	return static_cast<float>(zoom);
}

void Camera::updateCameraVectors() {
	front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
	front.y = sin(glm::radians(pitch));
	front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
	front = normalize(front);
	right = normalize(cross(front, worldUp));
	up = normalize(cross(right, front));
}
