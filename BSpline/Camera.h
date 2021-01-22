#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Camera {
public:
	enum class CameraMovement {
		FORWARD,
		BACKWARD,
		LEFT,
		RIGHT
	};

	explicit Camera(glm::vec3 _position = glm::vec3(0.0f, 0.0f, 0.0f));
	glm::mat4 getViewMatrix() const;
	void processKeyboard(CameraMovement direction);
	void processMouseMovement(float xoffset, float yoffset);
	void processMouseScroll(float yoffset);
	float getZoom() const;

private:
	void updateCameraVectors();

	glm::vec3 position{};
	glm::vec3 front{};
	glm::vec3 up{};
	glm::vec3 right{};
	glm::vec3 worldUp{0, 1, 0};
	float yaw = -90.f;
	float pitch = 0.f;
	float movementSpeed = 0.0125f;
	float mouseSensitivity = 0.025f;
	float zoom = 45.0f;
};
