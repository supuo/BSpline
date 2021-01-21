#version 330

layout (location = 0) in vec3 vPosition;
out vec3 world_normal;

void main() {
	gl_Position = vec4(vPosition, 1.0);
}