#version 330

layout (location = 0) in vec3 vPosition;

uniform mat4 transform;

void main() {
    gl_Position = transform * vec4(vPosition, 1.0);
}