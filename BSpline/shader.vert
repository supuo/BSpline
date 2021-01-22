#version 330

layout (location = 0) in vec3 vPosition;
out vec3 pos;

uniform mat4 transform;

void main() {
    pos = vPosition;
    gl_Position = transform * vec4(vPosition, 1.0);
}