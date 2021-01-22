#version 330

in vec3 pos;
out vec4 fragment_color;
uniform int type;

void main(){
	if (type == 0) {
		fragment_color = vec4(1, 0, 0, 1); 
	} else if(type == 1) {
		fragment_color=vec4(0, 1 , 0, 1);
	} else if (type == 2) {
		fragment_color = vec4(1-pos, 1);
	}
}