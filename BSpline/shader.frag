#version 330

out vec4 fragment_color;
uniform int shading_type;

void main(){
    if (shading_type == 2) {
		fragment_color = vec4(1, 1, 1, 1);
	}
    else if(shading_type == 1) {
		fragment_color=vec4(0, 1 , 0, 1);
	}
	else if (shading_type == 0) {
		fragment_color = vec4(1, 0, 0, 1); 
	}
}