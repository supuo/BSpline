#include <vector>
#include <iostream>

#include <glm/glm.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "BSpline.h"
#include "Fitter.h"
#include "Shader.h"

namespace {
	const int width = 800;
	const int height = 800;

	std::vector<glm::vec3> dataPoints;
	BSpline* bsp = nullptr;

	bool change = true;
	unsigned vao[4], vbo[4], ebo[4], indexSize[4];

	double du = 0.001;
}

using namespace std;
using glm::vec3;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}
}

void render(const Shader& shader) {
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// 	float cos_theta = std::cos(theta);
	// 	float cos_phi = std::cos(phi);
	// 	float sin_theta = std::sin(theta);
	// 	float sin_phi = std::sin(phi);
	// 	vec3 dir{sin_theta * sin_phi, cos_theta, sin_theta * cos_phi};
	// 	vec3 eye = bound_center + 2.f * bound_radius * dir;
	// 	float tnear = 0.5f * bound_radius, tfar = 3.5f * bound_radius;
	// 	auto view = glm::lookAt(eye, bound_center, {0.f, 1.f, 0.f});
	// 	auto perspective = glm::perspective(glm::radians(60.f), WIDTH * 1.f / HEIGHT, tnear, tfar);
	shader.use();
	// glUniformMatrix4fv(glGetUniformLocation(program, "P"), 1, GL_FALSE, &(perspective[0][0]));
	// glUniformMatrix4fv(glGetUniformLocation(program, "V"), 1, GL_FALSE, &(view[0][0]));
	// dataPoints
	shader.setInt("shading_type", 0);
	glBindVertexArray(vao[0]);
	glDrawElements(GL_POINTS, indexSize[0], GL_UNSIGNED_INT, 0);
	// controlPoints
	shader.setInt("shading_type", 1);
	glBindVertexArray(vao[1]);
	glDrawElements(GL_POINTS, indexSize[1], GL_UNSIGNED_INT, 0);
	// controlLines
	shader.setInt("shading_type", 1);
	glBindVertexArray(vao[3]);
	glDrawElements(GL_LINES, indexSize[3], GL_UNSIGNED_INT, 0);
	// curve
	shader.setInt("shading_type", 2);
	glBindVertexArray(vao[2]);
	glDrawElements(GL_LINES, indexSize[2], GL_UNSIGNED_INT, 0);
}

void init_data(int dataNum) {
	double radius = 0.5;
	double degreeStep = 360. / dataNum;
	for (int i = 0; i <= dataNum; ++i) {
		double radians = glm::radians(degreeStep * i);
		dataPoints.emplace_back(radius * cos(radians), radius * sin(radians), 0);
	}
}

void drawCurve() {
	bsp = Fitter::interpolateCurve(2, dataPoints);
	// bsp->wrap();
	int n = bsp->n;
	int dataNum = static_cast<int>(dataPoints.size());
	// dataPoints
	glBindVertexArray(vao[0]);
	std::vector<vec3> data(dataNum);
	for (int i = 0; i < dataNum; i++) {
		data[i] = dataPoints[i];
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(vec3), data.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), nullptr);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[0]);
	std::vector<int> pointIndex(dataNum);
	for (int i = 0; i < dataNum; i++) {
		pointIndex[i] = i;
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, pointIndex.size() * sizeof(int), pointIndex.data(), GL_STATIC_DRAW);
	indexSize[0] = dataNum;
	// controlPoints
	glBindVertexArray(vao[1]);
	std::vector<vec3> controlPoints(n + 1);
	for (int i = 0; i <= n; i++) {
		controlPoints[i] = bsp->controlPoints[i];
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, controlPoints.size() * sizeof(vec3), controlPoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[1]);
	std::vector<int> controlPointIndex;
	for (int i = 0; i <= n; i++) {
		controlPointIndex.push_back(i);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, controlPointIndex.size() * sizeof(int), controlPointIndex.data(), GL_STATIC_DRAW);
	indexSize[1] = static_cast<int>(controlPointIndex.size());
	// controlLines
	glBindVertexArray(vao[3]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[3]);
	std::vector<int> lineIndex;
	for (int i = 0; i <= n - 1; i++) {
		lineIndex.push_back(i);
		lineIndex.push_back(i + 1);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, lineIndex.size() * sizeof(int), lineIndex.data(), GL_STATIC_DRAW);
	indexSize[3] = static_cast<int>(lineIndex.size());
	// BSplineCurve
	glBindVertexArray(vao[2]);
	std::vector<vec3> curvePoints;
	for (double u = 0; u < 1; u += du) {
		curvePoints.emplace_back((*bsp)(u));
	}
	curvePoints.emplace_back((*bsp)(1));
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, curvePoints.size() * sizeof(vec3), curvePoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[2]);
	std::vector<int> curveIndex;
	for (int i = 0; i < curvePoints.size() - 1; i++) {
		curveIndex.push_back(i);
		curveIndex.push_back(i + 1);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, curveIndex.size() * sizeof(int), curveIndex.data(), GL_STATIC_DRAW);
	indexSize[2] = static_cast<int>(curveIndex.size());
}

int main() {
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(width, height, "BSpline surface fit", nullptr, nullptr);
	if (window == nullptr) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}
	glViewport(0, 0, width, height);
	glEnable(GL_DEPTH_TEST);
	glGenVertexArrays(4, vao);
	glGenBuffers(4, vbo);
	glGenBuffers(4, ebo);
	glPointSize(4);
	Shader shader("shader.vert", "shader.frag");
	init_data(4);
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		if (change) {
			drawCurve();
			change = false;
		}
		render(shader);
		glfwSwapBuffers(window);
	}
	glfwTerminate();
	return 0;
}

// BSplineSurface generateSurface() {
// 	vector<vector<vec3>> dataPoints;
// 	 for (int i = 0; i < 10; ++i) {
// 	 	vector<vec3> row;
// 	 	for (int j = 0; j < 10; ++j) {
// 	 		row.emplace_back(i, j, rand() / (1. * RAND_MAX));
// 	 	}
// 	 	dataPoints.push_back(row);
// 	 }
// 	 BSplineSurface bSplineSurface(dataPoints, 2, 2);
// 	 auto&& result = bSplineSurface.generateSurface();
// 	 vector<vector<vec3>> surface(result.size(),
// 	                               vector<vec3>(result[0].size()));
// 	 for (int i = 0; i < result.size(); ++i) {
// 	 	for (int j = 0; j < result[0].size(); ++j) {
// 	 		surface[i][j] = vec3(result[i][j]);
// 	 	}
// 	 }
// 	return BSplineSurface();
// }
