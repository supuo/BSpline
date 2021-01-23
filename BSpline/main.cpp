#include <vector>
#include <iostream>

#include <glm/glm.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "BSpline.h"
#include "Fitter.h"
#include "Shader.h"
#include "Camera.h"

using namespace std;
using glm::vec3;

namespace {
	enum Type {
		Curve, Surface
	} type;

	const int width = 800;
	const int height = 800;

	Fitter fitter;
	bool change = true;
	const double step = 0.01;
	unsigned vao[4], vbo[3], ebo[4], indexSize[4];

	Camera camera(vec3(0, 0, 3));
	double lastX;
	double lastY;
	bool firstMouse = true;

	double err = 1e-8;
}

void mouseCallback(GLFWwindow* window, double xpos, double ypos) {
	if (firstMouse) {
		lastX = (float)xpos;
		lastY = (float)ypos;
		firstMouse = false;
	}

	float xoffset = (float)xpos - lastX;
	float yoffset = lastY - (float)ypos;

	lastX = (float)xpos;
	lastY = (float)ypos;

	camera.processMouseMovement(xoffset, yoffset);
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	camera.processMouseScroll((float)yoffset);
}

GLFWwindow* initOpenGL() {
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(width, height, "BSpline surface fit", nullptr, nullptr);
	if (window == nullptr) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		exit(-1);
	}
	glfwMakeContextCurrent(window);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetScrollCallback(window, scrollCallback);

	if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
		std::cout << "Failed to initialize GLAD" << std::endl;
		exit(-1);
	}
	glViewport(0, 0, width, height);
	glEnable(GL_DEPTH_TEST);
	glGenVertexArrays(4, vao);
	glGenBuffers(3, vbo);
	glGenBuffers(4, ebo);
	glPointSize(5);
	return window;
}

void processInput(GLFWwindow* window) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Forward);
	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Backward);
	}
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Left);
	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Right);
	}
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Up);
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Down);
	}
}

void render(const Shader& shader) {
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	shader.use();
	glm::mat4 view = camera.getViewMatrix();
	glm::mat4 projection = glm::perspective(glm::radians(camera.getZoom()), 1.f * width / height, 0.01f, 100.0f);
	glm::mat4 transform = projection * view;
	shader.setMat4("transform", transform);
	// DataPoints
	shader.setInt("type", 0);
	glBindVertexArray(vao[0]);
	glDrawElements(GL_POINTS, indexSize[0], GL_UNSIGNED_INT, 0);
	// controlPoints
	shader.setInt("type", 1);
	glBindVertexArray(vao[1]);
	glDrawElements(GL_POINTS, indexSize[1], GL_UNSIGNED_INT, 0);
	// controlLines
	shader.setInt("type", 1);
	glBindVertexArray(vao[2]);
	glDrawElements(GL_LINES, indexSize[2], GL_UNSIGNED_INT, 0);
	// obj
	shader.setInt("type", 2);
	glBindVertexArray(vao[3]);
	if (type == Curve) {
		glDrawElements(GL_LINES, indexSize[3], GL_UNSIGNED_INT, 0);
	}
	if (type == Surface) {
		glDrawElements(GL_TRIANGLES, indexSize[3], GL_UNSIGNED_INT, 0);
	}
}

vector<vec3> generateCircle(int n) {
	vector<vec3> dataPoints;
	double radius = 0.5;
	double degreeStep = 360. / n;
	for (int i = 0; i < n; ++i) {
		double radians = glm::radians(degreeStep * i);
		dataPoints.emplace_back(radius * cos(radians), radius * sin(radians), 0);
	}
	return dataPoints;
}

vector<vector<vec3>> generateCylinder(int n) {
	vector<vector<vec3>> dataPoints;
	double radius = 0.5;
	double degreeStep = 360. / n, heightStep = 1. / n;
	for (int i = 0; i <= n; ++i) {
		vector<vec3> row;
		double radian = glm::radians(degreeStep * i);
		for (int j = 0; j <= n; ++j) {
			double height = -0.5 + heightStep * j;
			row.emplace_back(radius * sin(radian), radius * cos(radian), height);
		}
		dataPoints.emplace_back(row);
	}
	return dataPoints;

}

vector<vector<vec3>> generateSphere(int n) {
	vector<vector<vec3>> dataPoints;
	double radius = 0.5;
	double thetaStep = 360 / (1. * n);
	double phiStep = 180 / (1. * n);
	for (int i = 0; i <= n; ++i) {
		vector<vec3> row;
		double phi = glm::radians(phiStep * i);
		for (int j = 0; j <= n; ++j) {
			double theta = glm::radians(thetaStep * j);
			row.emplace_back(radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi));
		}
		dataPoints.emplace_back(row);
	}
	return dataPoints;
}

void bspCurve(int p, const vector<vec3>& dataPoints, BSpline::BSplineType bspType) {
	if (dataPoints.empty()) return;
	type = Curve;
	std::vector<vec3> controlPoints;
	auto bc = fitter.interpolateCurve(p, dataPoints, controlPoints, bspType);
	int n = bc.n;
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
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, controlPoints.size() * sizeof(vec3), controlPoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[1]);
	std::vector<int> controlPointIndex;
	for (int i = 0; i <= n; i++) {
		controlPointIndex.push_back(i);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	             controlPointIndex.size() * sizeof(int),
	             controlPointIndex.data(),
	             GL_STATIC_DRAW);
	indexSize[1] = static_cast<int>(controlPointIndex.size());
	// controlLines
	glBindVertexArray(vao[2]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[2]);
	std::vector<int> lineIndex;
	for (int i = 0; i <= n - 1; i++) {
		lineIndex.push_back(i);
		lineIndex.push_back(i + 1);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, lineIndex.size() * sizeof(int), lineIndex.data(), GL_STATIC_DRAW);
	indexSize[2] = static_cast<int>(lineIndex.size());
	// BSplineCurve
	glBindVertexArray(vao[3]);
	std::vector<vec3> curvePoints;
	for (double u = 0; u < 1 + err; u += step) {
		curvePoints.emplace_back(bc(controlPoints, u));
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, curvePoints.size() * sizeof(vec3), curvePoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[3]);
	std::vector<int> curveIndex;
	for (int i = 0; i < curvePoints.size() - 1; i++) {
		curveIndex.push_back(i);
		curveIndex.push_back(i + 1);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, curveIndex.size() * sizeof(int), curveIndex.data(), GL_STATIC_DRAW);
	indexSize[3] = static_cast<int>(curveIndex.size());
}

void bspSurface(int p,
                int q,
                const vector<vector<vec3>>& dataPoints,
                BSpline::BSplineType utype,
                BSpline::BSplineType vtype) {
	if (dataPoints.empty()) return;
	type = Surface;
	std::vector<vector<vec3>> controlPoints;
	auto bs = fitter.interpolateSurface(p, q, dataPoints, controlPoints, utype, vtype);
	int m = static_cast<int>(dataPoints.size()), n = static_cast<int>(dataPoints[0].size());
	// dataPoints
	glBindVertexArray(vao[0]);
	std::vector<vec3> data(m * n);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			data[i * n + j] = dataPoints[i][j];
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(vec3), data.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), nullptr);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[0]);
	std::vector<int> pointIndex(m * n);
	for (int i = 0; i < m * n; i++) {
		pointIndex[i] = i;
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, pointIndex.size() * sizeof(int), pointIndex.data(), GL_STATIC_DRAW);
	indexSize[0] = m * n;
	// controlPoints
	glBindVertexArray(vao[1]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	vector<vec3> ctrPoints;
	ctrPoints.reserve(m * n);
	for (auto& row : controlPoints) {
		for (auto& point : row) {
			ctrPoints.push_back(point);
		}
	}
	glBufferData(GL_ARRAY_BUFFER, ctrPoints.size() * sizeof(vec3), ctrPoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[1]);
	std::vector<int> controlPointIndex;
	for (int i = 0; i < m * n; i++) {
		controlPointIndex.push_back(i);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	             controlPointIndex.size() * sizeof(int),
	             controlPointIndex.data(),
	             GL_STATIC_DRAW);
	indexSize[1] = static_cast<int>(controlPointIndex.size());
	// controlLines
	glBindVertexArray(vao[2]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[2]);
	std::vector<int> lineIndex;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n - 1; ++j) {
			lineIndex.push_back(i * m + j);
			lineIndex.push_back(i * m + j + 1);
		}
	}
	for (int i = 0; i < m - 1; ++i) {
		for (int j = 0; j < n; ++j) {
			lineIndex.push_back(i * m + j);
			lineIndex.push_back((i + 1) * m + j);
		}
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, lineIndex.size() * sizeof(int), lineIndex.data(), GL_STATIC_DRAW);
	indexSize[2] = static_cast<int>(lineIndex.size());
	// BSplineSurface
	glBindVertexArray(vao[3]);
	std::vector<vec3> surfacePoints;
	int count = 0;
	for (double u = 0; u < 1 + err; u += step) {
		for (double v = 0; v < 1 + err; v += step) {
			surfacePoints.emplace_back(bs(controlPoints, u, v));
		}
		count++;
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, surfacePoints.size() * sizeof(vec3), surfacePoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[3]);
	std::vector<int> surfaceIndex;
	for (int i = 0; i < count - 1; ++i) {
		for (int j = 0; j < count - 1; ++j) {
			int idx = i * count + j;
			surfaceIndex.push_back(idx);
			surfaceIndex.push_back(idx + count);
			surfaceIndex.push_back(idx + 1);
			surfaceIndex.push_back(idx + count);
			surfaceIndex.push_back(idx + count + 1);
			surfaceIndex.push_back(idx + 1);
		}
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, surfaceIndex.size() * sizeof(int), surfaceIndex.data(), GL_STATIC_DRAW);
	indexSize[3] = static_cast<int>(surfaceIndex.size());
}

int main() {
	GLFWwindow* window = initOpenGL();
	Shader shader("shader.vert", "shader.frag");
	// auto curveDataPoints = generateCircle(10);
	// auto surfaceDataPoints = generateCylinder(8);
	auto surfaceDataPoints = generateSphere(10);
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		processInput(window);
		if (change) {
			change = false;
			// bspCurve(3, curveDataPoints);
			bspSurface(3, 3, surfaceDataPoints, BSpline::BSplineType::Clamped, BSpline::BSplineType::Closed);
		}
		render(shader);
		glfwSwapBuffers(window);
	}
	glfwTerminate();
	return 0;
}
