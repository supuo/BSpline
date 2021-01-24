// #pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" ) // Òþ²ØDOS

#include <vector>
#include <iostream>
#include <random>

#include <glm/glm.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "BSpline.h"
#include "Fitter.h"
#include "Shader.h"
#include "Camera.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

using namespace std;
using glm::vec3;

namespace {
	enum class GeometryType {
		Curve, Surface
	} gType;

	const int WIDTH = 800;
	const int HEIGHT = 800;
	const double PI = acos(-1);

	Fitter fitter(Fitter::ParametrizationMethod::Uniform, Fitter::KnotGenerationMethod::Uniform);
	vector<vec3>* curveDataPoints;
	vector<vector<vec3>>* surfaceDataPoints;

	int g_type = 0, p_type = 0, k_type = 0, bc_type = 0, bsu_type = 0, bsv_type = 0;
	int data_num = 10 , bc_p = 1, bs_p = 1, bs_q = 1, h = 3, e = 5, f = 5;
	BSpline::BSplineType bc = BSpline::BSplineType::Open, bsu = BSpline::BSplineType::Open, bsv =
		                     BSpline::BSplineType::Open;

	bool showData = true;
	bool showControl = true;
	bool showGeometry = true;

	const int maxData = 100;
	bool change = true;
	double radius = 0.5;
	const double step = 0.01;
	unsigned vao[4], vbo[3], ebo[4], indexSize[4];

	const int sampleNum = 100;

	Camera camera(vec3(0, 0, 3));
	bool firstMouse = true;
	int useCamera = 0;
	double lastX;
	double lastY;

	double error = 0;
	double err = 1e-8;
}

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_V && action == GLFW_PRESS) {
		firstMouse = true;
		useCamera ^= 1;
		if (!useCamera) {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		} else {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		}
	}
	if (key == GLFW_KEY_C && action == GLFW_PRESS) {
		showControl ^= 1;
	}
	if (key == GLFW_KEY_X && action == GLFW_PRESS) {
		showData ^= 1;
	}
	if (key == GLFW_KEY_Z && action == GLFW_PRESS) {
		showGeometry ^= 1;
	}
}

void mouseCallback(GLFWwindow* window, double xpos, double ypos) {
	if (!useCamera) return;

	if (firstMouse) {
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = static_cast<float>(xpos - lastX);
	float yoffset = static_cast<float>(lastY - ypos);

	lastX = xpos;
	lastY = ypos;

	camera.processMouseMovement(xoffset, yoffset);
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	if (!useCamera) return;
	camera.processMouseScroll((float)yoffset);
}

GLFWwindow* initOpenGL() {
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "BSpline fitting", nullptr, nullptr);
	if (window == nullptr) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		exit(-1);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetScrollCallback(window, scrollCallback);

	if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
		std::cout << "Failed to initialize GLAD" << std::endl;
		exit(-1);
	}
	glViewport(0, 0, WIDTH, HEIGHT);
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
	if (useCamera && glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Forward);
	}
	if (useCamera && glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Backward);
	}
	if (useCamera && glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Left);
	}
	if (useCamera && glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Right);
	}
	if (useCamera && glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Up);
	}
	if (useCamera && glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
		camera.processKeyboard(Camera::CameraMovement::Down);
	}
}

void render(const Shader& shader) {
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	shader.use();
	glm::mat4 view = camera.getViewMatrix();
	glm::mat4 projection = glm::perspective(glm::radians(camera.getZoom()), 1.f * WIDTH / HEIGHT, 0.01f, 100.0f);
	glm::mat4 transform = projection * view;
	shader.setMat4("transform", transform);
	if (showData) {
		// DataPoints
		shader.setInt("type", 0);
		glBindVertexArray(vao[0]);
		glDrawElements(GL_POINTS, indexSize[0], GL_UNSIGNED_INT, 0);
	}
	if (showControl) {
		// controlPoints
		shader.setInt("type", 1);
		glBindVertexArray(vao[1]);
		glDrawElements(GL_POINTS, indexSize[1], GL_UNSIGNED_INT, 0);
		// controlLines
		shader.setInt("type", 1);
		glBindVertexArray(vao[2]);
		glDrawElements(GL_LINES, indexSize[2], GL_UNSIGNED_INT, 0);
	}
	if (showGeometry) {
		shader.setInt("type", 2);
		glBindVertexArray(vao[3]);
		if (gType == GeometryType::Curve) {
			glDrawElements(GL_LINES, indexSize[3], GL_UNSIGNED_INT, 0);
		}
		if (gType == GeometryType::Surface) {
			glDrawElements(GL_TRIANGLES, indexSize[3], GL_UNSIGNED_INT, 0);
		}
	}
}

void generateCircle(int n) {
	delete curveDataPoints;
	curveDataPoints = new vector<vec3>;
	double degreeStep = 360. / n;
	for (int i = 0; i <= n; ++i) {
		double radians = glm::radians(degreeStep * i);
		curveDataPoints->emplace_back(radius * cos(radians), radius * sin(radians), 0);
	}
}

void generateCylinder(int n) {
	delete surfaceDataPoints;
	surfaceDataPoints = new vector<vector<vec3>>;
	double degreeStep = 360. / n, heightStep = 1. / n;
	for (int i = 0; i <= n; ++i) {
		vector<vec3> row;
		double radian = glm::radians(degreeStep * i);
		for (int j = 0; j <= n; ++j) {
			double height = -0.5 + heightStep * j;
			row.emplace_back(radius * sin(radian), radius * cos(radian), height);
		}
		surfaceDataPoints->emplace_back(row);
	}
}

void generateSphere(int n) {
	surfaceDataPoints = new vector<vector<vec3>>;
	double thetaStep = 360 / (1. * n);
	double phiStep = 360 / (1. * n);
	for (int i = 0; i <= n; ++i) {
		vector<vec3> row;
		double phi = glm::radians(phiStep * i);
		for (int j = 0; j <= n; ++j) {
			double theta = glm::radians(thetaStep * j);
			row.emplace_back(radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi));
		}
		surfaceDataPoints->emplace_back(row);
	}
}

void bspCurve(int p, const vector<vec3>& dataPoints, BSpline::BSplineType bspType) {
	std::vector<vec3> controlPoints;
	// auto bc = fitter.interpolateCurve(p, dataPoints, controlPoints, bspType);
	auto bc = fitter.approximateCurve(p, h, dataPoints, controlPoints);

	error = 0;
	std::default_random_engine random(static_cast<int>(time(nullptr)));
	uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
	for (int i = 0; i < sampleNum; ++i) {
		double u = uniformDistribution(random);
		auto samplePoint = vec3(radius * cos(u * 2 * PI), radius * sin(u * 2 * PI), 0);
		error += length(samplePoint - bc(controlPoints, u));
	}
	error /= radius * sampleNum;

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
	int h =static_cast<int>(controlPoints.size()) - 1;
	glBindVertexArray(vao[1]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, controlPoints.size() * sizeof(vec3), controlPoints.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[1]);
	std::vector<int> controlPointIndex;
	for (int i = 0; i <= h; i++) {
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
	for (int i = 0; i <= h - 1; i++) {
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
	std::vector<vector<vec3>> controlPoints;
	// auto bs = fitter.interpolateSurface(p, q, dataPoints, controlPoints, utype, vtype);
	auto bs = fitter.approximateSurface(p, q, e, f, dataPoints, controlPoints);
	
	error = 0;
	std::default_random_engine random(static_cast<int>(time(nullptr)));
	uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
	for (int i = 0; i < sampleNum; ++i) {
		for (int j = 0; j < sampleNum; ++j) {
			double u = uniformDistribution(random), v = uniformDistribution(random);
			double uradian = 2 * PI * u, vradian = 2 * PI * v;
			auto p1 = vec3(radius * sin(uradian) * cos(vradian),
			               radius * sin(uradian) * sin(vradian),
			               radius * cos(uradian));
			auto p2 = bs(controlPoints, u, v);
			error += static_cast<double>(glm::length(p1 - p2));
		}
	}
	error /= radius * sampleNum * sampleNum;

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

void processGUI() {
	ImGui::Begin("STATUS: ");
	string output("Average Error: ");
	output += to_string(error);
	ImGui::Text("Sampled from 100 points in each dimension");
	ImGui::Text(output.c_str());

	ImGui::Text("");
	ImGui::Text("Geometry Type: ");
	int og_type = g_type;
	ImGui::RadioButton("circle", &g_type, 0);
	ImGui::SameLine();
	ImGui::RadioButton("sphere", &g_type, 1);
	if (og_type != g_type) {
		change = true;
		if (g_type == 0) gType = GeometryType::Curve;
		else gType = GeometryType::Surface;
	}

	ImGui::Text("Parametrization method: ");
	int op_type = p_type;
	ImGui::RadioButton("uniform", &p_type, 0);
	ImGui::SameLine();
	ImGui::RadioButton("chordal", &p_type, 1);
	ImGui::SameLine();
	ImGui::RadioButton("centripetal", &p_type, 2);
	if (op_type != p_type) {
		change = true;
		if (p_type == 0) fitter.parametrizationType = Fitter::ParametrizationMethod::Uniform;
		else if (p_type == 1) fitter.parametrizationType = Fitter::ParametrizationMethod::Chordal;
		else fitter.parametrizationType = Fitter::ParametrizationMethod::Centripetal;
	}

	int odata_num = data_num;
	ImGui::SliderInt("dataNum", &data_num, 3, maxData);
	if (odata_num != data_num) {
		generateCircle(data_num);
		generateSphere(data_num);
		change = true;
	}

	if (g_type == 0) {
		int obc_type = bc_type;
		int obc_p = bc_p;
		ImGui::Text("b-spline type: ");
		ImGui::RadioButton("open", &bc_type, 0);
		ImGui::SameLine();
		ImGui::RadioButton("clamped", &bc_type, 1);
		ImGui::SameLine();
		ImGui::RadioButton("closed", &bc_type, 2);
		if (obc_type != bc_type) {
			change = true;
			if (bc_type == 0) bc = BSpline::BSplineType::Open;
			else if (bc_type == 1) bc = BSpline::BSplineType::Clamped;
			else {
				bc = BSpline::BSplineType::Closed;
				if (bc_p > data_num / 2) bc_p = min(data_num / 2, 20);

			}
		}
		bc_p = min(bc_p,
		           bc == BSpline::BSplineType::Closed
			           ? data_num / 2
			           : data_num);
		ImGui::SliderInt("P",
		                 &bc_p,
		                 1,
		                 min(bc_type == 2
			                     ? data_num / 2
			                     : data_num,
		                     20));
		if (obc_p != bc_p) {
			change = true;
		}
	} else {
		int obsu_type = bsu_type, obsv_type = bsv_type;
		ImGui::Text("uDir BSpline type: ");
		ImGui::RadioButton("uOpen", &bsu_type, 0);
		ImGui::SameLine();
		ImGui::RadioButton("uClamped", &bsu_type, 1);
		ImGui::SameLine();
		ImGui::RadioButton("uClosed", &bsu_type, 2);
		ImGui::Text("vDir bSpline type: ");
		ImGui::RadioButton("vOpen", &bsv_type, 0);
		ImGui::SameLine();
		ImGui::RadioButton("vClamped", &bsv_type, 1);
		ImGui::SameLine();
		ImGui::RadioButton("vClosed", &bsv_type, 2);
		if (obsu_type != bsu_type || obsv_type != bsv_type) {
			change = true;
			if (bsu_type == 0) bsu = BSpline::BSplineType::Open;
			else if (bsu_type == 1) bsu = BSpline::BSplineType::Clamped;
			else {
				bsu = BSpline::BSplineType::Closed;
				if (bs_p > data_num / 2) bs_p = min(data_num / 2, 20);
			}
			if (bsv_type == 0) bsv = BSpline::BSplineType::Open;
			else if (bsv_type == 1) bsv = BSpline::BSplineType::Clamped;
			else {
				bsv = BSpline::BSplineType::Closed;
				if (bs_q > data_num / 2) bs_q = min(data_num / 2, 20);
			}
		}
		bs_p = min(bs_p,
		           bsu == BSpline::BSplineType::Closed
			           ? data_num / 2
			           : data_num);
		bs_q = min(bs_q,
		           bsv == BSpline::BSplineType::Closed
			           ? data_num / 2
			           : data_num);
		int obs_p = bs_p, obs_q = bs_q;
		ImGui::SliderInt("uP",
		                 &bs_p,
		                 1,
		                 min(bsu_type == 2
			                     ? data_num / 2
			                     : data_num,
		                     20));
		ImGui::SliderInt("vQ",
		                 &bs_q,
		                 1,
		                 min(bsv_type == 2
			                     ? data_num / 2
			                     : data_num,
		                     20));
		if (obs_p != bs_p || obs_q != bs_q) {
			change = true;
		}
	}

	ImGui::End();
}

int main() {
	GLFWwindow* window = initOpenGL();

	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	(void)io;
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 130");

	Shader shader("shader.vert", "shader.frag");

	generateCircle(data_num);
	generateSphere(data_num);
	// auto surfaceDataPoints = generateCylinder(10);

	gType = GeometryType::Curve;

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		processInput(window);

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		if (!useCamera) {
			processGUI();
		}

		if (change) {
			change = false;
			if (gType == GeometryType::Curve) {
				bspCurve(bc_p, *curveDataPoints, bc);
			} else {
				bspSurface(bs_p, bs_q, *surfaceDataPoints, bsu, bsv);
			}
		}

		render(shader);
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
