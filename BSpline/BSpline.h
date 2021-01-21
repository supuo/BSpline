#pragma once
#include <vector>
#include <glm/glm.hpp>

struct BSpline {
	BSpline(int _p, const std::vector<double>& _knots) :
		p(_p),
		m(static_cast<int>(_knots.size()) - 1),
		n(m - p - 1),
		knots(_knots) {}

	std::vector<double> computeCoefficients(double u);
	glm::vec3 operator()(double u);
	void wrap();
	
	int p = 0, m = 0, n = 0;
	std::vector<glm::vec3> controlPoints;
	std::vector<double> knots;
};

class BSplineSurface {
	// public:
	// 	BSplineSurface(const std::vector<std::vector<Eigen::Vector3d>>& dataPoints, int p, int q);
	// 	Eigen::Vector3d calculatePoint(double u, double v);
	// 	std::vector<std::vector<Eigen::Vector3d>> generateSurface(double step = 0.01);
	// 	static void interpolateSurface(const std::vector<std::vector<Eigen::Vector3d>>& dataPoints);
	// private:
	// 	BSpline ubs, vbs;
	// 	int p = 0; // u向阶数
	// 	int q = 0; // v向阶数
	// 	std::vector<std::vector<Eigen::Vector3d>> controlPoints;
	// 	std::vector<double> knotsU;
	// 	std::vector<double> knotsV;
};
