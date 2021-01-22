#pragma once
#include <utility>
#include <vector>
#include <glm/glm.hpp>

#include "Fitter.h"

struct BSpline {
	BSpline(int _p, const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> computeCoefficients(int n, double u, const std::vector<double>& us);
	static glm::vec3 deBoor(const std::vector<glm::vec3>& controlPoints,
	                        const std::vector<double>& us,
	                        int n,
	                        int p,
	                        double u);
	glm::vec3 operator()(double u) const;
	void wrap();

	int p = 0, m = 0, n = 0;
	std::vector<double> us;
	std::vector<glm::vec3> controlPoints;
};

struct BSplineSurface {
	BSplineSurface(int _p, int _q, const std::vector<std::vector<glm::vec3>>& dataPoints);
	glm::vec3 operator()(double u, double v);

	int p = 0, q = 0;
	std::vector<double> us, vs;
	std::vector<std::vector<glm::vec3>> controlPoints;
};
