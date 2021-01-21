#pragma once
#include <vector>
#include <glm/glm.hpp>

#include "BSpline.h"

inline bool equal(const double a, const double b) {
	if (abs(a - b) < 1e-8) return true;
	return false;
}

class Fitter {
public:
	static double reparameterize(double t, double factor, double offset);
	static void reparameterize(std::vector<double>& ts, const std::vector<double>& us);
	static std::vector<double> computeKnots(int p, const std::vector<double>& ts);
	static std::vector<double> centriPetalParameterize(const std::vector<glm::vec3>& dataPoints, double alpha = 0.5);
	static BSpline* interpolateCurve(int p, const std::vector<glm::vec3>& dataPoints);
	static std::vector<std::vector<double>> parameterize(const std::vector<std::vector<glm::vec3>>& dataPoints);
};
