#pragma once
#include <vector>
#include <glm/glm.hpp>

inline bool equal(const double a, const double b) {
	if (abs(a - b) < 1e-8) return true;
	return false;
}

class Fitter {
public:
	static double reparameterize(double t, double factor, double offset);
	static void reparameterize(std::vector<double>& ts, const std::vector<double>& us);
	static std::vector<double> generateUniformKnots(int p, std::vector<double>& ts);
	static std::vector<double> generateAverageKnots(int p, std::vector<double>& ts);
	static std::vector<double> uniformParameterize(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> chordalParameterize(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> centriPetalParameterize(const std::vector<glm::vec3>& dataPoints, double alpha = 0.5);
	static std::vector<glm::vec3> interpolateCurve(int p,
	                                               const std::vector<glm::vec3>& dataPoints,
	                                               std::vector<double>& ts,
	                                               const std::vector<double>& us);
	static std::vector<std::vector<glm::vec3>> interpolateSurface(int p,
	                                                              int q,
	                                                              const std::vector<std::vector<glm::vec3>>& dataPoints,
	                                                              std::vector<std::vector<double>>& parameters,
	                                                              const std::vector<double>& us,
	                                                              const std::vector<double>& vs);
	static std::vector<std::vector<double>> parameterize(const std::vector<std::vector<glm::vec3>>& dataPoints);
};
