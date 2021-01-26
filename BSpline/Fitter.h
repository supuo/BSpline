#pragma once
#include <vector>
#include <glm/glm.hpp>

#include "BSpline.h"

inline bool equal(const double a, const double b) {
	if (abs(a - b) < 1e-8) return true;
	return false;
}

inline double reparameterize(const double t, const double factor, const double offset) {
	return offset + factor * t;
}

inline void reparameterize(std::vector<double>& ts, const std::vector<double>& us) {
	int n = static_cast<int>(ts.size()) - 1, m = static_cast<int>(us.size()) - 1, p = m - n - 1;
	for (auto& t : ts) {
		t = reparameterize(t, us[n + 1] - us[p], us[p]);
	}
}

class Fitter {
public:
	Fitter(int _pType = 0, int _kType = 0):
		pType(_pType),
		kType(_kType) {}

	// ParametrizationMethod
	std::vector<double> uniformParametrization(const std::vector<glm::vec3>& dataPoints) const;
	std::vector<double> chordalParametrization(const std::vector<glm::vec3>& dataPoints) const;
	std::vector<double> centripetalParametrization(const std::vector<glm::vec3>& dataPoints) const;

	// KnotGenerationMethod
	std::vector<double> generateKnots(int p, int n, std::vector<double>& ts, int bType) const;
	std::vector<double> generateUniformKnots(int p, int n, std::vector<double>& ts, int bType) const;
	std::vector<double> generateAverageKnots(int p, int n, std::vector<double>& ts, int bType) const;

	// Parametrization
	std::vector<double> curveParametrization(const std::vector<glm::vec3>& dataPoints) const;
	std::vector<std::vector<double>>
	surfaceParametrization(const std::vector<std::vector<glm::vec3>>& dataPoints) const;

	// Interpolation
	std::vector<glm::vec3> interpolateControlPoints(const BSpline& bsp,
	                                                const std::vector<glm::vec3>& dataPoints,
	                                                const std::vector<double>& ts) const;
	BSpline interpolateCurve(int p,
	                         const std::vector<glm::vec3>& dataPoints,
	                         std::vector<glm::vec3>& controlPoints,
	                         int bType) const;
	BSplineSurface interpolateSurface(int p,
	                                  int q,
	                                  const std::vector<std::vector<glm::vec3>>& dataPoints,
	                                  std::vector<std::vector<glm::vec3>>& controlPoints,
	                                  int uType,
	                                  int vType) const;

	// approximation
	std::vector<glm::vec3> approximateControlPoints(const BSpline& bsp,
	                                                int h,
	                                                const std::vector<glm::vec3>& dataPoints,
	                                                const std::vector<double>& ts) const;
	BSpline approximateCurve(int p,
	                         int h,
	                         const std::vector<glm::vec3>& dataPoints,
	                         std::vector<glm::vec3>& controlPoints) const;
	BSplineSurface approximateSurface(int p,
	                                  int q,
	                                  int e,
	                                  int f,
	                                  const std::vector<std::vector<glm::vec3>>& dataPoints,
	                                  std::vector<std::vector<glm::vec3>>& controlPoints) const;

	int pType; // 0 = uniform, 1 = chordal, 2 = centripetal;
	int kType; // 0 = uniform, 1 = average;
};
