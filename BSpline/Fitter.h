#pragma once
#include <vector>
#include <glm/glm.hpp>

#include "BSpline.h"

inline bool equal(const double a, const double b) {
	if (abs(a - b) < 1e-8) return true;
	return false;
}

inline double reparameterize(double t, double factor, double offset) {
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
	enum class ParametrizationMethod {
		Uniform, Chordal, Centripetal
	};

	enum class KnotGenerationMethod {
		Uniform, Average
	};

	// ParametrizationMethod
	static std::vector<double> uniformParametrization(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> chordalParametrization(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> centripetalParametrization(const std::vector<glm::vec3>& dataPoints);

	// KnotGenerationMethod
	std::vector<double> generateKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;
	std::vector<double> generateUniformKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;
	std::vector<double> generateAverageKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;

	// Parametrization
	std::vector<double> curveParametrization(const std::vector<glm::vec3>& dataPoints) const;
	std::vector<std::vector<double>>
	surfaceParametrization(const std::vector<std::vector<glm::vec3>>& dataPoints) const;

	// Interpolation
	static std::vector<glm::vec3> computeControlPoints(const BSpline& bsp,
	                                                   const std::vector<glm::vec3>& dataPoints,
	                                                   const std::vector<double>& ts);
	BSpline* interpolateCurve(int p,
	                          const std::vector<glm::vec3>& dataPoints,
	                          std::vector<glm::vec3>& controlPoints,
	                          BSpline::BSplineType type = BSpline::BSplineType::Closed) const;
	BSplineSurface* interpolateSurface(int p,
	                                   int q,
	                                   const std::vector<std::vector<glm::vec3>>& dataPoints,
	                                   std::vector<std::vector<glm::vec3>>& controlPoints,
	                                   BSpline::BSplineType utype = BSpline::BSplineType::Closed,
	                                   BSpline::BSplineType vtype = BSpline::BSplineType::Closed) const;

	ParametrizationMethod parametrizationType = ParametrizationMethod::Centripetal;
	KnotGenerationMethod knotGenerationType = KnotGenerationMethod::Uniform;
};
