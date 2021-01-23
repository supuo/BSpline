#pragma once

#include <vector>
#include <glm/glm.hpp>

struct BSpline {
	enum class BSplineType {
		Open, Clamped, Closed
	};

	BSpline(int _p, const std::vector<double>& us, BSplineType type);
	std::vector<double> computeCoefficients(double u) const;
	glm::vec3 operator()(const std::vector<glm::vec3>& controlPoints, double u) const;

	int n = 0, m = 0, p = 0;
	std::vector<double> us;
	BSplineType type;
};

struct BSplineSurface {
	BSplineSurface(int _p,
	               int _q,
	               const std::vector<double>& us,
	               const std::vector<double>& vs,
	               BSpline::BSplineType utype,
	               BSpline::BSplineType vtype);
	glm::vec3 operator()(const std::vector<std::vector<glm::vec3>>& controlPoints, double u, double v) const;

	BSpline ubsp, vbsp;
};
