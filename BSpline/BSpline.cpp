#include "BSpline.h"

#include <vector>

#include "Fitter.h"

using namespace std;
using namespace glm;

BSpline::BSpline(int _p, const std::vector<vec3>& dataPoints):
	p(_p) {
	std::vector<double> ts = Fitter::centriPetalParameterize(dataPoints);
	us = Fitter::generateUniformKnots(p, ts);
	m = static_cast<int>(us.size()) - 1;
	n = m - p - 1;
	controlPoints = Fitter::interpolateCurve(p, dataPoints, ts, us);
}

vector<double> BSpline::computeCoefficients(int n, double u, const vector<double>& us) {
	vector<double> N(n + 1);
	int m = static_cast<int>(us.size()) - 1;
	int p = m - n - 1;
	int k = 0;
	// if (equal(u, us[p])) {
	// 	N[0] = 1;
	// 	return N;
	// }
	// if (equal(u, us[n + 1])) {
	// 	N[n] = 1;
	// 	return N;
	// }
	for (; k < m; ++k) {
		if (us[k] <= u && u < us[k + 1]) {
			break;
		}
	}
	if (k > n) {
		k = n;
	}
	N[k] = 1;
	for (int d = 1; d <= p; ++d) {
		N[k - d] = (us[k + 1] - u) / (us[k + 1] - us[k - d + 1]) * N[(k - d) + 1];
		for (int i = k - d + 1; i <= k - 1; ++i) {
			N[i] = (u - us[i]) / (us[i + d] - us[i]) * N[i] + (us[i + d + 1] - u) / (
				us[i + d + 1] - us[i + 1]) * N[i + 1];
		}
		N[k] = (u - us[k]) / (us[k + d] - us[k]) * N[k];
	}
	return N;
}

vec3 BSpline::deBoor(const vector<vec3>& controlPoints, const vector<double>& us, int n, int p, double u) {
	u = Fitter::reparameterize(u, us[n + 1] - us[p], us[p]);
	vec3 point(0);
	auto coefficients = computeCoefficients(n, u, us);
	for (int i = 0; i <= n; i++) point += static_cast<float>(coefficients[i]) * controlPoints[i];
	return point;
}

vec3 BSpline::operator()(double u) const {
	return deBoor(controlPoints, us, n, p, u);
}

void BSpline::wrap() {
	for (int i = 1; i <= p; ++i) {
		controlPoints[n - p + i] = controlPoints[i - 1];
	}
}

BSplineSurface::BSplineSurface(int _p, int _q, const vector<vector<vec3>>& dataPoints):
	p(_p),
	q(_q) {
	vector<vector<double>> parameters = Fitter::parameterize(dataPoints);
	us = Fitter::generateUniformKnots(p, parameters[0]);
	vs = Fitter::generateUniformKnots(q, parameters[1]);
	controlPoints = Fitter::interpolateSurface(p, q, dataPoints, parameters, us, vs);
}

vec3 BSplineSurface::operator()(double u, double v) {
	vector<vec3> intermediate;
	int m = static_cast<int>(controlPoints.size()) - 1, n = static_cast<int>(controlPoints[0].size()) - 1;
	for (int i = 0; i <= m; ++i) {
		vector<vec3> rowControlPoints;
		for (int j = 0; j <= n; ++j) {
			rowControlPoints.emplace_back(controlPoints[i][j]);
		}
		intermediate.emplace_back(BSpline::deBoor(rowControlPoints, vs, n, q, v));
	}
	return BSpline::deBoor(intermediate, us, m, p, u);
}
