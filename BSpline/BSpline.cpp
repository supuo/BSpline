#include "BSpline.h"

#include <vector>

#include "Fitter.h"

using namespace std;

std::vector<double> BSpline::computeCoefficients(double u) {
	vector<double> N(n + 1);
	// if (equal(u, knots[p])) {
	// 	N[0] = 1;
	// 	return N;
	// }
	// if (equal(u, knots[n + 1])) {
	// 	N[n] = 1;
	// 	return N;
	// }
	int k = 0;
	for (; k < m; ++k) {
		if (knots[k] <= u && u < knots[k + 1]) {
			break;
		}
	}
	if (k > n) {
		k = n;
	}
	N[k] = 1;
	for (int d = 1; d <= p; ++d) {
		N[k - d] = (knots[k + 1] - u) / (knots[k + 1] - knots[k - d + 1]) * N[(k - d) + 1];
		for (int i = k - d + 1; i <= k - 1; ++i) {
			N[i] = (u - knots[i]) / (knots[i + d] - knots[i]) * N[i] + (knots[i + d + 1] - u) / (
				knots[i + d + 1] - knots[i + 1]) * N[i + 1];
		}
		N[k] = (u - knots[k]) / (knots[k + d] - knots[k]) * N[k];
	}
	return N;
}

glm::vec3 BSpline::operator()(double u) {
	u = Fitter::reparameterize(u, knots[n+1] - knots[p], knots[p]);
	glm::vec3 point(0);
	int n = static_cast<int>(controlPoints.size()) - 1;
	auto coefficients = computeCoefficients(u);
	for (int i = 0; i <= n; i++) point += static_cast<float>(coefficients[i]) * controlPoints[i];
	return point;
}

void BSpline::wrap() {
	for (int i = 1; i <= p; ++i) {
		controlPoints[n - p + i] = controlPoints[i - 1];
	}
}

// BSplineSurface::BSplineSurface(const vector<vector<Vector3d>>& dataPoints, int p, int q) : p(p),
// 	q(q) {
// 	interpolateSurface(dataPoints);
// }
//
// void BSplineSurface::interpolateSurface(const vector<vector<Vector3d>>& dataPoints) {
// 	vector<vector<double>> parameters = Fitter::parameterize(dataPoints);
// 	//p[0] is u-direction s, p[1] is v-direction t;
// 	knotsU = Fitter::computeKnots(parameters[0], p);
// 	knotsV = Fitter::computeKnots(parameters[1], q);
// 	int m = static_cast<int>(parameters[0].size()) - 1, n = static_cast<int>(parameters[1].size()) - 1;
// 	vector<vector<Vector3d>> Q(m + 1, vector<Vector3d>(n + 1));
// 	for (int d = 0; d <= n; ++d) {
// 		vector<Vector3d> columnDataPoints(m + 1);
// 		for (int i = 0; i <= m; ++i) {
// 			columnDataPoints[i] = dataPoints[i][d];
// 		}
// 		vector<Vector3d> intermediate = Fitter::interpolateCurve(columnDataPoints, p, parameters[0], knotsU);
// 		for (int i = 0; i <= m; ++i) {
// 			Q[i][d] = intermediate[i];
// 		}
// 	}
// 	controlPoints.resize(m + 1);
// 	for (int c = 0; c <= m; ++c) {
// 		controlPoints[c] = Fitter::interpolateCurve(Q[c], q, parameters[1], knotsV);
// 	}
// }
//
// vector<vector<Vector3d>> BSplineSurface::generateSurface(const double step) {
// 	int nu = static_cast<int>(controlPoints.size()) - 1, nv = static_cast<int>(controlPoints[0].size()) - 1;
// 	int m = static_cast<int>((knotsU[nu + 1] - knotsU[p]) / step);
// 	int n = static_cast<int>((knotsV[nv + 1] - knotsV[q]) / step);
// 	vector<vector<Vector3d>> vertices(m + 1, vector<Vector3d>(n + 1));
// 	for (int i = 0; i <= m; ++i) {
// 		for (int j = 0; j <= n; ++j) {
// 			double u, v;
// 			if (i == m) {
// 				u = knotsU[nu + 1];
// 				v = knotsV[q] + j * step;
// 			} else if (j == n) {
// 				u = knotsU[p] + i * step;
// 				v = knotsV[nv + 1];
// 			} else {
// 				u = knotsU[p] + i * step;
// 				v = knotsV[q] + j * step;
// 			}
// 			vertices[i][j] = calculatePoint(u, v);
// 		}
// 	}
// 	return vertices;
// }
//
// Vector3d BSplineSurface::calculatePoint(const double u, const double v) {
// 	vector<Vector3d> intermediate;
// 	int mcp = static_cast<int>(controlPoints.size()), ncp = static_cast<int>(controlPoints[0].size());
// 	for (int i = 0; i < mcp; ++i) {
// 		vector<Vector3d> rowControlPoints;
// 		for (int j = 0; j < ncp; ++j) {
// 			rowControlPoints.emplace_back(controlPoints[i][j]);
// 		}
// 		intermediate.emplace_back(BSpline::DeBoor(rowControlPoints, knotsV, q, v));
// 	}
// 	return BSpline::DeBoor(intermediate, knotsU, p, u);
// }
