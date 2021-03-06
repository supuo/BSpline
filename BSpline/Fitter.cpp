#include "Fitter.h"

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "BSpline.h"

using namespace std;
using glm::vec3;

vector<double> Fitter::uniformParametrization(const vector<vec3>& dataPoints) const {
	int size = static_cast<int>(dataPoints.size());
	vector<double> ret(size);
	double step = 1. / (size - 1);
	for (int i = 1; i < size; ++i) {
		ret[i] = ret[i - 1] + step;
	}
	return ret;
}

vector<double> Fitter::chordalParametrization(const std::vector<vec3>& dataPoints) const {
	int size = static_cast<int>(dataPoints.size());
	double accumulation = 0, total = 0;
	for (int i = 0; i < size - 1; ++i) {
		total += length(dataPoints[i + 1] - dataPoints[i]);
	}
	if (equal(total, 0)) {
		return uniformParametrization(dataPoints);
	}
	vector<double> ret(size);
	ret[0] = 0;
	for (int i = 1; i < size - 1; ++i) {
		accumulation += length(dataPoints[i + 1] - dataPoints[i]);
		ret[i] = accumulation / total;
	}
	ret[size - 1] = 1;
	return ret;
}

vector<double> Fitter::centripetalParametrization(const vector<vec3>& dataPoints) const {
	int size = static_cast<int>(dataPoints.size());
	double accumulation = 0, total = 0;
	for (int i = 0; i < size - 1; ++i) {
		total += sqrt(length(dataPoints[i + 1] - dataPoints[i]));
	}
	if (equal(total, 0)) {
		return uniformParametrization(dataPoints);
	}
	vector<double> ret(size);
	ret[0] = 0;
	for (int i = 1; i < size - 1; ++i) {
		accumulation += sqrt(length(dataPoints[i + 1] - dataPoints[i]));
		ret[i] = accumulation / total;
	}
	ret[size - 1] = 1;
	return ret;
}

std::vector<double> Fitter::generateKnots(int p, int n, std::vector<double>& ts, int bType) const {
	if (kType == 1) {
		return generateAverageKnots(p, n, ts, bType);
	}
	return generateUniformKnots(p, n, ts, bType);
}

vector<double> Fitter::generateUniformKnots(int p, int n, vector<double>& ts, int bType) const {
	int m = n + p + 1;
	vector<double> us(m + 1);
	if (bType == 1) {
		double du = 1. / (n - p + 1);
		for (int j = 1; j <= n - p; ++j) {
			us[j + p] = j * du;
		}
		fill(us.begin() + m - p, us.end(), 1);
		return us;
	}
	double du = 1. / m;
	for (int i = 0; i <= m; ++i) {
		us[i] = du * i;
	}
	reparameterize(ts, us);
	return us;
}

vector<double> Fitter::generateAverageKnots(int p, int n, vector<double>& ts, int bType) const {
	// todo: ��������open & closed
	if (bType != 1) return generateUniformKnots(p, n, ts, bType);
	int m = n + p + 1;
	vector<double> us(m + 1);
	for (int j = 1; j <= n - p; ++j) {
		for (int i = j; i <= j + p - 1; ++i) {
			us[j + p] += ts[i];
		}
		us[j + p] /= p;
	}
	for (int i = m - p; i <= m; ++i) {
		us[i] = 1;
	}
	return us;
}

vector<double> Fitter::curveParametrization(const std::vector<vec3>& dataPoints) const {
	if (pType == 0) {
		return uniformParametrization(dataPoints);
	}
	if (pType == 1) {
		return chordalParametrization(dataPoints);
	}
	return centripetalParametrization(dataPoints);
}

vector<vector<double>> Fitter::surfaceParametrization(const vector<vector<vec3>>& dataPoints) const {
	int m = static_cast<int>(dataPoints.size()) - 1;
	int n = static_cast<int>(dataPoints[0].size()) - 1;
	vector<vector<double>> uNet, vNet;
	vector<vector<double>> parameter(2);
	parameter[0].resize(m + 1);
	parameter[1].resize(n + 1);
	for (int j = 0; j <= n; ++j) {
		vector<vec3> columnData(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnData[i] = dataPoints[i][j];
		}
		uNet.emplace_back(curveParametrization(columnData));
	}
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= m; ++j) {
			parameter[0][j] += uNet[i][j] / (n + 1);
		}
	}
	for (int i = 0; i <= m; ++i) {
		vector<vec3> rowData(n + 1);
		for (int j = 0; j <= n; ++j) {
			rowData[j] = dataPoints[i][j];
		}
		vNet.emplace_back(curveParametrization(rowData));
	}
	for (int i = 0; i <= m; ++i) {
		for (int j = 0; j <= n; ++j) {
			parameter[1][j] += vNet[i][j] / (m + 1);
		}
	}
	return parameter;
}

std::vector<vec3> Fitter::interpolateControlPoints(const BSpline& bsp,
                                                   const vector<vec3>& dataPoints,
                                                   const std::vector<double>& ts) const {
	int p = bsp.p;
	int n = static_cast<int>(dataPoints.size()) - 1;
	vector<vec3> controlPoints;
	Eigen::MatrixXd N(n + 1, n + 1);
	for (int i = 0; i <= n; ++i) {
		auto coefficients = bsp.computeCoefficients(ts[i]);
		for (int j = 0; j <= n; ++j) {
			N(i, j) = coefficients[j];
		}
	}
	Eigen::MatrixXd D(n + 1, 3);
	for (int i = 0; i <= n; ++i) {
		D(i, 0) = static_cast<double>(dataPoints[i][0]);
		D(i, 1) = static_cast<double>(dataPoints[i][1]);
		D(i, 2) = static_cast<double>(dataPoints[i][2]);
	}
	if (bsp.bType == 2) {
		Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n - p + 1);
		for (int i = 0; i <= n - p; ++i) {
			M(i, i) = 1;
		}
		for (int i = 1; i <= p; ++i) {
			M(n - p + i, i - 1) = 1;
		}

		Eigen::MatrixXd P = (N * M).householderQr().solve(D);
		for (int i = 0; i <= n - p; ++i) {
			controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
		}
		for (int i = 0; i <= p - 1; ++i) {
			controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
		}
	} else {
		Eigen::MatrixXd P = N.householderQr().solve(D);
		for (int i = 0; i <= n; ++i) {
			controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
		}
	}
	return controlPoints;
}

BSpline* Fitter::interpolateCurve(int p,
                                  const vector<vec3>& dataPoints,
                                  vector<vec3>& controlPoints,
                                  int bType) const {
	std::vector<double> ts = curveParametrization(dataPoints);
	std::vector<double> us = generateKnots(p, static_cast<int>(ts.size()) - 1, ts, bType);
	BSpline* bsp = new BSpline(p, us, bType);
	controlPoints = interpolateControlPoints(*bsp, dataPoints, ts);
	return bsp;
}

BSplineSurface* Fitter::interpolateSurface(int p,
                                           int q,
                                           const std::vector<std::vector<vec3>>& dataPoints,
                                           std::vector<std::vector<vec3>>& controlPoints,
                                           int uType,
                                           int vType) const {
	vector<vector<double>> ts = surfaceParametrization(dataPoints);
	vector<double> us = generateKnots(p, static_cast<int>(ts[0].size()) - 1, ts[0], uType);
	vector<double> vs = generateKnots(q, static_cast<int>(ts[1].size()) - 1, ts[1], vType);
	BSplineSurface* bs = new BSplineSurface(p, q, us, vs, uType, vType);
	int m = static_cast<int>(dataPoints.size()) - 1, n = static_cast<int>(dataPoints[0].size()) - 1;
	vector<vector<vec3>> Q(m + 1, vector<vec3>(n + 1));
	for (int d = 0; d <= n; ++d) {
		vector<vec3> columnDataPoints(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnDataPoints[i] = dataPoints[i][d];
		}
		vector<vec3> intermediate = interpolateControlPoints(bs->ubsp, columnDataPoints, ts[0]);
		for (int i = 0; i <= m; ++i) {
			Q[i][d] = intermediate[i];
		}
	}
	controlPoints.resize(m + 1);
	for (int c = 0; c <= m; ++c) {
		controlPoints[c] = interpolateControlPoints(bs->vbsp, Q[c], ts[1]);
	}
	return bs;
}

std::vector<vec3> Fitter::approximateControlPoints(const BSpline& bsp,
                                                   int h,
                                                   const std::vector<glm::vec3>& dataPoints,
                                                   const std::vector<double>& ts) const {
	int n = static_cast<int>(dataPoints.size()) - 1;

	vector<vector<double>> N(n + 1);
	for (int k = 0; k <= n; ++k) {
		N[k] = bsp.computeCoefficients(ts[k]);
	}
	Eigen::MatrixXd matrixN(n - 1, h - 1);
	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < h - 1; ++j) {
			matrixN(i, j) = N[i + 1][j + 1];
		}
	}
	Eigen::MatrixXd D(n + 1, 3);
	for (int i = 0; i <= n; ++i) {
		D(i, 0) = static_cast<double>(dataPoints[i][0]);
		D(i, 1) = static_cast<double>(dataPoints[i][1]);
		D(i, 2) = static_cast<double>(dataPoints[i][2]);
	}

	Eigen::MatrixXd Qk(n - 1, 3);
	for (int i = 0; i < n - 1; ++i) {
		Qk.row(i) = D.row(i + 1) - N[i + 1][0] * D.row(0) - N[i + 1][h] * D.row(n);
	}
	Eigen::MatrixXd Q(h - 1, 3);
	for (int i = 0; i < h - 1; ++i) {
		Q = matrixN.transpose() * Qk;
	}
	Eigen::MatrixXd M = matrixN.transpose() * matrixN;
	Eigen::MatrixXd P = M.householderQr().solve(Q);
	vector<vec3> controlPoints(h + 1);
	controlPoints[0] = dataPoints[0];
	for (int i = 1; i < h; ++i) {
		auto p = P.row(i - 1);
		controlPoints[i] = vec3(p(0), p(1), p(2));
	}
	controlPoints[h] = dataPoints[n];
	return controlPoints;
}

BSpline* Fitter::approximateCurve(int p,
                                  int h,
                                  const std::vector<glm::vec3>& dataPoints,
                                  std::vector<glm::vec3>& controlPoints) const {
	std::vector<double> ts = curveParametrization(dataPoints);
	std::vector<double> us = generateKnots(p, h, ts, 1);
	BSpline* bsp = new BSpline(p, us, 1);
	controlPoints = approximateControlPoints(*bsp, h, dataPoints, ts);
	return bsp;
}

BSplineSurface* Fitter::approximateSurface(int p,
                                           int q,
                                           int e,
                                           int f,
                                           const std::vector<std::vector<glm::vec3>>& dataPoints,
                                           std::vector<std::vector<glm::vec3>>& controlPoints) const {
	vector<vector<double>> ts = surfaceParametrization(dataPoints);
	vector<double> us = generateKnots(p, e, ts[0], 1);
	vector<double> vs = generateKnots(q, f, ts[1], 1);
	BSplineSurface* bs = new BSplineSurface(p, q, us, vs, 1, 1);
	int m = static_cast<int>(dataPoints.size()) - 1, n = static_cast<int>(dataPoints[0].size()) - 1;
	vector<vector<vec3>> Q(e + 1, vector<vec3>(n + 1));
	for (int d = 0; d <= n; ++d) {
		vector<vec3> columnDataPoints(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnDataPoints[i] = dataPoints[i][d];
		}
		vector<vec3> intermediate = approximateControlPoints(bs->ubsp, e, columnDataPoints, ts[0]);
		for (int i = 0; i <= e; ++i) {
			Q[i][d] = intermediate[i];
		}
	}
	controlPoints.resize(e + 1);
	for (int c = 0; c <= e; ++c) {
		controlPoints[c] = approximateControlPoints(bs->vbsp, f, Q[c], ts[1]);
	}
	return bs;
}
