#include "Fitter.h"

#include <vector>
#include <Eigen/Dense>
#include "BSpline.h"

using namespace std;
using glm::vec3;

vector<double> Fitter::centriPetalParameterize(const vector<vec3>& dataPoints, const double alpha) {
	int size = static_cast<int>(dataPoints.size());
	vector<double> ret(size);
	double accumulation = 0, total = 0;
	for (int i = 0; i < size - 1; ++i) {
		total += pow(length(dataPoints[i + 1] - dataPoints[i]), alpha);
	}
	ret[0] = 0;
	for (int i = 1; i < size - 1; ++i) {
		accumulation += pow(length(dataPoints[i + 1] - dataPoints[i]), alpha);
		ret[i] = accumulation / total;
	}
	ret[size - 1] = 1;
	return ret;
}

vector<double> Fitter::computeKnots(int p, const vector<double>& ts) {
	// int n = static_cast<int>(ts.size()) - 1;
	// int m = n + p + 1;
	// vector<double> us(m + 1);
	// for (int j = 1; j <= n - p; ++j) {
	// 	for (int i = j; i <= j + p - 1; ++i) {
	// 		us[j + p] += ts[i];
	// 	}
	// 	us[j + p] /= p;
	// }
	// for (int i = m - p; i <= m; ++i) {
	// 	us[i] = 1;
	// }
	// return us;

	int n = static_cast<int>(ts.size()) - 1;
	int m = n + p + 1;
	vector<double> us(m + 1);
	double du = 1. / m;
	for (int i = 0; i <= m; ++i) {
		us[i] = du * i;
	}
	return us;
}

double Fitter::reparameterize(double t, double factor, double offset) {
	return offset + factor * t;
}

void Fitter::reparameterize(vector<double>& ts, const vector<double>& us) {
	int n = static_cast<int>(ts.size()) - 1, m = static_cast<int>(us.size()) - 1, p = m - n - 1;
	for (auto& t : ts) {
		t = reparameterize(t, us[n + 1] - us[p], us[p]);
	}
}

BSpline* Fitter::interpolateCurve(int p, const vector<vec3>& dataPoints) {
	vector<double> ts = centriPetalParameterize(dataPoints);
	vector<double> us = computeKnots(p, ts);
	reparameterize(ts, us);
	BSpline* bsp = new BSpline(p, us);
	int n = static_cast<int>(dataPoints.size()) - 1;
	vector<vec3> controlPoints;
	Eigen::MatrixXf N(n + 1, n + 1);
	for (int i = 0; i <= n; ++i) {
		auto coefficients = bsp->computeCoefficients(ts[i]);
		for (int j = 0; j <= n; ++j) {
			N(i, j) = static_cast<float>(coefficients[j]);
		}
	}
	Eigen::MatrixXf D(n + 1, 3), P(n + 1, 3);
	for (int i = 0; i <= n; ++i) {
		D(i, 0) = dataPoints[i][0];
		D(i, 1) = dataPoints[i][1];
		D(i, 2) = dataPoints[i][2];
	}
	P = N.colPivHouseholderQr().solve(D);
	for (int i = 0; i <= n; ++i) {
		bsp->controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
	}
	return bsp;
}

vector<vector<double>> Fitter::parameterize(const vector<vector<vec3>>& dataPoints) {
	// int m = static_cast<int>(dataPoints.size()) - 1;
	// int n = static_cast<int>(dataPoints[0].size()) - 1;
	// vector<vector<double>> uNet, vNet;
	// vector<vector<double>> parameter(2);
	// parameter[0].resize(m + 1);
	// parameter[1].resize(n + 1);
	// for (int j = 0; j <= n; ++j) {
	// 	vector<double> u(m + 1);
	// 	vector<double> w(m + 1);
	// 	for (int i = 0; i <= m; ++i) {
	// 		u[i] = dataPoints[i][j][1];
	// 		w[i] = dataPoints[i][j][2];
	// 	}
	// 	uNet.emplace_back(centriPetalParameterize(u, w));
	// }
	// for (int i = 0; i <= n; ++i) {
	// 	for (int j = 0; j <= m; ++j) {
	// 		parameter[0][j] += uNet[i][j] / (n + 1);
	// 	}
	// }
	// for (int i = 0; i <= m; ++i) {
	// 	vector<double> v(n + 1);
	// 	vector<double> w(n + 1);
	// 	for (int j = 0; j <= n; ++j) {
	// 		v[j] = dataPoints[i][j][0];
	// 		w[j] = dataPoints[i][j][2];
	// 	}
	// 	vNet.emplace_back(centriPetalParameterize(v, w));
	// }
	// for (int i = 0; i <= m; ++i) {
	// 	for (int j = 0; j <= n; ++j) {
	// 		parameter[1][j] += vNet[i][j] / (m + 1);
	// 	}
	// }
	// return parameter;
	return vector<vector<double>>();
}
