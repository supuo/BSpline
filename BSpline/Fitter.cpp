#include "Fitter.h"

#include <vector>
#include <Eigen/Dense>

#include "BSpline.h"

using namespace std;
using glm::vec3;

vector<double> Fitter::uniformParameterize(const vector<vec3>& dataPoints) {
	int size = static_cast<int>(dataPoints.size());
	vector<double> ret(size);
	double step = 1. / (size - 1);
	for (int i = 1; i < size; ++i) {
		ret[i] = ret[i - 1] + step;
	}
	return ret;
}

vector<double> Fitter::chordalParameterize(const std::vector<glm::vec3>& dataPoints) {
	int size = static_cast<int>(dataPoints.size());
	double accumulation = 0, total = 0;
	for (int i = 0; i < size - 1; ++i) {
		total += length(dataPoints[i + 1] - dataPoints[i]);
	}
	if (equal(total, 0)) {
		return uniformParameterize(dataPoints);
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

vector<double> Fitter::centriPetalParameterize(const vector<vec3>& dataPoints, const double alpha) {
	int size = static_cast<int>(dataPoints.size());
	double accumulation = 0, total = 0;
	for (int i = 0; i < size - 1; ++i) {
		total += pow(length(dataPoints[i + 1] - dataPoints[i]), alpha);
	}
	if (equal(total, 0)) {
		return uniformParameterize(dataPoints);
	}
	vector<double> ret(size);
	ret[0] = 0;
	for (int i = 1; i < size - 1; ++i) {
		accumulation += pow(length(dataPoints[i + 1] - dataPoints[i]), alpha);
		ret[i] = accumulation / total;
	}
	ret[size - 1] = 1;
	return ret;
}

vector<double> Fitter::generateUniformKnots(int p, vector<double>& ts) {
	int n = static_cast<int>(ts.size()) - 1;
	int m = n + p + 1;
	vector<double> us(m + 1);
	double du = 1. / m;
	for (int i = 0; i <= m; ++i) {
		us[i] = du * i;
	}
	reparameterize(ts, us);
	return us;
}

vector<double> Fitter::generateAverageKnots(int p, vector<double>& ts) {
	int n = static_cast<int>(ts.size()) - 1;
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
	reparameterize(ts, us);
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

vector<vec3> Fitter::interpolateCurve(int p,
                                      const vector<vec3>& dataPoints,
                                      vector<double>& ts,
                                      const vector<double>& us) {
	int n = static_cast<int>(dataPoints.size()) - 1;
	vector<vec3> controlPoints;
	Eigen::MatrixXd N(n + 1, n + 1);
	for (int i = 0; i <= n; ++i) {
		auto coefficients = BSpline::computeCoefficients(n, ts[i], us);
		for (int j = 0; j <= n; ++j) {
			N(i, j) = (coefficients[j]);
		}
	}
	Eigen::MatrixXd D(n + 1, 3);
	for (int i = 0; i <= n; ++i) {
		D(i, 0) = dataPoints[i][0];
		D(i, 1) = dataPoints[i][1];
		D(i, 2) = dataPoints[i][2];
	}
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n - p + 1, 3);
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n - p + 1);
	for (int i = 0; i <= n - p; ++i) {
		M(i, i) = 1;
	}
	for (int i = 1; i <= p; ++i) {
		M(n - p + i, i - 1) = 1;
	}
	
	P = (N * M).colPivHouseholderQr().solve(D);
	
	for (int i = 0; i <= n - p; ++i) {
		controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
	}
	for (int i = 0; i <= p - 1; ++i) {
		controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
	}
	//
	// Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n + 1, 3);
	// P = N.colPivHouseholderQr().solve(D);
	// for (int i = 0; i <= n; ++i) {
	// 	controlPoints.emplace_back(P(i, 0), P(i, 1), P(i, 2));
	// }
	//
	return controlPoints;
}

std::vector<std::vector<vec3>> Fitter::interpolateSurface(int p,
                                                          int q,
                                                          const std::vector<std::vector<vec3>>& dataPoints,
                                                          std::vector<std::vector<double>>& parameters,
                                                          const std::vector<double>& us,
                                                          const std::vector<double>& vs) {
	vector<vector<vec3>> controlPoints;
	int m = dataPoints.size() - 1, n = dataPoints[0].size() - 1;
	vector<vector<vec3>> Q(m + 1, vector<vec3>(n + 1));
	for (int d = 0; d <= n; ++d) {
		vector<vec3> columnDataPoints(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnDataPoints[i] = dataPoints[i][d];
		}
		vector<vec3> intermediate = interpolateCurve(p, columnDataPoints, parameters[0], us);
		for (int i = 0; i <= m; ++i) {
			Q[i][d] = intermediate[i];
		}
	}
	controlPoints.resize(m + 1);
	for (int c = 0; c <= m; ++c) {
		controlPoints[c] = interpolateCurve(q, Q[c], parameters[1], vs);
	}
	return controlPoints;
}

vector<vector<double>> Fitter::parameterize(const vector<vector<vec3>>& dataPoints) {
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
		uNet.emplace_back(centriPetalParameterize(columnData));
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
		vNet.emplace_back(centriPetalParameterize(rowData));
	}
	for (int i = 0; i <= m; ++i) {
		for (int j = 0; j <= n; ++j) {
			parameter[1][j] += vNet[i][j] / (m + 1);
		}
	}
	return parameter;
}
