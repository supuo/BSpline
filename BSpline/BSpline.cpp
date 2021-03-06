#include "BSpline.h"

#include <vector>

#include "Fitter.h"

using namespace std;
using namespace glm;

BSpline::BSpline(int _p, const vector<double>& _us, int _type):
	p(_p),
	us(_us),
	bType(_type) {
	m = static_cast<int>(us.size()) - 1;
	n = m - p - 1;
}

vector<double> BSpline::computeCoefficients(double u) const {
	vector<double> N(n + 1);
	int k = 0;
	if (bType == 1) {
		if (equal(u, us[p])) {
			N[0] = 1;
			return N;
		}
		if (equal(u, us[n + 1])) {
			N[n] = 1;
			return N;
		}
	}
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

vec3 BSpline::operator()(const std::vector<vec3>& controlPoints, double u) const {
	u = reparameterize(u, us[n + 1] - us[p], us[p]);
	vec3 point(0);
	auto coefficients = computeCoefficients(u);
	int h = static_cast<int>(controlPoints.size()) - 1;
	for (int i = 0; i <= n && i <= h; i++) point += static_cast<float>(coefficients[i]) * controlPoints[i];
	return point;
}

BSplineSurface::BSplineSurface(int _p,
                               int _q,
                               const vector<double>& us,
                               const vector<double>& vs,
                               int utype,
                               int vtype):
	ubsp(_p, us, utype),
	vbsp(_q, vs, vtype) {}

vec3 BSplineSurface::operator()(const vector<vector<vec3>>& controlPoints, double u, double v) const {
	int m = static_cast<int>(controlPoints.size()) - 1, n = static_cast<int>(controlPoints[0].size()) - 1;
	vector<vec3> intermediate(m + 1);
	for (int i = 0; i <= m; ++i) {
		vector<vec3> rowControlPoints(n + 1);
		for (int j = 0; j <= n; ++j) {
			rowControlPoints[j] = controlPoints[i][j];
		}
		intermediate[i] = vbsp(rowControlPoints, v);
	}
	return ubsp(intermediate, u);
}
