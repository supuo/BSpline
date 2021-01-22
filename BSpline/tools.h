#pragma once

namespace tools {
	bool equal(const double a, const double b) {
		if (abs(a - b) < 1e-8) return true;
		return false;
	}
}
