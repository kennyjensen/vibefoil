// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XUTILS_H
#define WASM_XUTILS_H

#include <vector>

void setexp(std::vector<double> &s, double ds1, double smax, int nn);

double atanc(double y, double x, double thold);

#endif  // WASM_XUTILS_H
