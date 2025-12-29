// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XGDES_H
#define WASM_XGDES_H

#include <tuple>
#include <vector>

std::tuple<double, double, double, double> getxyf(const std::vector<double> &x, const std::vector<double> &xp,
                                                 const std::vector<double> &y, const std::vector<double> &yp,
                                                 const std::vector<double> &s, int n, double tops, double bots,
                                                 double xf, double yf);

#endif  // WASM_XGDES_H
