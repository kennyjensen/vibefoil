// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XSOLVE_H
#define WASM_XSOLVE_H

#include <vector>

struct XFoilState;

void gauss(int nsiz, int nn, std::vector<std::vector<double>> &z, std::vector<std::vector<double>> &r, int nrhs);
void gauss(int nsiz, int nn, std::vector<std::vector<double>> &z, std::vector<double> &r, int nrhs);
void ludcmp(int nsiz, int n, std::vector<std::vector<double>> &a, std::vector<int> &indx);
void baksub(int nsiz, int n, std::vector<std::vector<double>> &a, const std::vector<int> &indx, std::vector<double> &b);

void blsolv(XFoilState &ctx);

#endif  // WASM_XSOLVE_H
