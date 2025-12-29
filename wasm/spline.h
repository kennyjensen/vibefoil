// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_SPLINE_H
#define WASM_SPLINE_H

#include <vector>

typedef std::vector<double> DVec;

typedef std::vector<std::vector<double>> DVec2;

void trisol(DVec &a, DVec &b, DVec &c, DVec &d, int kk);

void splind(const DVec &x, DVec &xs, const DVec &s, int n, double xs1, double xs2);

double seval(double ss, const DVec &x, const DVec &xs, const DVec &s, int n);

double deval(double ss, const DVec &x, const DVec &xs, const DVec &s, int n);

double sinvrt(double si, double xi, const DVec &x, const DVec &xs, const DVec &s, int n);

void scalc(const DVec &x, const DVec &y, DVec &s, int n);

void segspl(const DVec &x, DVec &xs, const DVec &s, int n);

double d2val(double ss, const DVec &x, const DVec &xs, const DVec &s, int n);

double curv(double ss, const DVec &x, const DVec &xs, const DVec &y, const DVec &ys, const DVec &s, int n);

#endif  // WASM_SPLINE_H
