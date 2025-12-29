// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XBLSYS_H
#define WASM_XBLSYS_H

#include <tuple>

struct XBlState;

void trchek(XBlState &bl);
void blsys(XBlState &bl);
void tesys(XBlState &bl, double cte, double tte, double dte);
void blprv(XBlState &bl, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei);
void blkin(XBlState &bl);
void blvar(XBlState &bl, int ityp);
void blmid(XBlState &bl, int ityp);

std::tuple<double, double, double> hkin(double h, double msq);

#endif  // WASM_XBLSYS_H
