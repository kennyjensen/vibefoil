// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XOPER_H
#define WASM_XOPER_H

struct XFoilState;
struct XBlState;

void mhinge(XFoilState &ctx);
void viscal(XFoilState &ctx, XBlState &bl, int niter1);

#endif  // WASM_XOPER_H
