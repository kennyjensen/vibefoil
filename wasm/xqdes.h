// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XQDES_H
#define WASM_XQDES_H

struct XFoilState;

void QDES(XFoilState &ctx);
void NEWPLOTQ(XFoilState &ctx);
void QPLINI(XFoilState &ctx, bool ldef);
void QSPLOT(XFoilState &ctx);
void QSPPLT(XFoilState &ctx, int iqspl1, int iqspl2, int kqsp, int nt);
void IQSGET(XFoilState &ctx);
void SPLQSP(XFoilState &ctx, int kqsp);
void SMOOQ(XFoilState &ctx, int kq1, int kq2, int kqsp);
double QINCOM(double qc, double qinf, double tklam);
void GAMQSP(XFoilState &ctx, int kqsp);
void SYMQSP(XFoilState &ctx, int kqsp);
void MIXED(XFoilState &ctx, int kqsp, int niterq);
void GAMLIN(XFoilState &ctx, int i, int j, double coef);

#endif  // WASM_XQDES_H
