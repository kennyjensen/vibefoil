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
