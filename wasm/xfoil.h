#ifndef WASM_XFOIL_H
#define WASM_XFOIL_H

#include <tuple>
#include <vector>

struct XFoilState;

std::pair<double, double> mrcl(XFoilState &ctx, double cls);
void comset(XFoilState &ctx);
void cpcalc(int n, const std::vector<double> &q, double qinf, double minf, std::vector<double> &cp);
std::tuple<double, double, double, double, double> clcalc(int n, const std::vector<double> &x, const std::vector<double> &y,
                                                          const std::vector<double> &gam, const std::vector<double> &gam_a,
                                                          double alfa, double minf, double qinf, double xref, double yref);
void cdcalc(XFoilState &ctx);
void tecalc(XFoilState &ctx);
void naca(XFoilState &ctx, int ides1);
void pangen(XFoilState &ctx, bool shopar);

#endif  // WASM_XFOIL_H
