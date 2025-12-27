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
