#ifndef WASM_XUTILS_H
#define WASM_XUTILS_H

#include <vector>

void setexp(std::vector<double> &s, double ds1, double smax, int nn);

double atanc(double y, double x, double thold);

#endif  // WASM_XUTILS_H
