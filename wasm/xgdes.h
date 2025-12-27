#ifndef WASM_XGDES_H
#define WASM_XGDES_H

#include <tuple>
#include <vector>

std::tuple<double, double, double, double> getxyf(const std::vector<double> &x, const std::vector<double> &xp,
                                                 const std::vector<double> &y, const std::vector<double> &yp,
                                                 const std::vector<double> &s, int n, double tops, double bots,
                                                 double xf, double yf);

#endif  // WASM_XGDES_H
