#ifndef WASM_NACA_H
#define WASM_NACA_H

#include <string>
#include <utility>
#include <vector>

std::pair<int, std::string> naca4(int ides, std::vector<double> &xx, std::vector<double> &yt, std::vector<double> &yc,
                                 int nside, std::vector<double> &xb, std::vector<double> &yb);

std::pair<int, std::string> naca5(int ides, std::vector<double> &xx, std::vector<double> &yt, std::vector<double> &yc,
                                 int nside, std::vector<double> &xb, std::vector<double> &yb);

#endif  // WASM_NACA_H
