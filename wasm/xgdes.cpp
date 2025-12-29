// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#include "xgdes.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "spline.h"

std::tuple<double, double, double, double> getxyf(const std::vector<double> &x, const std::vector<double> &xp,
                                                 const std::vector<double> &y, const std::vector<double> &yp,
                                                 const std::vector<double> &s, int n, double tops, double bots,
                                                 double xf, double yf) {
    if (xf == -999.0) {
        throw std::runtime_error("GETXYF: hinge x not specified.");
    }

    tops = s[1] + (x[1] - xf);
    bots = s[n] - (x[n] - xf);
    tops = sinvrt(tops, xf, x, xp, s, n);
    bots = sinvrt(bots, xf, x, xp, s, n);
    const double topy = seval(tops, y, yp, s, n);
    const double boty = seval(bots, y, yp, s, n);

    std::cout << std::endl;
    std::cout << "  Top    surface:  y =" << std::setw(8) << std::fixed << std::setprecision(4) << topy << "     y/t = 1.0"
              << std::endl;
    std::cout << "  Bottom surface:  y =" << std::setw(8) << std::fixed << std::setprecision(4) << boty << "     y/t = 0.0"
              << std::endl;

    if (yf == -999.0) {
        throw std::runtime_error("GETXYF: hinge y not specified.");
    }

    if (yf == 999.0) {
        throw std::runtime_error("GETXYF: y/t hinge input not supported.");
    }

    return {tops, bots, xf, yf};
}
