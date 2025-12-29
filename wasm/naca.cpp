// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#include "naca.h"

#include <cmath>
#include <iostream>

std::pair<int, std::string> naca4(int ides, std::vector<double> &xx, std::vector<double> &yt, std::vector<double> &yc,
                                 int nside, std::vector<double> &xb, std::vector<double> &yb) {
    const std::string digits = "0123456789";
    const double an = 1.5;

    const int n4 = ides / 1000;
    const int n3 = (ides - n4 * 1000) / 100;
    const int n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    const int n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;

    const double m = static_cast<double>(n4) / 100.0;
    const double p = static_cast<double>(n3) / 10.0;
    const double t = static_cast<double>(n2 * 10 + n1) / 100.0;

    const double anp = an + 1.0;
    for (int i = 1; i <= nside; ++i) {
        const double frac = static_cast<double>(i - 1) / static_cast<double>(nside - 1);
        if (i == nside) {
            xx[i] = 1.0;
        } else {
            xx[i] = 1.0 - anp * frac * std::pow(1.0 - frac, an) - std::pow(1.0 - frac, anp);
        }

        yt[i] = (0.29690 * std::sqrt(xx[i]) - 0.12600 * xx[i] - 0.35160 * xx[i] * xx[i] + 0.28430 * xx[i] * xx[i] * xx[i]
                - 0.10150 * xx[i] * xx[i] * xx[i] * xx[i])
                * t / 0.20;

        if (xx[i] < p) {
            yc[i] = m / (p * p) * (2.0 * p * xx[i] - xx[i] * xx[i]);
        } else {
            yc[i] = m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * xx[i] - xx[i] * xx[i]);
        }
    }

    int ib = 0;
    for (int i = nside; i >= 1; --i) {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
    }
    for (int i = 2; i <= nside; ++i) {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
    }

    std::string name = "NACA";
    name.append(10 - static_cast<int>(name.size()), ' ');
    name = name.substr(0, 5) + digits.substr(static_cast<size_t>(n4), 1) + digits.substr(static_cast<size_t>(n3), 1)
           + digits.substr(static_cast<size_t>(n2), 1) + digits.substr(static_cast<size_t>(n1), 1);
    return {ib, name};
}

std::pair<int, std::string> naca5(int ides, std::vector<double> &xx, std::vector<double> &yt, std::vector<double> &yc,
                                 int nside, std::vector<double> &xb, std::vector<double> &yb) {
    const std::string digits = "0123456789";
    const double an = 1.5;

    const int n5 = ides / 10000;
    const int n4 = (ides - n5 * 10000) / 1000;
    const int n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    const int n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    const int n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;

    const int n543 = 100 * n5 + 10 * n4 + n3;

    double m = 0.0;
    double c = 0.0;
    if (n543 == 210) {
        m = 0.0580;
        c = 361.4;
    } else if (n543 == 220) {
        m = 0.1260;
        c = 51.64;
    } else if (n543 == 230) {
        m = 0.2025;
        c = 15.957;
    } else if (n543 == 240) {
        m = 0.2900;
        c = 6.643;
    } else if (n543 == 250) {
        m = 0.3910;
        c = 3.230;
    } else {
        std::cout << "Illegal 5-digit designation" << std::endl;
        std::cout << "First three digits must be 210, 220, ... 250" << std::endl;
        return {0, ""};
    }

    const double t = static_cast<double>(n2 * 10 + n1) / 100.0;

    const double anp = an + 1.0;
    for (int i = 1; i <= nside; ++i) {
        const double frac = static_cast<double>(i - 1) / static_cast<double>(nside - 1);
        if (i == nside) {
            xx[i] = 1.0;
        } else {
            xx[i] = 1.0 - anp * frac * std::pow(1.0 - frac, an) - std::pow(1.0 - frac, anp);
        }

        yt[i] = (0.29690 * std::sqrt(xx[i]) - 0.12600 * xx[i] - 0.35160 * xx[i] * xx[i] + 0.28430 * xx[i] * xx[i] * xx[i]
                - 0.10150 * xx[i] * xx[i] * xx[i] * xx[i])
                * t / 0.20;

        if (xx[i] < m) {
            yc[i] = (c / 6.0) * (xx[i] * xx[i] * xx[i] - 3.0 * m * xx[i] * xx[i] + m * m * (3.0 - m) * xx[i]);
        } else {
            yc[i] = (c / 6.0) * (m * m * m) * (1.0 - xx[i]);
        }
    }

    int ib = 0;
    for (int i = nside; i >= 1; --i) {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
    }
    for (int i = 2; i <= nside; ++i) {
        ib += 1;
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
    }

    std::string name = "NACA";
    name.append(10 - static_cast<int>(name.size()), ' ');
    name = name.substr(0, 5) + digits.substr(static_cast<size_t>(n5), 1) + digits.substr(static_cast<size_t>(n4), 1)
           + digits.substr(static_cast<size_t>(n3), 1) + digits.substr(static_cast<size_t>(n2), 1)
           + digits.substr(static_cast<size_t>(n1), 1);
    return {ib, name};
}
