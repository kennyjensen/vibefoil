#include "xutils.h"

#include <cmath>
#include <algorithm>
#include <iostream>

void setexp(std::vector<double> &s, double ds1, double smax, int nn) {
    const double sigma = smax / ds1;
    const int nex = nn - 1;
    const double rnex = static_cast<double>(nex);
    const double rni = 1.0 / rnex;

    const double aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    const double bbb = rnex * (rnex - 1.0) / 2.0;
    const double ccc = rnex - sigma;

    double disc = bbb * bbb - 4.0 * aaa * ccc;
    disc = std::max(0.0, disc);

    double ratio = 0.0;
    if (nex <= 1) {
        throw std::runtime_error("SETEXP: Cannot fill array.  N too small.");
    } else if (nex == 2) {
        ratio = -ccc / bbb + 1.0;
    } else {
        ratio = (-bbb + std::sqrt(disc)) / (2.0 * aaa) + 1.0;
    }

    if (ratio != 1.0) {
        for (int iter = 1; iter <= 100; ++iter) {
            const double sigman = (std::pow(ratio, nex) - 1.0) / (ratio - 1.0);
            const double res = std::pow(sigman, rni) - std::pow(sigma, rni);
            const double dresdr = rni * std::pow(sigman, rni) * (rnex * std::pow(ratio, nex - 1) - sigman) /
                                  (std::pow(ratio, nex) - 1.0);
            const double dratio = -res / dresdr;
            ratio = ratio + dratio;
            if (std::abs(dratio) < 1.0e-5) {
                break;
            }
            if (iter == 100) {
                std::cout << "SETEXP: Convergence failed.  Continuing anyway ..." << std::endl;
            }
        }
    }

    s[1] = 0.0;
    double ds = ds1;
    for (int n = 2; n <= nn; ++n) {
        s[n] = s[n - 1] + ds;
        ds = ds * ratio;
    }
}

double atanc(double y, double x, double thold) {
    const double pi = std::acos(-1.0);
    const double tpi = 2.0 * pi;

    const double thnew = std::atan2(y, x);
    const double dthet = thnew - thold;
    const double dtcorr = dthet - tpi * static_cast<int>((dthet + std::copysign(pi, dthet)) / tpi);
    return thold + dtcorr;
}
