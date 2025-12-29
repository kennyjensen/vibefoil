// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#include "xgeom.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "spline.h"
#include "xutils.h"

double lefind(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
              const std::vector<double> &yp, const std::vector<double> &s, int n) {
    const double dseps = (s[n] - s[1]) * 1.0e-5;

    const double xte = 0.5 * (x[1] + x[n]);
    const double yte = 0.5 * (y[1] + y[n]);

    int i = 3;
    for (i = 3; i <= n - 2; ++i) {
        const double dxte = x[i] - xte;
        const double dyte = y[i] - yte;
        const double dx = x[i + 1] - x[i];
        const double dy = y[i + 1] - y[i];
        const double dotp = dxte * dx + dyte * dy;
        if (dotp < 0.0) {
            break;
        }
    }

    double sle = s[i];

    if (s[i] == s[i - 1]) {
        return sle;
    }

    for (int iter = 1; iter <= 50; ++iter) {
        const double xle = seval(sle, x, xp, s, n);
        const double yle = seval(sle, y, yp, s, n);
        const double dxds = deval(sle, x, xp, s, n);
        const double dyds = deval(sle, y, yp, s, n);
        const double dxdd = d2val(sle, x, xp, s, n);
        const double dydd = d2val(sle, y, yp, s, n);

        const double xchord = xle - xte;
        const double ychord = yle - yte;

        const double res = xchord * dxds + ychord * dyds;
        const double ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        double dsle = -res / ress;

        dsle = std::max(dsle, -0.02 * std::abs(xchord + ychord));
        dsle = std::min(dsle, 0.02 * std::abs(xchord + ychord));
        sle = sle + dsle;
        if (std::abs(dsle) < dseps) {
            return sle;
        }
    }

    std::cout << "LEFIND:  LE point not found.  Continuing..." << std::endl;
    return s[i];
}

double sopps(double si, const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
             const std::vector<double> &yp, const std::vector<double> &s, int n, double sle) {
    const double slen = s[n] - s[1];

    const double xle = seval(sle, x, xp, s, n);
    const double yle = seval(sle, y, yp, s, n);
    const double xte = 0.5 * (x[1] + x[n]);
    const double yte = 0.5 * (y[1] + y[n]);
    const double chord = std::sqrt((xte - xle) * (xte - xle) + (yte - yle) * (yte - yle));
    const double dxc = (xte - xle) / chord;
    const double dyc = (yte - yle) / chord;

    int inp = 1;
    int inopp = n;
    if (si >= sle) {
        inp = n;
        inopp = 1;
    }
    const double sfrac = (si - sle) / (s[inp] - sle);
    double sopp = sle + sfrac * (s[inopp] - sle);

    if (std::abs(sfrac) <= 1.0e-5) {
        return sle;
    }

    const double xi = seval(si, x, xp, s, n);
    const double yi = seval(si, y, yp, s, n);
    const double xbar = (xi - xle) * dxc + (yi - yle) * dyc;

    for (int iter = 1; iter <= 12; ++iter) {
        const double xopp = seval(sopp, x, xp, s, n);
        const double yopp = seval(sopp, y, yp, s, n);
        const double xoppd = deval(sopp, x, xp, s, n);
        const double yoppd = deval(sopp, y, yp, s, n);

        const double res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar;
        const double resd = xoppd * dxc + yoppd * dyc;

        if (std::abs(res) / slen < 1.0e-5) {
            return sopp;
        }
        if (resd == 0.0) {
            break;
        }

        const double dsopp = -res / resd;
        sopp = sopp + dsopp;

        if (std::abs(dsopp) / slen < 1.0e-5) {
            return sopp;
        }
    }

    std::cout << std::endl;
    std::cout << "SOPPS: Opposite-point location failed. Continuing..." << std::endl;
    return sle + sfrac * (s[inopp] - sle);
}

void aecalc(int n, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &t, int itype,
            double &area, double &xcen, double &ycen, double &ei11, double &ei22, double &apx1, double &apx2) {
    const double pi = std::acos(-1.0);

    double sint = 0.0;
    double aint = 0.0;
    double xint = 0.0;
    double yint = 0.0;
    double xxint = 0.0;
    double xyint = 0.0;
    double yyint = 0.0;

    for (int io = 1; io <= n; ++io) {
        const int ip = (io == n) ? 1 : io + 1;

        const double dx = x[io] - x[ip];
        const double dy = y[io] - y[ip];
        const double xa = (x[io] + x[ip]) * 0.50;
        const double ya = (y[io] + y[ip]) * 0.50;
        const double ta = (t[io] + t[ip]) * 0.50;

        const double ds = std::sqrt(dx * dx + dy * dy);
        sint = sint + ds;

        if (itype == 1) {
            const double da = ya * dx;
            aint = aint + da;
            xint = xint + xa * da;
            yint = yint + ya * da / 2.0;
            xxint = xxint + xa * xa * da;
            xyint = xyint + xa * ya * da / 2.0;
            yyint = yyint + ya * ya * da / 3.0;
        } else {
            const double da = ta * ds;
            aint = aint + da;
            xint = xint + xa * da;
            yint = yint + ya * da;
            xxint = xxint + xa * xa * da;
            xyint = xyint + xa * ya * da;
            yyint = yyint + ya * ya * da;
        }
    }

    area = aint;

    if (aint == 0.0) {
        xcen = 0.0;
        ycen = 0.0;
        ei11 = 0.0;
        ei22 = 0.0;
        apx1 = 0.0;
        apx2 = std::atan2(1.0, 0.0);
        return;
    }

    xcen = xint / aint;
    ycen = yint / aint;

    const double eixx = yyint - ycen * ycen * aint;
    const double eixy = xyint - xcen * ycen * aint;
    const double eiyy = xxint - xcen * xcen * aint;

    const double eisq = 0.25 * (eixx - eiyy) * (eixx - eiyy) + eixy * eixy;
    const double sgn = std::copysign(1.0, eiyy - eixx);
    ei11 = 0.5 * (eixx + eiyy) - sgn * std::sqrt(eisq);
    ei22 = 0.5 * (eixx + eiyy) + sgn * std::sqrt(eisq);

    if (ei11 == 0.0 || ei22 == 0.0) {
        apx1 = 0.0;
        apx2 = std::atan2(1.0, 0.0);
    } else if (eisq / (ei11 * ei22) < std::pow(0.001 * sint, 4)) {
        apx1 = 0.0;
        apx2 = std::atan2(1.0, 0.0);
    } else {
        const double c1 = eixy;
        const double s1 = eixx - ei11;

        const double c2 = eixy;
        const double s2 = eixx - ei22;

        if (std::abs(s1) > std::abs(s2)) {
            apx1 = std::atan2(s1, c1);
            apx2 = apx1 + 0.5 * pi;
        } else {
            apx2 = std::atan2(s2, c2);
            apx1 = apx2 - 0.5 * pi;
        }

        if (apx1 < -0.5 * pi) {
            apx1 = apx1 + pi;
        }
        if (apx1 > 0.5 * pi) {
            apx1 = apx1 - pi;
        }
        if (apx2 < -0.5 * pi) {
            apx2 = apx2 + pi;
        }
        if (apx2 > 0.5 * pi) {
            apx2 = apx2 - pi;
        }
    }
}

void tccalc(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
            const std::vector<double> &yp, const std::vector<double> &s, int n, double &thick, double &xthick,
            double &cambr, double &xcambr) {
    const double sle = lefind(x, xp, y, yp, s, n);
    const double xle = seval(sle, x, xp, s, n);
    const double yle = seval(sle, y, yp, s, n);
    const double xte = 0.5 * (x[1] + x[n]);
    const double yte = 0.5 * (y[1] + y[n]);
    const double chord = std::sqrt((xte - xle) * (xte - xle) + (yte - yle) * (yte - yle));

    const double dxc = (xte - xle) / chord;
    const double dyc = (yte - yle) / chord;

    thick = 0.0;
    xthick = 0.0;
    cambr = 0.0;
    xcambr = 0.0;

    for (int i = 1; i <= n; ++i) {
        const double xbar = (x[i] - xle) * dxc + (y[i] - yle) * dyc;
        const double ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc;

        const double sopp = sopps(s[i], x, xp, y, yp, s, n, sle);
        const double xopp = seval(sopp, x, xp, s, n);
        const double yopp = seval(sopp, y, yp, s, n);

        const double ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc;

        const double yc = 0.5 * (ybar + ybarop);
        const double yt = std::abs(ybar - ybarop);

        if (std::abs(yc) > std::abs(cambr)) {
            cambr = yc;
            xcambr = xopp;
        }
        if (std::abs(yt) > std::abs(thick)) {
            thick = yt;
            xthick = xopp;
        }
    }
}

GeoparResults geopar(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
                     const std::vector<double> &yp, const std::vector<double> &s, int n, std::vector<double> &t) {
    const double sle = lefind(x, xp, y, yp, s, n);

    const double xle = seval(sle, x, xp, s, n);
    const double yle = seval(sle, y, yp, s, n);
    const double xte = 0.5 * (x[1] + x[n]);
    const double yte = 0.5 * (y[1] + y[n]);

    const double chsq = (xte - xle) * (xte - xle) + (yte - yle) * (yte - yle);
    const double chord = std::sqrt(chsq);

    const double curvle = curv(sle, x, xp, y, yp, s, n);

    double radle = 0.0;
    if (std::abs(curvle) > 0.001 * (s[n] - s[1])) {
        radle = 1.0 / curvle;
    }

    const double ang1 = std::atan2(-yp[1], -xp[1]);
    const double ang2 = atanc(yp[n], xp[n], ang1);
    const double angte = ang2 - ang1;

    for (int i = 1; i <= n; ++i) {
        t[i] = 1.0;
    }

    double area = 0.0;
    double xcena = 0.0;
    double ycena = 0.0;
    double ei11a = 0.0;
    double ei22a = 0.0;
    double apx1a = 0.0;
    double apx2a = 0.0;
    aecalc(n, x, y, t, 1, area, xcena, ycena, ei11a, ei22a, apx1a, apx2a);

    double slen = 0.0;
    double xcent = 0.0;
    double ycent = 0.0;
    double ei11t = 0.0;
    double ei22t = 0.0;
    double apx1t = 0.0;
    double apx2t = 0.0;
    aecalc(n, x, y, t, 2, slen, xcent, ycent, ei11t, ei22t, apx1t, apx2t);

    double thick = 0.0;
    double xthick = 0.0;
    double cambr = 0.0;
    double xcambr = 0.0;
    tccalc(x, xp, y, yp, s, n, thick, xthick, cambr, xcambr);

    std::cout << " Max thickness = " << std::setw(12) << std::fixed << std::setprecision(6) << thick << "  at x = "
              << std::setw(7) << std::setprecision(3) << xthick << std::endl;
    std::cout << " Max camber    = " << std::setw(12) << std::fixed << std::setprecision(6) << cambr << "  at x = "
              << std::setw(7) << std::setprecision(3) << xcambr << std::endl;

    GeoparResults out{};
    out.sle = sle;
    out.chord = chord;
    out.area = area;
    out.radle = radle;
    out.angte = angte;
    out.ei11a = ei11a;
    out.ei22a = ei22a;
    out.apx1a = apx1a;
    out.apx2a = apx2a;
    out.ei11t = ei11t;
    out.ei22t = ei22t;
    out.apx1t = apx1t;
    out.apx2t = apx2t;
    out.thick = thick;
    out.cambr = cambr;
    return out;
}
