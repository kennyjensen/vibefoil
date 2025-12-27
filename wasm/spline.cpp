#include "spline.h"

#include <cmath>
#include <iostream>

void trisol(DVec &a, DVec &b, DVec &c, DVec &d, int kk) {
    for (int k = 2; k <= kk; ++k) {
        const int km = k - 1;
        c[km] = c[km] / a[km];
        d[km] = d[km] / a[km];
        a[k] = a[k] - b[k] * c[km];
        d[k] = d[k] - b[k] * d[km];
    }

    d[kk] = d[kk] / a[kk];

    for (int k = kk - 1; k >= 1; --k) {
        d[k] = d[k] - c[k] * d[k + 1];
    }
}

void splind(const DVec &x, DVec &xs, const DVec &s, int n, double xs1, double xs2) {
    const int nmax = 1000;
    if (n > nmax) {
        throw std::runtime_error("SPLIND: array overflow, increase NMAX");
    }

    DVec a(static_cast<size_t>(n + 1), 0.0);
    DVec b(static_cast<size_t>(n + 1), 0.0);
    DVec c(static_cast<size_t>(n + 1), 0.0);

    for (int i = 2; i <= n - 1; ++i) {
        const double dsm = s[i] - s[i - 1];
        const double dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }

    if (xs1 == 999.0) {
        a[1] = 2.0;
        c[1] = 1.0;
        xs[1] = 3.0 * (x[2] - x[1]) / (s[2] - s[1]);
    } else if (xs1 == -999.0) {
        a[1] = 1.0;
        c[1] = 1.0;
        xs[1] = 2.0 * (x[2] - x[1]) / (s[2] - s[1]);
    } else {
        a[1] = 1.0;
        c[1] = 0.0;
        xs[1] = xs1;
    }

    if (xs2 == 999.0) {
        b[n] = 1.0;
        a[n] = 2.0;
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    } else if (xs2 == -999.0) {
        b[n] = 1.0;
        a[n] = 1.0;
        xs[n] = 2.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    } else {
        a[n] = 1.0;
        b[n] = 0.0;
        xs[n] = xs2;
    }

    if (n == 2 && xs1 == -999.0 && xs2 == -999.0) {
        b[n] = 1.0;
        a[n] = 2.0;
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1]);
    }

    trisol(a, b, c, xs, n);
}

static int spline_locate(double ss, const DVec &s, int n) {
    int ilow = 1;
    int i = n;

    while (i - ilow > 1) {
        const int imid = (i + ilow) / 2;
        if (ss < s[imid]) {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    return i;
}

double seval(double ss, const DVec &x, const DVec &xs, const DVec &s, int n) {
    const int i = spline_locate(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    return t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
}

double deval(double ss, const DVec &x, const DVec &xs, const DVec &s, int n) {
    const int i = spline_locate(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    const double deval_val = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    return deval_val / ds;
}

double sinvrt(double si, double xi, const DVec &x, const DVec &xs, const DVec &s, int n) {
    const double sisav = si;
    for (int iter = 0; iter < 10; ++iter) {
        const double res = seval(si, x, xs, s, n) - xi;
        const double resp = deval(si, x, xs, s, n);
        const double ds = -res / resp;
        si = si + ds;
        if (std::abs(ds / (s[n] - s[1])) < 1.0e-5) {
            return si;
        }
    }
    std::cout << std::endl;
    std::cout << "SINVRT: spline inversion failed. Input value returned." << std::endl;
    return sisav;
}

void scalc(const DVec &x, const DVec &y, DVec &s, int n) {
    s[1] = 0.0;
    for (int i = 2; i <= n; ++i) {
        s[i] = s[i - 1] + std::sqrt((x[i] - x[i - 1]) * (x[i] - x[i - 1]) + (y[i] - y[i - 1]) * (y[i] - y[i - 1]));
    }
}

void segspl(const DVec &x, DVec &xs, const DVec &s, int n) {
    if (s[1] == s[2]) {
        throw std::runtime_error("SEGSPL:  First input point duplicated");
    }
    if (s[n] == s[n - 1]) {
        throw std::runtime_error("SEGSPL:  Last  input point duplicated");
    }

    int iseg0 = 1;
    for (int iseg = 2; iseg <= n - 2; ++iseg) {
        if (s[iseg] == s[iseg + 1]) {
            const int nseg = iseg - iseg0 + 1;
            DVec xt(static_cast<size_t>(nseg + 1), 0.0);
            DVec xst(static_cast<size_t>(nseg + 1), 0.0);
            DVec st(static_cast<size_t>(nseg + 1), 0.0);
            for (int i = 1; i <= nseg; ++i) {
                xt[i] = x[iseg0 + i - 1];
                xst[i] = xs[iseg0 + i - 1];
                st[i] = s[iseg0 + i - 1];
            }
            splind(xt, xst, st, nseg, -999.0, -999.0);
            for (int i = 1; i <= nseg; ++i) {
                xs[iseg0 + i - 1] = xst[i];
            }
            iseg0 = iseg + 1;
        }
    }

    const int nseg = n - iseg0 + 1;
    DVec xt(static_cast<size_t>(nseg + 1), 0.0);
    DVec xst(static_cast<size_t>(nseg + 1), 0.0);
    DVec st(static_cast<size_t>(nseg + 1), 0.0);
    for (int i = 1; i <= nseg; ++i) {
        xt[i] = x[iseg0 + i - 1];
        xst[i] = xs[iseg0 + i - 1];
        st[i] = s[iseg0 + i - 1];
    }
    splind(xt, xst, st, nseg, -999.0, -999.0);
    for (int i = 1; i <= nseg; ++i) {
        xs[iseg0 + i - 1] = xst[i];
    }
}

double d2val(double ss, const DVec &x, const DVec &xs, const DVec &s, int n) {
    const int i = spline_locate(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    const double d2val_val = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    return d2val_val / (ds * ds);
}

double curv(double ss, const DVec &x, const DVec &xs, const DVec &y, const DVec &ys, const DVec &s, int n) {
    const int i = spline_locate(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;

    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    const double xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    const double xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

    const double cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    const double cy2 = ds * ys[i] - y[i] + y[i - 1];
    const double yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    const double ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

    double sd = std::sqrt(xd * xd + yd * yd);
    sd = std::max(sd, 0.001 * ds);

    return (xd * ydd - yd * xdd) / (sd * sd * sd);
}
