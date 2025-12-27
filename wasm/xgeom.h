#ifndef WASM_XGEOM_H
#define WASM_XGEOM_H

#include <vector>

struct GeoparResults {
    double sle;
    double chord;
    double area;
    double radle;
    double angte;
    double ei11a;
    double ei22a;
    double apx1a;
    double apx2a;
    double ei11t;
    double ei22t;
    double apx1t;
    double apx2t;
    double thick;
    double cambr;
};

double lefind(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
              const std::vector<double> &yp, const std::vector<double> &s, int n);

double sopps(double si, const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
             const std::vector<double> &yp, const std::vector<double> &s, int n, double sle);

void aecalc(int n, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &t, int itype,
            double &area, double &xcen, double &ycen, double &ei11, double &ei22, double &apx1, double &apx2);

void tccalc(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
            const std::vector<double> &yp, const std::vector<double> &s, int n, double &thick, double &xthick,
            double &cambr, double &xcambr);

GeoparResults geopar(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &y,
                     const std::vector<double> &yp, const std::vector<double> &s, int n, std::vector<double> &t);

#endif  // WASM_XGEOM_H
