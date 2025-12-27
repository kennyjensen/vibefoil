#ifndef WASM_XPANEL_H
#define WASM_XPANEL_H

#include <utility>
#include <vector>

struct XFoilState;

void apcalc(XFoilState &ctx);
void ncalc(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, int n,
           std::vector<double> &xn, std::vector<double> &yn);
struct PsinResult {
    double psi;
    double psi_ni;
    double qt1;
    double qt2;
    double qtanm;
};

PsinResult psilin(XFoilState &ctx, int i, double xi, double yi, double nxi, double nyi, bool geolin, bool siglin);
std::pair<double, double> pswlin(XFoilState &ctx, int i, double xi, double yi, double nxi, double nyi);
void ggcalc(XFoilState &ctx);
void qwcalc(XFoilState &ctx);
void qdcalc(XFoilState &ctx);
void xywake(XFoilState &ctx);
void stfind(XFoilState &ctx);
void iblpan(XFoilState &ctx);
void xicalc(XFoilState &ctx);
void uicalc(XFoilState &ctx);
void qvfue(XFoilState &ctx);
void qiset(XFoilState &ctx);
void gamqv(XFoilState &ctx);
void stmove(XFoilState &ctx);
void ueset(XFoilState &ctx);
void dsset(XFoilState &ctx);

#endif  // WASM_XPANEL_H
