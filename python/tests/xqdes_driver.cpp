#include <cmath>
#include <iomanip>
#include <iostream>
#include "../../wasm/xbl.h"
#include "../../wasm/xqdes.h"
#include "../../wasm/xgeom.h"
#include "../../wasm/spline.h"

int main() {
    XFoilState ctx;
    int iterations, lcp, image, sharp;
    std::cin >> ctx.N >> ctx.NSP >> ctx.IQ1 >> ctx.IQ2 >> iterations;
    std::cin >> ctx.ALFA >> ctx.MINF >> ctx.QINF >> ctx.XCMREF >> ctx.YCMREF
             >> ctx.PSIO >> ctx.QDOF0 >> ctx.QDOF1 >> ctx.QDOF2 >> ctx.QDOF3;
    std::cin >> lcp >> image >> sharp;
    ctx.LCPXX = lcp;
    ctx.LIMAGE = image;
    ctx.SHARP = sharp;
    std::cin >> ctx.ANTE >> ctx.ASTE >> ctx.DSTE >> ctx.XTE >> ctx.YTE >> ctx.YIMAGE;
    for (auto *field : {&ctx.X, &ctx.Y, &ctx.S, &ctx.XP, &ctx.YP, &ctx.NX, &ctx.NY, &ctx.APANEL})
        for (int i = 1; i <= ctx.N; ++i) std::cin >> (*field)[i];
    for (int i = 1; i <= ctx.NSP; ++i) std::cin >> ctx.SSPEC[i];
    for (int i = 1; i <= ctx.NSP; ++i) std::cin >> ctx.QSPEC[i][1];
    for (int i = 1; i <= ctx.N; ++i) std::cin >> ctx.GAM[i];
    for (int i = 1; i <= ctx.N; ++i) std::cin >> ctx.SIG[i];
    ctx.SLE = lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N);
    ctx.XLE = seval(ctx.SLE, ctx.X, ctx.XP, ctx.S, ctx.N);
    ctx.YLE = seval(ctx.SLE, ctx.Y, ctx.YP, ctx.S, ctx.N);
    const double dx = ctx.XTE - ctx.XLE, dy = ctx.YTE - ctx.YLE;
    ctx.CHORD = std::sqrt(dx * dx + dy * dy);
    MIXED(ctx, 1, iterations);
    std::cout << "REFERENCE_RESULT\n" << std::scientific << std::setprecision(17);
    for (int i = 1; i <= ctx.N; ++i)
        std::cout << ctx.X[i] << ' ' << ctx.Y[i] << ' ' << ctx.GAM[i] << '\n';
    std::cout << ctx.PSIO << ' ' << ctx.QDOF0 << ' ' << ctx.QDOF1 << ' ' << ctx.QDOF2 << ' '
              << ctx.QDOF3 << ' ' << ctx.CL << ' ' << ctx.CM << ' ' << ctx.CDP << ' '
              << ctx.CL_ALF << ' ' << ctx.CL_MSQ << '\n';
}
