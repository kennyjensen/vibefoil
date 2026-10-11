#include <iomanip>
#include <iostream>
#include "../../wasm/naca.h"
#include "../../wasm/spline.h"
#include "../../wasm/xgeom.h"
#include "../../wasm/xfoil.h"
#include "../../wasm/xpanel.h"
// The correlation is file-private; include its implementation in this test TU.
#include "../../wasm/xblsys.cpp"

int main() {
    std::cout << std::scientific << std::setprecision(17);
    int mode;
    std::cin >> mode;
    if (mode == 1) {
        double hk, rt, msq, factor;
        std::cin >> hk >> rt >> msq >> factor;
        const auto [cf, cf_hk, cf_rt, cf_msq] = cft(hk, rt, msq, factor);
        std::cout << cf << ' ' << cf_hk << ' ' << cf_rt << ' ' << cf_msq << '\n';
    } else if (mode == 2) {
        int ides, n;
        char sharp;
        std::cin >> ides >> n >> sharp;
        DVec xx(n + 1), yt(n + 1), yc(n + 1), x(2 * n + 1), y(2 * n + 1);
        auto result = ides <= 9999 ? naca4(ides, xx, yt, yc, n, x, y, sharp == 'T')
                                  : naca5(ides, xx, yt, yc, n, x, y, sharp == 'T');
        for (int i = 1; i <= result.first; ++i) std::cout << x[i] << ' ' << y[i] << '\n';
    } else if (mode == 3) {
        int n;
        double first, last;
        std::cin >> n >> first >> last;
        DVec x(n + 1), s(n + 1), xp(n + 1);
        for (int i = 1; i <= n; ++i) std::cin >> s[i] >> x[i];
        segspld(x, xp, s, n, first, last);
        for (int i = 1; i <= n; ++i) std::cout << xp[i] << '\n';
    } else if (mode == 4) {
        int n;
        std::cin >> n;
        DVec x(n + 1), y(n + 1), s(n + 1), xp(n + 1), yp(n + 1);
        for (int i = 1; i <= n; ++i) std::cin >> x[i] >> y[i];
        scalc(x, y, s, n);
        segspl(x, xp, s, n);
        segspl(y, yp, s, n);
        std::cout << lefind(x, xp, y, yp, s, n) << '\n';
    } else if (mode == 5) {
        XFoilState ctx;
        int ip, sharp, image;
        double xi, yi, nxi, nyi;
        std::cin >> ctx.N >> ip >> sharp >> image >> ctx.ALFA >> ctx.YIMAGE >> xi >> yi >> nxi >> nyi;
        ctx.LIMAGE = image != 0;
        ctx.CHORD = 1.0;
        for (int j = 1; j <= ctx.N; ++j) {
            std::cin >> ctx.X[j] >> ctx.Y[j];
            ctx.GAM[j] = 0.1 * (j - 3);
            ctx.GAMU[j][1] = 0.2 * j;
            ctx.GAMU[j][2] = -0.1 * j;
            ctx.SIG[j] = 0.005 * j;
            ctx.QF0[j] = 0.03 * j;
            ctx.QF1[j] = -0.02 * j;
            ctx.QF2[j] = 0.01 * j;
            ctx.QF3[j] = -0.04 * j;
        }
        scalc(ctx.X, ctx.Y, ctx.S, ctx.N);
        segspl(ctx.X, ctx.XP, ctx.S, ctx.N);
        segspl(ctx.Y, ctx.YP, ctx.S, ctx.N);
        ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY);
        tecalc(ctx);
        ctx.SHARP = sharp != 0;
        apcalc(ctx);
        auto result = psilin(ctx, ip, xi, yi, nxi, nyi, true, true);
        std::cout << result.psi << ' ' << result.psi_ni << ' ' << result.qt1 << ' ' << result.qt2 << ' '
                  << ctx.Z_QINF << ' ' << ctx.Z_ALFA << ' ' << ctx.Z_QDOF0 << ' ' << ctx.Z_QDOF1 << ' '
                  << ctx.Z_QDOF2 << ' ' << ctx.Z_QDOF3 << '\n';
        for (int j = 1; j <= ctx.N; ++j)
            std::cout << ctx.DZDG[j] << ' ' << ctx.DZDN[j] << ' ' << ctx.DQDG[j] << ' '
                      << ctx.DZDM[j] << ' ' << ctx.DQDM[j] << '\n';
    } else if (mode == 6) {
        XFoilState ctx;
        double cl;
        std::cin >> ctx.MATYP >> ctx.RETYP >> ctx.MINF1 >> ctx.REINF1 >> cl;
        mrcl(ctx, cl);
        comset(ctx);
        std::cout << "REFERENCE_RESULT\n" << ctx.MATYP << ' ' << ctx.RETYP << ' '
                  << ctx.MINF << ' ' << ctx.REINF << ' ' << ctx.MINF_CL << ' ' << ctx.REINF_CL << ' '
                  << ctx.TKLAM << ' ' << ctx.TKL_MSQ << ' ' << ctx.CPSTAR << ' ' << ctx.QSTAR << '\n';
    } else if (mode == 7) {
        XFoilState ctx;
        std::cin >> ctx.N;
        for (int j = 1; j <= ctx.N; ++j) std::cin >> ctx.X[j] >> ctx.Y[j];
        scalc(ctx.X, ctx.Y, ctx.S, ctx.N);
        segspl(ctx.X, ctx.XP, ctx.S, ctx.N);
        segspl(ctx.Y, ctx.YP, ctx.S, ctx.N);
        const double sle = lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N);
        const double dx = .5 * (ctx.X[1] + ctx.X[ctx.N]) - seval(sle, ctx.X, ctx.XP, ctx.S, ctx.N);
        const double dy = .5 * (ctx.Y[1] + ctx.Y[ctx.N]) - seval(sle, ctx.Y, ctx.YP, ctx.S, ctx.N);
        ctx.CHORD = std::sqrt(dx * dx + dy * dy);
        tecalc(ctx);
        std::cout << ctx.ANTE << ' ' << ctx.ASTE << ' ' << ctx.DSTE << ' ' << ctx.SHARP << '\n';
    }
}
