#include "xsolve.h"

#include <cmath>
#include <stdexcept>

#include "xbl.h"

void gauss(int /*nsiz*/, int nn, std::vector<std::vector<double>> &z, std::vector<std::vector<double>> &r, int nrhs) {
    for (int np = 1; np <= nn - 1; ++np) {
        const int np1 = np + 1;
        int nx = np;
        for (int n = np1; n <= nn; ++n) {
            if (std::abs(z[n][np]) > std::abs(z[nx][np])) {
                nx = n;
            }
        }
        const double pivot = 1.0 / z[nx][np];

        z[nx][np] = z[np][np];

        for (int l = np1; l <= nn; ++l) {
            const double temp = z[nx][l] * pivot;
            z[nx][l] = z[np][l];
            z[np][l] = temp;
        }

        for (int l = 1; l <= nrhs; ++l) {
            const double temp = r[nx][l] * pivot;
            r[nx][l] = r[np][l];
            r[np][l] = temp;
        }

        for (int k = np1; k <= nn; ++k) {
            const double ztmp = z[k][np];
            for (int l = np1; l <= nn; ++l) {
                z[k][l] = z[k][l] - ztmp * z[np][l];
            }
            for (int l = 1; l <= nrhs; ++l) {
                r[k][l] = r[k][l] - ztmp * r[np][l];
            }
        }
    }

    for (int l = 1; l <= nrhs; ++l) {
        r[nn][l] = r[nn][l] / z[nn][nn];
    }

    for (int np = nn - 1; np >= 1; --np) {
        const int np1 = np + 1;
        for (int l = 1; l <= nrhs; ++l) {
            for (int k = np1; k <= nn; ++k) {
                r[np][l] = r[np][l] - z[np][k] * r[k][l];
            }
        }
    }
}

void gauss(int /*nsiz*/, int nn, std::vector<std::vector<double>> &z, std::vector<double> &r, int /*nrhs*/) {
    for (int np = 1; np <= nn - 1; ++np) {
        const int np1 = np + 1;
        int nx = np;
        for (int n = np1; n <= nn; ++n) {
            if (std::abs(z[n][np]) > std::abs(z[nx][np])) {
                nx = n;
            }
        }
        const double pivot = 1.0 / z[nx][np];

        z[nx][np] = z[np][np];

        for (int l = np1; l <= nn; ++l) {
            const double temp = z[nx][l] * pivot;
            z[nx][l] = z[np][l];
            z[np][l] = temp;
        }

        const double temp = r[nx] * pivot;
        r[nx] = r[np];
        r[np] = temp;

        for (int k = np1; k <= nn; ++k) {
            const double ztmp = z[k][np];
            for (int l = np1; l <= nn; ++l) {
                z[k][l] = z[k][l] - ztmp * z[np][l];
            }
            r[k] = r[k] - ztmp * r[np];
        }
    }

    r[nn] = r[nn] / z[nn][nn];

    for (int np = nn - 1; np >= 1; --np) {
        const int np1 = np + 1;
        for (int k = np1; k <= nn; ++k) {
            r[np] = r[np] - z[np][k] * r[k];
        }
    }
}

void ludcmp(int /*nsiz*/, int n, std::vector<std::vector<double>> &a, std::vector<int> &indx) {
    const int nvx = 500;
    if (n > nvx) {
        throw std::runtime_error("LUDCMP: Array overflow. Increase NVX.");
    }

    std::vector<double> vv(static_cast<size_t>(nvx + 1), 0.0);

    for (int i = 1; i <= n; ++i) {
        double aamax = 0.0;
        for (int j = 1; j <= n; ++j) {
            aamax = std::max(std::abs(a[i][j]), aamax);
        }
        vv[i] = 1.0 / aamax;
    }

    for (int j = 1; j <= n; ++j) {
        for (int i = 1; i <= j - 1; ++i) {
            double summ = a[i][j];
            for (int k = 1; k <= i - 1; ++k) {
                summ = summ - a[i][k] * a[k][j];
            }
            a[i][j] = summ;
        }

        double aamax = 0.0;
        int imax = j;
        for (int i = j; i <= n; ++i) {
            double summ = a[i][j];
            for (int k = 1; k <= j - 1; ++k) {
                summ = summ - a[i][k] * a[k][j];
            }
            a[i][j] = summ;
            const double dum = vv[i] * std::abs(summ);
            if (dum >= aamax) {
                imax = i;
                aamax = dum;
            }
        }

        if (j != imax) {
            for (int k = 1; k <= n; ++k) {
                const double dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (a[j][j] == 0.0) {
            a[j][j] = 1.0e-20;
        }

        if (j != n) {
            const double dum = 1.0 / a[j][j];
            for (int i = j + 1; i <= n; ++i) {
                a[i][j] = a[i][j] * dum;
            }
        }
    }
}

void baksub(int /*nsiz*/, int n, std::vector<std::vector<double>> &a, const std::vector<int> &indx, std::vector<double> &b) {
    int ii = 0;
    for (int i = 1; i <= n; ++i) {
        const int ll = indx[i];
        double summ = b[ll];
        b[ll] = b[i];
        if (ii != 0) {
            for (int j = ii; j <= i - 1; ++j) {
                summ = summ - a[i][j] * b[j];
            }
        } else if (summ != 0.0) {
            ii = i;
        }
        b[i] = summ;
    }

    for (int i = n; i >= 1; --i) {
        double summ = b[i];
        if (i < n) {
            for (int j = i + 1; j <= n; ++j) {
                summ = summ - a[i][j] * b[j];
            }
        }
        b[i] = summ / a[i][i];
    }
}

void blsolv(XFoilState &ctx) {
    const int ivte1 = ctx.ISYS[ctx.IBLTE[1]][1];

    const double vacc1 = ctx.VACCEL;
    const double vacc2 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1]);
    const double vacc3 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1]);

    for (int iv = 1; iv <= ctx.NSYS; ++iv) {
        const int ivp = iv + 1;

        double pivot = 1.0 / ctx.VA[1][1][iv];
        ctx.VA[1][2][iv] = ctx.VA[1][2][iv] * pivot;
        for (int l = iv; l <= ctx.NSYS; ++l) {
            ctx.VM[1][l][iv] = ctx.VM[1][l][iv] * pivot;
        }
        ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] * pivot;
        ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] * pivot;

        for (int k = 2; k <= 3; ++k) {
            const double vtmp = ctx.VA[k][1][iv];
            ctx.VA[k][2][iv] = ctx.VA[k][2][iv] - vtmp * ctx.VA[1][2][iv];
            for (int l = iv; l <= ctx.NSYS; ++l) {
                ctx.VM[k][l][iv] = ctx.VM[k][l][iv] - vtmp * ctx.VM[1][l][iv];
            }
            ctx.VDEL[k][1][iv] = ctx.VDEL[k][1][iv] - vtmp * ctx.VDEL[1][1][iv];
            ctx.VDEL[k][2][iv] = ctx.VDEL[k][2][iv] - vtmp * ctx.VDEL[1][2][iv];
        }

        pivot = 1.0 / ctx.VA[2][2][iv];
        for (int l = iv; l <= ctx.NSYS; ++l) {
            ctx.VM[2][l][iv] = ctx.VM[2][l][iv] * pivot;
        }
        ctx.VDEL[2][1][iv] = ctx.VDEL[2][1][iv] * pivot;
        ctx.VDEL[2][2][iv] = ctx.VDEL[2][2][iv] * pivot;

        {
            const int k = 3;
            const double vtmp = ctx.VA[k][2][iv];
            for (int l = iv; l <= ctx.NSYS; ++l) {
                ctx.VM[k][l][iv] = ctx.VM[k][l][iv] - vtmp * ctx.VM[2][l][iv];
            }
            ctx.VDEL[k][1][iv] = ctx.VDEL[k][1][iv] - vtmp * ctx.VDEL[2][1][iv];
            ctx.VDEL[k][2][iv] = ctx.VDEL[k][2][iv] - vtmp * ctx.VDEL[2][2][iv];
        }

        pivot = 1.0 / ctx.VM[3][iv][iv];
        for (int l = ivp; l <= ctx.NSYS; ++l) {
            ctx.VM[3][l][iv] = ctx.VM[3][l][iv] * pivot;
        }
        ctx.VDEL[3][1][iv] = ctx.VDEL[3][1][iv] * pivot;
        ctx.VDEL[3][2][iv] = ctx.VDEL[3][2][iv] * pivot;

        const double vtmp1 = ctx.VM[1][iv][iv];
        const double vtmp2 = ctx.VM[2][iv][iv];
        for (int l = ivp; l <= ctx.NSYS; ++l) {
            ctx.VM[1][l][iv] = ctx.VM[1][l][iv] - vtmp1 * ctx.VM[3][l][iv];
            ctx.VM[2][l][iv] = ctx.VM[2][l][iv] - vtmp2 * ctx.VM[3][l][iv];
        }
        ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] - vtmp1 * ctx.VDEL[3][1][iv];
        ctx.VDEL[2][1][iv] = ctx.VDEL[2][1][iv] - vtmp2 * ctx.VDEL[3][1][iv];
        ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] - vtmp1 * ctx.VDEL[3][2][iv];
        ctx.VDEL[2][2][iv] = ctx.VDEL[2][2][iv] - vtmp2 * ctx.VDEL[3][2][iv];

        {
            const double vtmp = ctx.VA[1][2][iv];
            for (int l = ivp; l <= ctx.NSYS; ++l) {
                ctx.VM[1][l][iv] = ctx.VM[1][l][iv] - vtmp * ctx.VM[2][l][iv];
            }
            ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] - vtmp * ctx.VDEL[2][1][iv];
            ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] - vtmp * ctx.VDEL[2][2][iv];
        }

        if (iv == ctx.NSYS) {
            continue;
        }

        for (int k = 1; k <= 3; ++k) {
            const double vtmp1 = ctx.VB[k][1][ivp];
            const double vtmp2 = ctx.VB[k][2][ivp];
            const double vtmp3 = ctx.VM[k][iv][ivp];
            for (int l = ivp; l <= ctx.NSYS; ++l) {
                ctx.VM[k][l][ivp] = ctx.VM[k][l][ivp]
                                    - (vtmp1 * ctx.VM[1][l][iv] + vtmp2 * ctx.VM[2][l][iv] + vtmp3 * ctx.VM[3][l][iv]);
            }
            ctx.VDEL[k][1][ivp] = ctx.VDEL[k][1][ivp]
                                  - (vtmp1 * ctx.VDEL[1][1][iv] + vtmp2 * ctx.VDEL[2][1][iv] + vtmp3 * ctx.VDEL[3][1][iv]);
            ctx.VDEL[k][2][ivp] = ctx.VDEL[k][2][ivp]
                                  - (vtmp1 * ctx.VDEL[1][2][iv] + vtmp2 * ctx.VDEL[2][2][iv] + vtmp3 * ctx.VDEL[3][2][iv]);
        }

        if (iv == ivte1) {
            const int ivz = ctx.ISYS[ctx.IBLTE[2] + 1][2];
            for (int k = 1; k <= 3; ++k) {
                const double vtmp1 = ctx.VZ[k][1];
                const double vtmp2 = ctx.VZ[k][2];
                for (int l = ivp; l <= ctx.NSYS; ++l) {
                    ctx.VM[k][l][ivz] = ctx.VM[k][l][ivz] - (vtmp1 * ctx.VM[1][l][iv] + vtmp2 * ctx.VM[2][l][iv]);
                }
                ctx.VDEL[k][1][ivz] = ctx.VDEL[k][1][ivz] - (vtmp1 * ctx.VDEL[1][1][iv] + vtmp2 * ctx.VDEL[2][1][iv]);
                ctx.VDEL[k][2][ivz] = ctx.VDEL[k][2][ivz] - (vtmp1 * ctx.VDEL[1][2][iv] + vtmp2 * ctx.VDEL[2][2][iv]);
            }
        }

        if (ivp != ctx.NSYS) {
            for (int kv = iv + 2; kv <= ctx.NSYS; ++kv) {
                const double vtmp1 = ctx.VM[1][iv][kv];
                const double vtmp2 = ctx.VM[2][iv][kv];
                const double vtmp3 = ctx.VM[3][iv][kv];

                if (std::abs(vtmp1) > vacc1) {
                    for (int l = ivp; l <= ctx.NSYS; ++l) {
                        ctx.VM[1][l][kv] = ctx.VM[1][l][kv] - vtmp1 * ctx.VM[3][l][iv];
                    }
                    ctx.VDEL[1][1][kv] = ctx.VDEL[1][1][kv] - vtmp1 * ctx.VDEL[3][1][iv];
                    ctx.VDEL[1][2][kv] = ctx.VDEL[1][2][kv] - vtmp1 * ctx.VDEL[3][2][iv];
                }

                if (std::abs(vtmp2) > vacc2) {
                    for (int l = ivp; l <= ctx.NSYS; ++l) {
                        ctx.VM[2][l][kv] = ctx.VM[2][l][kv] - vtmp2 * ctx.VM[3][l][iv];
                    }
                    ctx.VDEL[2][1][kv] = ctx.VDEL[2][1][kv] - vtmp2 * ctx.VDEL[3][1][iv];
                    ctx.VDEL[2][2][kv] = ctx.VDEL[2][2][kv] - vtmp2 * ctx.VDEL[3][2][iv];
                }

                if (std::abs(vtmp3) > vacc3) {
                    for (int l = ivp; l <= ctx.NSYS; ++l) {
                        ctx.VM[3][l][kv] = ctx.VM[3][l][kv] - vtmp3 * ctx.VM[3][l][iv];
                    }
                    ctx.VDEL[3][1][kv] = ctx.VDEL[3][1][kv] - vtmp3 * ctx.VDEL[3][1][iv];
                    ctx.VDEL[3][2][kv] = ctx.VDEL[3][2][kv] - vtmp3 * ctx.VDEL[3][2][iv];
                }
            }
        }
    }

    for (int iv = ctx.NSYS; iv >= 2; --iv) {
        double vtmp = ctx.VDEL[3][1][iv];
        for (int kv = iv - 1; kv >= 1; --kv) {
            ctx.VDEL[1][1][kv] = ctx.VDEL[1][1][kv] - ctx.VM[1][iv][kv] * vtmp;
            ctx.VDEL[2][1][kv] = ctx.VDEL[2][1][kv] - ctx.VM[2][iv][kv] * vtmp;
            ctx.VDEL[3][1][kv] = ctx.VDEL[3][1][kv] - ctx.VM[3][iv][kv] * vtmp;
        }

        vtmp = ctx.VDEL[3][2][iv];
        for (int kv = iv - 1; kv >= 1; --kv) {
            ctx.VDEL[1][2][kv] = ctx.VDEL[1][2][kv] - ctx.VM[1][iv][kv] * vtmp;
            ctx.VDEL[2][2][kv] = ctx.VDEL[2][2][kv] - ctx.VM[2][iv][kv] * vtmp;
            ctx.VDEL[3][2][kv] = ctx.VDEL[3][2][kv] - ctx.VM[3][iv][kv] * vtmp;
        }
    }
}
