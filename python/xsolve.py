def gauss(nsiz, nn, z, r, nrhs):
    for np in range(1, nn):
        np1 = np + 1
        nx = np
        for n in range(np1, nn + 1):
            if abs(z[n][np]) > abs(z[nx][np]):
                nx = n
        pivot = 1.0 / z[nx][np]

        z[nx][np] = z[np][np]

        for l in range(np1, nn + 1):
            temp = z[nx][l] * pivot
            z[nx][l] = z[np][l]
            z[np][l] = temp

        for l in range(1, nrhs + 1):
            temp = r[nx][l] * pivot
            r[nx][l] = r[np][l]
            r[np][l] = temp

        for k in range(np1, nn + 1):
            ztmp = z[k][np]
            for l in range(np1, nn + 1):
                z[k][l] = z[k][l] - ztmp * z[np][l]
            for l in range(1, nrhs + 1):
                r[k][l] = r[k][l] - ztmp * r[np][l]

    for l in range(1, nrhs + 1):
        r[nn][l] = r[nn][l] / z[nn][nn]

    for np in range(nn - 1, 0, -1):
        np1 = np + 1
        for l in range(1, nrhs + 1):
            for k in range(np1, nn + 1):
                r[np][l] = r[np][l] - z[np][k] * r[k][l]


def ludcmp(nsiz, n, a, indx):
    nvx = 500
    if n > nvx:
        raise RuntimeError("LUDCMP: Array overflow. Increase NVX.")

    vv = [0.0] * (nvx + 1)

    for i in range(1, n + 1):
        aamax = 0.0
        for j in range(1, n + 1):
            aamax = max(abs(a[i][j]), aamax)
        vv[i] = 1.0 / aamax

    for j in range(1, n + 1):
        for i in range(1, j):
            summ = a[i][j]
            for k in range(1, i):
                summ = summ - a[i][k] * a[k][j]
            a[i][j] = summ

        aamax = 0.0
        imax = j
        for i in range(j, n + 1):
            summ = a[i][j]
            for k in range(1, j):
                summ = summ - a[i][k] * a[k][j]
            a[i][j] = summ
            dum = vv[i] * abs(summ)
            if dum >= aamax:
                imax = i
                aamax = dum

        if j != imax:
            for k in range(1, n + 1):
                dum = a[imax][k]
                a[imax][k] = a[j][k]
                a[j][k] = dum
            vv[imax] = vv[j]

        indx[j] = imax
        if a[j][j] == 0.0:
            a[j][j] = 1.0e-20

        if j != n:
            dum = 1.0 / a[j][j]
            for i in range(j + 1, n + 1):
                a[i][j] = a[i][j] * dum


def baksub(nsiz, n, a, indx, b):
    ii = 0
    for i in range(1, n + 1):
        ll = indx[i]
        summ = b[ll]
        b[ll] = b[i]
        if ii != 0:
            for j in range(ii, i):
                summ = summ - a[i][j] * b[j]
        elif summ != 0.0:
            ii = i
        b[i] = summ

    for i in range(n, 0, -1):
        summ = b[i]
        if i < n:
            for j in range(i + 1, n + 1):
                summ = summ - a[i][j] * b[j]
        b[i] = summ / a[i][i]


def blsolv(ctx):
    ivte1 = ctx.ISYS[ctx.IBLTE[1]][1]

    vacc1 = ctx.VACCEL
    vacc2 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1])
    vacc3 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1])

    for iv in range(1, ctx.NSYS + 1):
        ivp = iv + 1

        pivot = 1.0 / ctx.VA[1][1][iv]
        ctx.VA[1][2][iv] = ctx.VA[1][2][iv] * pivot
        for l in range(iv, ctx.NSYS + 1):
            ctx.VM[1][l][iv] = ctx.VM[1][l][iv] * pivot
        ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] * pivot
        ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] * pivot

        for k in range(2, 4):
            vtmp = ctx.VA[k][1][iv]
            ctx.VA[k][2][iv] = ctx.VA[k][2][iv] - vtmp * ctx.VA[1][2][iv]
            for l in range(iv, ctx.NSYS + 1):
                ctx.VM[k][l][iv] = ctx.VM[k][l][iv] - vtmp * ctx.VM[1][l][iv]
            ctx.VDEL[k][1][iv] = ctx.VDEL[k][1][iv] - vtmp * ctx.VDEL[1][1][iv]
            ctx.VDEL[k][2][iv] = ctx.VDEL[k][2][iv] - vtmp * ctx.VDEL[1][2][iv]

        pivot = 1.0 / ctx.VA[2][2][iv]
        for l in range(iv, ctx.NSYS + 1):
            ctx.VM[2][l][iv] = ctx.VM[2][l][iv] * pivot
        ctx.VDEL[2][1][iv] = ctx.VDEL[2][1][iv] * pivot
        ctx.VDEL[2][2][iv] = ctx.VDEL[2][2][iv] * pivot

        k = 3
        vtmp = ctx.VA[k][2][iv]
        for l in range(iv, ctx.NSYS + 1):
            ctx.VM[k][l][iv] = ctx.VM[k][l][iv] - vtmp * ctx.VM[2][l][iv]
        ctx.VDEL[k][1][iv] = ctx.VDEL[k][1][iv] - vtmp * ctx.VDEL[2][1][iv]
        ctx.VDEL[k][2][iv] = ctx.VDEL[k][2][iv] - vtmp * ctx.VDEL[2][2][iv]

        pivot = 1.0 / ctx.VM[3][iv][iv]
        for l in range(ivp, ctx.NSYS + 1):
            ctx.VM[3][l][iv] = ctx.VM[3][l][iv] * pivot
        ctx.VDEL[3][1][iv] = ctx.VDEL[3][1][iv] * pivot
        ctx.VDEL[3][2][iv] = ctx.VDEL[3][2][iv] * pivot

        vtmp1 = ctx.VM[1][iv][iv]
        vtmp2 = ctx.VM[2][iv][iv]
        for l in range(ivp, ctx.NSYS + 1):
            ctx.VM[1][l][iv] = ctx.VM[1][l][iv] - vtmp1 * ctx.VM[3][l][iv]
            ctx.VM[2][l][iv] = ctx.VM[2][l][iv] - vtmp2 * ctx.VM[3][l][iv]
        ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] - vtmp1 * ctx.VDEL[3][1][iv]
        ctx.VDEL[2][1][iv] = ctx.VDEL[2][1][iv] - vtmp2 * ctx.VDEL[3][1][iv]
        ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] - vtmp1 * ctx.VDEL[3][2][iv]
        ctx.VDEL[2][2][iv] = ctx.VDEL[2][2][iv] - vtmp2 * ctx.VDEL[3][2][iv]

        vtmp = ctx.VA[1][2][iv]
        for l in range(ivp, ctx.NSYS + 1):
            ctx.VM[1][l][iv] = ctx.VM[1][l][iv] - vtmp * ctx.VM[2][l][iv]
        ctx.VDEL[1][1][iv] = ctx.VDEL[1][1][iv] - vtmp * ctx.VDEL[2][1][iv]
        ctx.VDEL[1][2][iv] = ctx.VDEL[1][2][iv] - vtmp * ctx.VDEL[2][2][iv]

        if iv == ctx.NSYS:
            continue

        for k in range(1, 4):
            vtmp1 = ctx.VB[k][1][ivp]
            vtmp2 = ctx.VB[k][2][ivp]
            vtmp3 = ctx.VM[k][iv][ivp]
            for l in range(ivp, ctx.NSYS + 1):
                ctx.VM[k][l][ivp] = ctx.VM[k][l][ivp] - (
                    vtmp1 * ctx.VM[1][l][iv] + vtmp2 * ctx.VM[2][l][iv] + vtmp3 * ctx.VM[3][l][iv]
                )
            ctx.VDEL[k][1][ivp] = ctx.VDEL[k][1][ivp] - (
                vtmp1 * ctx.VDEL[1][1][iv] + vtmp2 * ctx.VDEL[2][1][iv] + vtmp3 * ctx.VDEL[3][1][iv]
            )
            ctx.VDEL[k][2][ivp] = ctx.VDEL[k][2][ivp] - (
                vtmp1 * ctx.VDEL[1][2][iv] + vtmp2 * ctx.VDEL[2][2][iv] + vtmp3 * ctx.VDEL[3][2][iv]
            )

        if iv == ivte1:
            ivz = ctx.ISYS[ctx.IBLTE[2] + 1][2]
            for k in range(1, 4):
                vtmp1 = ctx.VZ[k][1]
                vtmp2 = ctx.VZ[k][2]
                for l in range(ivp, ctx.NSYS + 1):
                    ctx.VM[k][l][ivz] = ctx.VM[k][l][ivz] - (
                        vtmp1 * ctx.VM[1][l][iv] + vtmp2 * ctx.VM[2][l][iv]
                    )
                ctx.VDEL[k][1][ivz] = ctx.VDEL[k][1][ivz] - (
                    vtmp1 * ctx.VDEL[1][1][iv] + vtmp2 * ctx.VDEL[2][1][iv]
                )
                ctx.VDEL[k][2][ivz] = ctx.VDEL[k][2][ivz] - (
                    vtmp1 * ctx.VDEL[1][2][iv] + vtmp2 * ctx.VDEL[2][2][iv]
                )

        if ivp != ctx.NSYS:
            for kv in range(iv + 2, ctx.NSYS + 1):
                vtmp1 = ctx.VM[1][iv][kv]
                vtmp2 = ctx.VM[2][iv][kv]
                vtmp3 = ctx.VM[3][iv][kv]

                if abs(vtmp1) > vacc1:
                    for l in range(ivp, ctx.NSYS + 1):
                        ctx.VM[1][l][kv] = ctx.VM[1][l][kv] - vtmp1 * ctx.VM[3][l][iv]
                    ctx.VDEL[1][1][kv] = ctx.VDEL[1][1][kv] - vtmp1 * ctx.VDEL[3][1][iv]
                    ctx.VDEL[1][2][kv] = ctx.VDEL[1][2][kv] - vtmp1 * ctx.VDEL[3][2][iv]

                if abs(vtmp2) > vacc2:
                    for l in range(ivp, ctx.NSYS + 1):
                        ctx.VM[2][l][kv] = ctx.VM[2][l][kv] - vtmp2 * ctx.VM[3][l][iv]
                    ctx.VDEL[2][1][kv] = ctx.VDEL[2][1][kv] - vtmp2 * ctx.VDEL[3][1][iv]
                    ctx.VDEL[2][2][kv] = ctx.VDEL[2][2][kv] - vtmp2 * ctx.VDEL[3][2][iv]

                if abs(vtmp3) > vacc3:
                    for l in range(ivp, ctx.NSYS + 1):
                        ctx.VM[3][l][kv] = ctx.VM[3][l][kv] - vtmp3 * ctx.VM[3][l][iv]
                    ctx.VDEL[3][1][kv] = ctx.VDEL[3][1][kv] - vtmp3 * ctx.VDEL[3][1][iv]
                    ctx.VDEL[3][2][kv] = ctx.VDEL[3][2][kv] - vtmp3 * ctx.VDEL[3][2][iv]

    for iv in range(ctx.NSYS, 1, -1):
        vtmp = ctx.VDEL[3][1][iv]
        for kv in range(iv - 1, 0, -1):
            ctx.VDEL[1][1][kv] = ctx.VDEL[1][1][kv] - ctx.VM[1][iv][kv] * vtmp
            ctx.VDEL[2][1][kv] = ctx.VDEL[2][1][kv] - ctx.VM[2][iv][kv] * vtmp
            ctx.VDEL[3][1][kv] = ctx.VDEL[3][1][kv] - ctx.VM[3][iv][kv] * vtmp

        vtmp = ctx.VDEL[3][2][iv]
        for kv in range(iv - 1, 0, -1):
            ctx.VDEL[1][2][kv] = ctx.VDEL[1][2][kv] - ctx.VM[1][iv][kv] * vtmp
            ctx.VDEL[2][2][kv] = ctx.VDEL[2][2][kv] - ctx.VM[2][iv][kv] * vtmp
            ctx.VDEL[3][2][kv] = ctx.VDEL[3][2][kv] - ctx.VM[3][iv][kv] * vtmp
