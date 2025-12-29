# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math


def trisol(a, b, c, d, kk):
    for k in range(2, kk + 1):
        km = k - 1
        c[km] = c[km] / a[km]
        d[km] = d[km] / a[km]
        a[k] = a[k] - b[k] * c[km]
        d[k] = d[k] - b[k] * d[km]

    d[kk] = d[kk] / a[kk]

    for k in range(kk - 1, 0, -1):
        d[k] = d[k] - c[k] * d[k + 1]


def splind(x, xs, s, n, xs1, xs2):
    nmax = 1000
    if n > nmax:
        raise RuntimeError("SPLIND: array overflow, increase NMAX")

    a = [0.0] * (n + 1)
    b = [0.0] * (n + 1)
    c = [0.0] * (n + 1)

    for i in range(2, n):
        dsm = s[i] - s[i - 1]
        dsp = s[i + 1] - s[i]
        b[i] = dsp
        a[i] = 2.0 * (dsm + dsp)
        c[i] = dsm
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm)

    if xs1 == 999.0:
        a[1] = 2.0
        c[1] = 1.0
        xs[1] = 3.0 * (x[2] - x[1]) / (s[2] - s[1])
    elif xs1 == -999.0:
        a[1] = 1.0
        c[1] = 1.0
        xs[1] = 2.0 * (x[2] - x[1]) / (s[2] - s[1])
    else:
        a[1] = 1.0
        c[1] = 0.0
        xs[1] = xs1

    if xs2 == 999.0:
        b[n] = 1.0
        a[n] = 2.0
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1])
    elif xs2 == -999.0:
        b[n] = 1.0
        a[n] = 1.0
        xs[n] = 2.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1])
    else:
        a[n] = 1.0
        b[n] = 0.0
        xs[n] = xs2

    if n == 2 and xs1 == -999.0 and xs2 == -999.0:
        b[n] = 1.0
        a[n] = 2.0
        xs[n] = 3.0 * (x[n] - x[n - 1]) / (s[n] - s[n - 1])

    trisol(a, b, c, xs, n)


def seval(ss, x, xs, s, n):
    ilow = 1
    i = n

    while i - ilow > 1:
        imid = (i + ilow) // 2
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid

    ds = s[i] - s[i - 1]
    t = (ss - s[i - 1]) / ds
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1]
    cx2 = ds * xs[i] - x[i] + x[i - 1]
    return t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)


def deval(ss, x, xs, s, n):
    ilow = 1
    i = n

    while i - ilow > 1:
        imid = (i + ilow) // 2
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid

    ds = s[i] - s[i - 1]
    t = (ss - s[i - 1]) / ds
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1]
    cx2 = ds * xs[i] - x[i] + x[i - 1]
    deval_val = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2
    return deval_val / ds


def sinvrt(si, xi, x, xs, s, n):
    sisav = si
    for _ in range(10):
        res = seval(si, x, xs, s, n) - xi
        resp = deval(si, x, xs, s, n)
        ds = -res / resp
        si = si + ds
        if abs(ds / (s[n] - s[1])) < 1.0e-5:
            return si
    print()
    print("SINVRT: spline inversion failed. Input value returned.")
    return sisav


def scalc(x, y, s, n):
    s[1] = 0.0
    for i in range(2, n + 1):
        s[i] = s[i - 1] + math.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2)


def segspl(x, xs, s, n):
    if s[1] == s[2]:
        raise RuntimeError("SEGSPL:  First input point duplicated")
    if s[n] == s[n - 1]:
        raise RuntimeError("SEGSPL:  Last  input point duplicated")

    iseg0 = 1
    for iseg in range(2, n - 1):
        if s[iseg] == s[iseg + 1]:
            nseg = iseg - iseg0 + 1
            xt = [0.0] * (nseg + 1)
            xst = [0.0] * (nseg + 1)
            st = [0.0] * (nseg + 1)
            for i in range(1, nseg + 1):
                xt[i] = x[iseg0 + i - 1]
                xst[i] = xs[iseg0 + i - 1]
                st[i] = s[iseg0 + i - 1]
            splind(xt, xst, st, nseg, -999.0, -999.0)
            for i in range(1, nseg + 1):
                xs[iseg0 + i - 1] = xst[i]
            iseg0 = iseg + 1

    nseg = n - iseg0 + 1
    xt = [0.0] * (nseg + 1)
    xst = [0.0] * (nseg + 1)
    st = [0.0] * (nseg + 1)
    for i in range(1, nseg + 1):
        xt[i] = x[iseg0 + i - 1]
        xst[i] = xs[iseg0 + i - 1]
        st[i] = s[iseg0 + i - 1]
    splind(xt, xst, st, nseg, -999.0, -999.0)
    for i in range(1, nseg + 1):
        xs[iseg0 + i - 1] = xst[i]
def d2val(ss, x, xs, s, n):
    ilow = 1
    i = n

    while i - ilow > 1:
        imid = (i + ilow) // 2
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid

    ds = s[i] - s[i - 1]
    t = (ss - s[i - 1]) / ds
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1]
    cx2 = ds * xs[i] - x[i] + x[i - 1]
    d2val_val = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2
    return d2val_val / ds**2


def curv(ss, x, xs, y, ys, s, n):
    ilow = 1
    i = n

    while i - ilow > 1:
        imid = (i + ilow) // 2
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid

    ds = s[i] - s[i - 1]
    t = (ss - s[i - 1]) / ds

    cx1 = ds * xs[i - 1] - x[i] + x[i - 1]
    cx2 = ds * xs[i] - x[i] + x[i - 1]
    xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2
    xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2

    cy1 = ds * ys[i - 1] - y[i] + y[i - 1]
    cy2 = ds * ys[i] - y[i] + y[i - 1]
    yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2
    ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2

    sd = math.sqrt(xd * xd + yd * yd)
    sd = max(sd, 0.001 * ds)

    return (xd * ydd - yd * xdd) / sd**3
