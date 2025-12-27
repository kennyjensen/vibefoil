import math

from .spline import curv, d2val, deval, seval, sinvrt
from .xutils import atanc


def lefind(x, xp, y, yp, s, n):
    dseps = (s[n] - s[1]) * 1.0e-5

    xte = 0.5 * (x[1] + x[n])
    yte = 0.5 * (y[1] + y[n])

    i = 3
    for i in range(3, n - 1):
        dxte = x[i] - xte
        dyte = y[i] - yte
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        dotp = dxte * dx + dyte * dy
        if dotp < 0.0:
            break

    sle = s[i]

    if s[i] == s[i - 1]:
        return sle

    for _ in range(1, 51):
        xle = seval(sle, x, xp, s, n)
        yle = seval(sle, y, yp, s, n)
        dxds = deval(sle, x, xp, s, n)
        dyds = deval(sle, y, yp, s, n)
        dxdd = d2val(sle, x, xp, s, n)
        dydd = d2val(sle, y, yp, s, n)

        xchord = xle - xte
        ychord = yle - yte

        res = xchord * dxds + ychord * dyds
        ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd

        dsle = -res / ress

        dsle = max(dsle, -0.02 * abs(xchord + ychord))
        dsle = min(dsle, 0.02 * abs(xchord + ychord))
        sle = sle + dsle
        if abs(dsle) < dseps:
            return sle

    print("LEFIND:  LE point not found.  Continuing...")
    return s[i]


def sopps(si, x, xp, y, yp, s, n, sle):
    slen = s[n] - s[1]

    xle = seval(sle, x, xp, s, n)
    yle = seval(sle, y, yp, s, n)
    xte = 0.5 * (x[1] + x[n])
    yte = 0.5 * (y[1] + y[n])
    chord = math.sqrt((xte - xle) ** 2 + (yte - yle) ** 2)
    dxc = (xte - xle) / chord
    dyc = (yte - yle) / chord

    if si < sle:
        inp = 1
        inopp = n
    else:
        inp = n
        inopp = 1
    sfrac = (si - sle) / (s[inp] - sle)
    sopp = sle + sfrac * (s[inopp] - sle)

    if abs(sfrac) <= 1.0e-5:
        return sle

    xi = seval(si, x, xp, s, n)
    yi = seval(si, y, yp, s, n)
    xle = seval(sle, x, xp, s, n)
    yle = seval(sle, y, yp, s, n)
    xbar = (xi - xle) * dxc + (yi - yle) * dyc

    for _ in range(1, 13):
        xopp = seval(sopp, x, xp, s, n)
        yopp = seval(sopp, y, yp, s, n)
        xoppd = deval(sopp, x, xp, s, n)
        yoppd = deval(sopp, y, yp, s, n)

        res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar
        resd = xoppd * dxc + yoppd * dyc

        if abs(res) / slen < 1.0e-5:
            return sopp
        if resd == 0.0:
            break

        dsopp = -res / resd
        sopp = sopp + dsopp

        if abs(dsopp) / slen < 1.0e-5:
            return sopp

    print()
    print("SOPPS: Opposite-point location failed. Continuing...")
    return sle + sfrac * (s[inopp] - sle)


def aecalc(n, x, y, t, itype):
    pi = math.pi

    sint = 0.0
    aint = 0.0
    xint = 0.0
    yint = 0.0
    xxint = 0.0
    xyint = 0.0
    yyint = 0.0

    for io in range(1, n + 1):
        if io == n:
            ip = 1
        else:
            ip = io + 1

        dx = x[io] - x[ip]
        dy = y[io] - y[ip]
        xa = (x[io] + x[ip]) * 0.50
        ya = (y[io] + y[ip]) * 0.50
        ta = (t[io] + t[ip]) * 0.50

        ds = math.sqrt(dx * dx + dy * dy)
        sint = sint + ds

        if itype == 1:
            da = ya * dx
            aint = aint + da
            xint = xint + xa * da
            yint = yint + ya * da / 2.0
            xxint = xxint + xa * xa * da
            xyint = xyint + xa * ya * da / 2.0
            yyint = yyint + ya * ya * da / 3.0
        else:
            da = ta * ds
            aint = aint + da
            xint = xint + xa * da
            yint = yint + ya * da
            xxint = xxint + xa * xa * da
            xyint = xyint + xa * ya * da
            yyint = yyint + ya * ya * da

    area = aint

    if aint == 0.0:
        xcen = 0.0
        ycen = 0.0
        ei11 = 0.0
        ei22 = 0.0
        apx1 = 0.0
        apx2 = math.atan2(1.0, 0.0)
        return area, xcen, ycen, ei11, ei22, apx1, apx2

    xcen = xint / aint
    ycen = yint / aint

    eixx = yyint - ycen * ycen * aint
    eixy = xyint - xcen * ycen * aint
    eiyy = xxint - xcen * xcen * aint

    eisq = 0.25 * (eixx - eiyy) ** 2 + eixy**2
    sgn = math.copysign(1.0, eiyy - eixx)
    ei11 = 0.5 * (eixx + eiyy) - sgn * math.sqrt(eisq)
    ei22 = 0.5 * (eixx + eiyy) + sgn * math.sqrt(eisq)

    if ei11 == 0.0 or ei22 == 0.0:
        apx1 = 0.0
        apx2 = math.atan2(1.0, 0.0)
    elif eisq / (ei11 * ei22) < (0.001 * sint) ** 4:
        apx1 = 0.0
        apx2 = math.atan2(1.0, 0.0)
    else:
        c1 = eixy
        s1 = eixx - ei11

        c2 = eixy
        s2 = eixx - ei22

        if abs(s1) > abs(s2):
            apx1 = math.atan2(s1, c1)
            apx2 = apx1 + 0.5 * pi
        else:
            apx2 = math.atan2(s2, c2)
            apx1 = apx2 - 0.5 * pi

        if apx1 < -0.5 * pi:
            apx1 = apx1 + pi
        if apx1 > 0.5 * pi:
            apx1 = apx1 - pi
        if apx2 < -0.5 * pi:
            apx2 = apx2 + pi
        if apx2 > 0.5 * pi:
            apx2 = apx2 - pi

    return area, xcen, ycen, ei11, ei22, apx1, apx2


def tccalc(x, xp, y, yp, s, n):
    sle = lefind(x, xp, y, yp, s, n)
    xle = seval(sle, x, xp, s, n)
    yle = seval(sle, y, yp, s, n)
    xte = 0.5 * (x[1] + x[n])
    yte = 0.5 * (y[1] + y[n])
    chord = math.sqrt((xte - xle) ** 2 + (yte - yle) ** 2)

    dxc = (xte - xle) / chord
    dyc = (yte - yle) / chord

    thick = 0.0
    xthick = 0.0
    cambr = 0.0
    xcambr = 0.0

    for i in range(1, n + 1):
        xbar = (x[i] - xle) * dxc + (y[i] - yle) * dyc
        ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc

        sopp = sopps(s[i], x, xp, y, yp, s, n, sle)
        xopp = seval(sopp, x, xp, s, n)
        yopp = seval(sopp, y, yp, s, n)

        ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc

        yc = 0.5 * (ybar + ybarop)
        yt = abs(ybar - ybarop)

        if abs(yc) > abs(cambr):
            cambr = yc
            xcambr = xopp
        if abs(yt) > abs(thick):
            thick = yt
            xthick = xopp

    return thick, xthick, cambr, xcambr


def geopar(x, xp, y, yp, s, n, t):
    sle = lefind(x, xp, y, yp, s, n)

    xle = seval(sle, x, xp, s, n)
    yle = seval(sle, y, yp, s, n)
    xte = 0.5 * (x[1] + x[n])
    yte = 0.5 * (y[1] + y[n])

    chsq = (xte - xle) ** 2 + (yte - yle) ** 2
    chord = math.sqrt(chsq)

    curvle = curv(sle, x, xp, y, yp, s, n)

    radle = 0.0
    if abs(curvle) > 0.001 * (s[n] - s[1]):
        radle = 1.0 / curvle

    ang1 = math.atan2(-yp[1], -xp[1])
    ang2 = atanc(yp[n], xp[n], ang1)
    angte = ang2 - ang1

    for i in range(1, n + 1):
        t[i] = 1.0

    area, xcena, ycena, ei11a, ei22a, apx1a, apx2a = aecalc(n, x, y, t, 1)
    slen, xcent, ycent, ei11t, ei22t, apx1t, apx2t = aecalc(n, x, y, t, 2)

    thick, xthick, cambr, xcambr = tccalc(x, xp, y, yp, s, n)

    print(f" Max thickness = {thick:12.6f}  at x = {xthick:7.3f}")
    print(f" Max camber    = {cambr:12.6f}  at x = {xcambr:7.3f}")

    return {
        "sle": sle,
        "chord": chord,
        "area": area,
        "radle": radle,
        "angte": angte,
        "ei11a": ei11a,
        "ei22a": ei22a,
        "apx1a": apx1a,
        "apx2a": apx2a,
        "ei11t": ei11t,
        "ei22t": ei22t,
        "apx1t": apx1t,
        "apx2t": apx2t,
        "thick": thick,
        "cambr": cambr,
    }
