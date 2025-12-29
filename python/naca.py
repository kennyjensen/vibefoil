# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math


def naca4(ides, xx, yt, yc, nside, xb, yb):
    digits = "0123456789"
    an = 1.5

    n4 = ides // 1000
    n3 = (ides - n4 * 1000) // 100
    n2 = (ides - n4 * 1000 - n3 * 100) // 10
    n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10

    m = float(n4) / 100.0
    p = float(n3) / 10.0
    t = float(n2 * 10 + n1) / 100.0

    anp = an + 1.0
    for i in range(1, nside + 1):
        frac = float(i - 1) / float(nside - 1)
        if i == nside:
            xx[i] = 1.0
        else:
            xx[i] = 1.0 - anp * frac * (1.0 - frac) ** an - (1.0 - frac) ** anp

        yt[i] = (
            0.29690 * math.sqrt(xx[i])
            - 0.12600 * xx[i]
            - 0.35160 * xx[i] ** 2
            + 0.28430 * xx[i] ** 3
            - 0.10150 * xx[i] ** 4
        ) * t / 0.20

        if xx[i] < p:
            yc[i] = m / p**2 * (2.0 * p * xx[i] - xx[i] ** 2)
        else:
            yc[i] = m / (1.0 - p) ** 2 * ((1.0 - 2.0 * p) + 2.0 * p * xx[i] - xx[i] ** 2)

    ib = 0
    for i in range(nside, 0, -1):
        ib += 1
        xb[ib] = xx[i]
        yb[ib] = yc[i] + yt[i]
    for i in range(2, nside + 1):
        ib += 1
        xb[ib] = xx[i]
        yb[ib] = yc[i] - yt[i]

    name = "NACA"
    name = name.ljust(10)
    name = name[:5] + digits[n4 : n4 + 1] + digits[n3 : n3 + 1] + digits[n2 : n2 + 1] + digits[n1 : n1 + 1]
    return ib, name


def naca5(ides, xx, yt, yc, nside, xb, yb):
    digits = "0123456789"
    an = 1.5

    n5 = ides // 10000
    n4 = (ides - n5 * 10000) // 1000
    n3 = (ides - n5 * 10000 - n4 * 1000) // 100
    n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) // 10
    n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10

    n543 = 100 * n5 + 10 * n4 + n3

    if n543 == 210:
        m = 0.0580
        c = 361.4
    elif n543 == 220:
        m = 0.1260
        c = 51.64
    elif n543 == 230:
        m = 0.2025
        c = 15.957
    elif n543 == 240:
        m = 0.2900
        c = 6.643
    elif n543 == 250:
        m = 0.3910
        c = 3.230
    else:
        print("Illegal 5-digit designation")
        print("First three digits must be 210, 220, ... 250")
        return 0, ""

    t = float(n2 * 10 + n1) / 100.0

    anp = an + 1.0
    for i in range(1, nside + 1):
        frac = float(i - 1) / float(nside - 1)
        if i == nside:
            xx[i] = 1.0
        else:
            xx[i] = 1.0 - anp * frac * (1.0 - frac) ** an - (1.0 - frac) ** anp

        yt[i] = (
            0.29690 * math.sqrt(xx[i])
            - 0.12600 * xx[i]
            - 0.35160 * xx[i] ** 2
            + 0.28430 * xx[i] ** 3
            - 0.10150 * xx[i] ** 4
        ) * t / 0.20

        if xx[i] < m:
            yc[i] = (c / 6.0) * (xx[i] ** 3 - 3.0 * m * xx[i] ** 2 + m * m * (3.0 - m) * xx[i])
        else:
            yc[i] = (c / 6.0) * m**3 * (1.0 - xx[i])

    ib = 0
    for i in range(nside, 0, -1):
        ib += 1
        xb[ib] = xx[i]
        yb[ib] = yc[i] + yt[i]
    for i in range(2, nside + 1):
        ib += 1
        xb[ib] = xx[i]
        yb[ib] = yc[i] - yt[i]

    name = "NACA"
    name = name.ljust(10)
    name = (
        name[:5]
        + digits[n5 : n5 + 1]
        + digits[n4 : n4 + 1]
        + digits[n3 : n3 + 1]
        + digits[n2 : n2 + 1]
        + digits[n1 : n1 + 1]
    )
    return ib, name
