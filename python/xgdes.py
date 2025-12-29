# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

from .spline import sinvrt, seval


def getxyf(x, xp, y, yp, s, n, tops, bots, xf, yf):
    if xf == -999.0:
        raise RuntimeError("GETXYF: hinge x not specified.")

    tops = s[1] + (x[1] - xf)
    bots = s[n] - (x[n] - xf)
    tops = sinvrt(tops, xf, x, xp, s, n)
    bots = sinvrt(bots, xf, x, xp, s, n)
    topy = seval(tops, y, yp, s, n)
    boty = seval(bots, y, yp, s, n)

    print()
    print(f"  Top    surface:  y ={topy:8.4f}     y/t = 1.0")
    print(f"  Bottom surface:  y ={boty:8.4f}     y/t = 0.0")

    if yf == -999.0:
        raise RuntimeError("GETXYF: hinge y not specified.")

    if yf == 999.0:
        raise RuntimeError("GETXYF: y/t hinge input not supported.")

    return tops, bots, xf, yf
