import math


def setexp(s, ds1, smax, nn):
    sigma = smax / ds1
    nex = nn - 1
    rnex = float(nex)
    rni = 1.0 / rnex

    aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0
    bbb = rnex * (rnex - 1.0) / 2.0
    ccc = rnex - sigma

    disc = bbb**2 - 4.0 * aaa * ccc
    disc = max(0.0, disc)

    if nex <= 1:
        raise RuntimeError("SETEXP: Cannot fill array.  N too small.")
    elif nex == 2:
        ratio = -ccc / bbb + 1.0
    else:
        ratio = (-bbb + math.sqrt(disc)) / (2.0 * aaa) + 1.0

    if ratio != 1.0:
        for _ in range(1, 101):
            sigman = (ratio**nex - 1.0) / (ratio - 1.0)
            res = sigman**rni - sigma**rni
            dresdr = rni * sigman**rni * (rnex * ratio ** (nex - 1) - sigman) / (ratio**nex - 1.0)
            dratio = -res / dresdr
            ratio = ratio + dratio
            if abs(dratio) < 1.0e-5:
                break
        else:
            print("SETEXP: Convergence failed.  Continuing anyway ...")

    s[1] = 0.0
    ds = ds1
    for n in range(2, nn + 1):
        s[n] = s[n - 1] + ds
        ds = ds * ratio


def atanc(y, x, thold):
    pi = math.pi
    tpi = 2.0 * pi

    thnew = math.atan2(y, x)
    dthet = thnew - thold
    dtcorr = dthet - tpi * int((dthet + math.copysign(pi, dthet)) / tpi)
    return thold + dtcorr
