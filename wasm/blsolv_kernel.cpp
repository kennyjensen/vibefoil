// Port of the flat BLSOLV path in js/xsolve.js, derived from XFOIL by Mark Drela.
// GNU General Public License v2 or later.
// Preserve double precision, operation order, pivot checks, and VACCEL cutoffs.
// Build with build_blsolv.py: no reassociation or fused multiply-add.
// The source is row iv of plane 3; destinations are row kv >= iv + 2
// in distinct planes. Every span covers columns iv + 1 through nsys.
// They cannot overlap, so SIMD can update independent coefficients without
// changing the multiply/subtract order within any coefficient.
static inline void subtractScaledRows1(
    const double* __restrict source,
    double* __restrict dest1,
    double factor1, int count) {
  for (int j = 0; j < count; ++j) {
    const double value = source[j];
    dest1[j] -= factor1 * value;
  }
}

static inline void subtractScaledRows2(
    const double* __restrict source,
    double* __restrict dest1, double* __restrict dest2,
    double factor1, double factor2, int count) {
  for (int j = 0; j < count; ++j) {
    const double value = source[j];
    dest1[j] -= factor1 * value;
    dest2[j] -= factor2 * value;
  }
}

static inline void subtractScaledRows3(
    const double* __restrict source,
    double* __restrict dest1, double* __restrict dest2, double* __restrict dest3,
    double factor1, double factor2, double factor3, int count) {
  for (int j = 0; j < count; ++j) {
    const double value = source[j];
    dest1[j] -= factor1 * value;
    dest2[j] -= factor2 * value;
    dest3[j] -= factor3 * value;
  }
}

static void transposeMassCoefficients(double* vm, int nsys) {
  const int stride = nsys + 1;
  const int plane = stride * stride;
  for (int k = 1; k <= 3; k += 1) {
    const int base = k * plane;
    for (int row = 1; row < nsys; row += 1) {
      const int rowBase = base + row * stride;
      for (int col = row + 1; col <= nsys; col += 1) {
        const int a = rowBase + col;
        const int b = base + col * stride + row;
        const double value = vm[a];
        vm[a] = vm[b];
        vm[b] = value;
      }
    }
  }
}

static bool solve(int nsys, int ivte1, int ivz, double vaccel, double s1, double sN,
                  double* VAF, const double* VBF, double* VMF, double* VDELF, const double* VZ) {
    const double vacc1 = vaccel;
    const double vaccScale = vaccel * 2.0 / (sN - s1);
    const double vacc2 = vaccScale;
    const double vacc3 = vaccScale;
    const int vStride3 = nsys + 1;
    const int vStride2 = 3 * vStride3;
    const int vmStride3 = vStride3;
    const int vmStride2 = (nsys + 1) * vmStride3;
    const int vaBase1 = vStride2;
    const int vaBase2 = 2 * vStride2;
    const int vaBase3 = 3 * vStride2;
    const int vbBase1 = vStride2;
    const int vbBase2 = 2 * vStride2;
    const int vbBase3 = 3 * vStride2;
    const int vdelBase1 = vStride2;
    const int vdelBase2 = 2 * vStride2;
    const int vdelBase3 = 3 * vStride2;
    const int vmBase1 = vmStride2;
    const int vmBase2 = 2 * vmStride2;
    const int vmBase3 = 3 * vmStride2;
    const int va11Base = vaBase1 + vStride3;
    const int va12Base = vaBase1 + 2 * vStride3;
    const int va21Base = vaBase2 + vStride3;
    const int va22Base = vaBase2 + 2 * vStride3;
    const int va31Base = vaBase3 + vStride3;
    const int va32Base = vaBase3 + 2 * vStride3;
    const int vb11Base = vbBase1 + vStride3;
    const int vb12Base = vbBase1 + 2 * vStride3;
    const int vb21Base = vbBase2 + vStride3;
    const int vb22Base = vbBase2 + 2 * vStride3;
    const int vb31Base = vbBase3 + vStride3;
    const int vb32Base = vbBase3 + 2 * vStride3;
    const int vdel11Base = vdelBase1 + vStride3;
    const int vdel12Base = vdelBase1 + 2 * vStride3;
    const int vdel21Base = vdelBase2 + vStride3;
    const int vdel22Base = vdelBase2 + 2 * vStride3;
    const int vdel31Base = vdelBase3 + vStride3;
    const int vdel32Base = vdelBase3 + 2 * vStride3;

    for (int iv = 1; iv <= nsys; iv += 1) {
      const int ivp = iv + 1;

      const double va11 = VAF[va11Base + iv];
      const double va22 = VAF[va22Base + iv];
      const double vm33 = VMF[vmBase3 + iv * vmStride3 + iv];
      if (!__builtin_isfinite(va11) || !__builtin_isfinite(va22) || !__builtin_isfinite(vm33)
        || va11 == 0.0 || va22 == 0.0 || vm33 == 0.0) {
        return false;
      }

      double pivot = 1.0 / va11;
      VAF[va12Base + iv] *= pivot;
      for (int l = iv; l <= nsys; l += 1) {
        VMF[vmBase1 + iv * vmStride3 + l] *= pivot;
      }
      VDELF[vdel11Base + iv] *= pivot;
      VDELF[vdel12Base + iv] *= pivot;

      for (int k = 2; k <= 3; k += 1) {
        const int vaBase = k == 2 ? vaBase2 : vaBase3;
        const int vdelBase = k == 2 ? vdelBase2 : vdelBase3;
        const int vmBase = k == 2 ? vmBase2 : vmBase3;
        const double vtmp = VAF[vaBase + vStride3 + iv];
        VAF[vaBase + 2 * vStride3 + iv] -= vtmp * VAF[va12Base + iv];
        for (int l = iv; l <= nsys; l += 1) {
          const int vmkl = vmBase + iv * vmStride3 + l;
          const int vm1l = vmBase1 + iv * vmStride3 + l;
          VMF[vmkl] -= vtmp * VMF[vm1l];
        }
        VDELF[vdelBase + vStride3 + iv] -= vtmp * VDELF[vdel11Base + iv];
        VDELF[vdelBase + 2 * vStride3 + iv] -= vtmp * VDELF[vdel12Base + iv];
      }

      pivot = 1.0 / VAF[va22Base + iv];
      for (int l = iv; l <= nsys; l += 1) {
        VMF[vmBase2 + iv * vmStride3 + l] *= pivot;
      }
      VDELF[vdel21Base + iv] *= pivot;
      VDELF[vdel22Base + iv] *= pivot;

      {
        const double vtmp = VAF[va32Base + iv];
        for (int l = iv; l <= nsys; l += 1) {
          const int vm3l = vmBase3 + iv * vmStride3 + l;
          const int vm2l = vmBase2 + iv * vmStride3 + l;
          VMF[vm3l] -= vtmp * VMF[vm2l];
        }
        VDELF[vdel31Base + iv] -= vtmp * VDELF[vdel21Base + iv];
        VDELF[vdel32Base + iv] -= vtmp * VDELF[vdel22Base + iv];
      }

      pivot = 1.0 / VMF[vmBase3 + iv * vmStride3 + iv];
      for (int l = ivp; l <= nsys; l += 1) {
        VMF[vmBase3 + iv * vmStride3 + l] *= pivot;
      }
      VDELF[vdel31Base + iv] *= pivot;
      VDELF[vdel32Base + iv] *= pivot;

      {
        const double vtmp1 = VMF[vmBase1 + iv * vmStride3 + iv];
        const double vtmp2 = VMF[vmBase2 + iv * vmStride3 + iv];
        for (int l = ivp; l <= nsys; l += 1) {
          const int vm1l = vmBase1 + iv * vmStride3 + l;
          const int vm2l = vmBase2 + iv * vmStride3 + l;
          const int vm3l = vmBase3 + iv * vmStride3 + l;
          VMF[vm1l] -= vtmp1 * VMF[vm3l];
          VMF[vm2l] -= vtmp2 * VMF[vm3l];
        }
        VDELF[vdel11Base + iv] -= vtmp1 * VDELF[vdel31Base + iv];
        VDELF[vdel21Base + iv] -= vtmp2 * VDELF[vdel31Base + iv];
        VDELF[vdel12Base + iv] -= vtmp1 * VDELF[vdel32Base + iv];
        VDELF[vdel22Base + iv] -= vtmp2 * VDELF[vdel32Base + iv];
      }

      {
        const double vtmp = VAF[va12Base + iv];
        for (int l = ivp; l <= nsys; l += 1) {
          const int vm1l = vmBase1 + iv * vmStride3 + l;
          const int vm2l = vmBase2 + iv * vmStride3 + l;
          VMF[vm1l] -= vtmp * VMF[vm2l];
        }
        VDELF[vdel11Base + iv] -= vtmp * VDELF[vdel21Base + iv];
        VDELF[vdel12Base + iv] -= vtmp * VDELF[vdel22Base + iv];
      }

      if (iv == nsys) continue;

      for (int k = 1; k <= 3; k += 1) {
        const int vbBase = k == 1 ? vbBase1 : (k == 2 ? vbBase2 : vbBase3);
        const int vdelBase = k == 1 ? vdelBase1 : (k == 2 ? vdelBase2 : vdelBase3);
        const int vmBase = k == 1 ? vmBase1 : (k == 2 ? vmBase2 : vmBase3);
        const double vtmp1 = VBF[vbBase + vStride3 + ivp];
        const double vtmp2 = VBF[vbBase + 2 * vStride3 + ivp];
        const double vtmp3 = VMF[vmBase + ivp * vmStride3 + iv];
        for (int l = ivp; l <= nsys; l += 1) {
          const int vmkl = vmBase + ivp * vmStride3 + l;
          const int vm1l = vmBase1 + iv * vmStride3 + l;
          const int vm2l = vmBase2 + iv * vmStride3 + l;
          const int vm3l = vmBase3 + iv * vmStride3 + l;
          VMF[vmkl] -= vtmp1 * VMF[vm1l]
            + vtmp2 * VMF[vm2l]
            + vtmp3 * VMF[vm3l];
        }
        VDELF[vdelBase + vStride3 + ivp] -= vtmp1 * VDELF[vdel11Base + iv]
          + vtmp2 * VDELF[vdel21Base + iv]
          + vtmp3 * VDELF[vdel31Base + iv];
        VDELF[vdelBase + 2 * vStride3 + ivp] -= vtmp1 * VDELF[vdel12Base + iv]
          + vtmp2 * VDELF[vdel22Base + iv]
          + vtmp3 * VDELF[vdel32Base + iv];
      }

      if (iv == ivte1) {
        for (int k = 1; k <= 3; k += 1) {
          const double vtmp1 = VZ[k * 3 + 1];
          const double vtmp2 = VZ[k * 3 + 2];
          const int vdelBase = k == 1 ? vdelBase1 : (k == 2 ? vdelBase2 : vdelBase3);
          const int vmBase = k == 1 ? vmBase1 : (k == 2 ? vmBase2 : vmBase3);
          for (int l = ivp; l <= nsys; l += 1) {
            const int vmkl = vmBase + ivz * vmStride3 + l;
            const int vm1l = vmBase1 + iv * vmStride3 + l;
            const int vm2l = vmBase2 + iv * vmStride3 + l;
            VMF[vmkl] -= vtmp1 * VMF[vm1l]
              + vtmp2 * VMF[vm2l];
          }
          VDELF[vdelBase + vStride3 + ivz] -= vtmp1 * VDELF[vdel11Base + iv]
            + vtmp2 * VDELF[vdel21Base + iv];
          VDELF[vdelBase + 2 * vStride3 + ivz] -= vtmp1 * VDELF[vdel12Base + iv]
            + vtmp2 * VDELF[vdel22Base + iv];
        }
      }

      if (ivp == nsys) continue;

      for (int kv = iv + 2; kv <= nsys; kv += 1) {
        const double vtmp1 = VMF[vmBase1 + kv * vmStride3 + iv];
        const double vtmp2 = VMF[vmBase2 + kv * vmStride3 + iv];
        const double vtmp3 = VMF[vmBase3 + kv * vmStride3 + iv];

        const bool active1 = __builtin_fabs(vtmp1) > vacc1;
        const bool active2 = __builtin_fabs(vtmp2) > vacc2;
        const bool active3 = __builtin_fabs(vtmp3) > vacc3;
        // These destination rows are independent and use the same source row.
        // Fuse the common two/three-row updates without reordering arithmetic
        // within any coefficient or changing the acceleration thresholds.
        if (active2 && active3) {
          if (active1) {
            subtractScaledRows3(VMF + (vmBase3 + iv * vmStride3 + ivp),
              VMF + (vmBase1 + kv * vmStride3 + ivp),
              VMF + (vmBase2 + kv * vmStride3 + ivp),
              VMF + (vmBase3 + kv * vmStride3 + ivp),
              vtmp1, vtmp2, vtmp3, nsys - iv);
            VDELF[vdel11Base + kv] -= vtmp1 * VDELF[vdel31Base + iv];
            VDELF[vdel12Base + kv] -= vtmp1 * VDELF[vdel32Base + iv];
          } else {
            subtractScaledRows2(VMF + (vmBase3 + iv * vmStride3 + ivp),
              VMF + (vmBase2 + kv * vmStride3 + ivp),
              VMF + (vmBase3 + kv * vmStride3 + ivp),
              vtmp2, vtmp3, nsys - iv);
          }
          VDELF[vdel21Base + kv] -= vtmp2 * VDELF[vdel31Base + iv];
          VDELF[vdel22Base + kv] -= vtmp2 * VDELF[vdel32Base + iv];
          VDELF[vdel31Base + kv] -= vtmp3 * VDELF[vdel31Base + iv];
          VDELF[vdel32Base + kv] -= vtmp3 * VDELF[vdel32Base + iv];
          continue;
        }

        if (active1) {
          subtractScaledRows1(VMF + (vmBase3 + iv * vmStride3 + ivp),
            VMF + (vmBase1 + kv * vmStride3 + ivp),
            vtmp1, nsys - iv);
          VDELF[vdel11Base + kv] -= vtmp1 * VDELF[vdel31Base + iv];
          VDELF[vdel12Base + kv] -= vtmp1 * VDELF[vdel32Base + iv];
        }
        if (active2) {
          subtractScaledRows1(VMF + (vmBase3 + iv * vmStride3 + ivp),
            VMF + (vmBase2 + kv * vmStride3 + ivp),
            vtmp2, nsys - iv);
          VDELF[vdel21Base + kv] -= vtmp2 * VDELF[vdel31Base + iv];
          VDELF[vdel22Base + kv] -= vtmp2 * VDELF[vdel32Base + iv];
        }
        if (active3) {
          subtractScaledRows1(VMF + (vmBase3 + iv * vmStride3 + ivp),
            VMF + (vmBase3 + kv * vmStride3 + ivp),
            vtmp3, nsys - iv);
          VDELF[vdel31Base + kv] -= vtmp3 * VDELF[vdel31Base + iv];
          VDELF[vdel32Base + kv] -= vtmp3 * VDELF[vdel32Base + iv];
        }
      }
    }

    for (int iv = nsys; iv >= 2; iv -= 1) {
      double vtmp = VDELF[vdel31Base + iv];
      for (int kv = iv - 1; kv >= 1; kv -= 1) {
        VDELF[vdel11Base + kv] -= VMF[vmBase1 + kv * vmStride3 + iv] * vtmp;
        VDELF[vdel21Base + kv] -= VMF[vmBase2 + kv * vmStride3 + iv] * vtmp;
        VDELF[vdel31Base + kv] -= VMF[vmBase3 + kv * vmStride3 + iv] * vtmp;
      }
      vtmp = VDELF[vdel32Base + iv];
      for (int kv = iv - 1; kv >= 1; kv -= 1) {
        VDELF[vdel12Base + kv] -= VMF[vmBase1 + kv * vmStride3 + iv] * vtmp;
        VDELF[vdel22Base + kv] -= VMF[vmBase2 + kv * vmStride3 + iv] * vtmp;
        VDELF[vdel32Base + kv] -= VMF[vmBase3 + kv * vmStride3 + iv] * vtmp;
      }
    }

    return true;
}

static bool finiteArray(const double* values, int count) {
  for (int i = 0; i < count; ++i) if (!__builtin_isfinite(values[i])) return false;
  return true;
}

extern "C" int blsolv_flat(int nsys, int ivte1, int ivz, double vaccel, double s1, double sN,
                           double* va, const double* vb, double* vm, double* vdel, const double* vz) {
  transposeMassCoefficients(vm, nsys);
  const bool ok = solve(nsys, ivte1, ivz, vaccel, s1, sN, va, vb, vm, vdel, vz);
  transposeMassCoefficients(vm, nsys);
  if (!ok) return 0;
  const int smallSize = 12 * (nsys + 1);
  const int massSize = 4 * (nsys + 1) * (nsys + 1);
  return finiteArray(va, smallSize) && finiteArray(vm, massSize) && finiteArray(vdel, smallSize);
}
