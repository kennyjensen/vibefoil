// GPL-2.0-or-later. XFOIL LUDCMP and BAKSUB, matching js/xsolve.js.
// Strict Float64, original pivot ties and sequential dot-product reductions.
// Caller supplies private scratch arrays. No fast math or fused multiply-add.
#include <wasm_simd128.h>
static bool finite(double x) { return __builtin_isfinite(x); }

extern "C" int panel_lu(int n, double* a, int* pivots, double* scale) {
  for (int i = 0; i < n; ++i) {
    double maximum = 0.0;
    for (int j = 0; j < n; ++j) {
      const double value = __builtin_fabs(a[i*n+j]);
      if (value > maximum) maximum = value;
    }
    scale[i] = 1.0 / maximum;
    if (!finite(scale[i])) return 0;
  }
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      double sum = a[i*n+j];
      for (int k = 0; k < i; ++k) sum -= a[i*n+k] * a[k*n+j];
      a[i*n+j] = sum;
    }
    double maximum = 0.0;
    int imax = j;
    for (int i = j; i < n; ++i) {
      double sum = a[i*n+j];
      for (int k = 0; k < j; ++k) sum -= a[i*n+k] * a[k*n+j];
      a[i*n+j] = sum;
      const double value = scale[i] * __builtin_fabs(sum);
      if (value >= maximum) { imax = i; maximum = value; }
    }
    if (j != imax) {
      for (int k = 0; k < n; ++k) {
        const double value = a[imax*n+k];
        a[imax*n+k] = a[j*n+k];
        a[j*n+k] = value;
      }
      scale[imax] = scale[j];
    }
    pivots[j] = imax;
    if (j != n-1) {
      const double reciprocal = 1.0 / a[j*n+j];
      for (int i = j+1; i < n; ++i) a[i*n+j] *= reciprocal;
    }
  }
  for (int i = 0; i < n*n; ++i) if (!finite(a[i])) return 0;
  return 1;
}

extern "C" int panel_solve(int n, const double* a, const int* pivots,
                            double* columns, int count) {
  // Independent columns occupy SIMD lanes. Each lane retains the original
  // sequential subtraction order; no horizontal reductions or relaxed SIMD.
  for (int column = 0; column < count; column += 4) {
    double* b = columns + column*n;
    int ii0 = -1, ii1 = -1, ii2 = -1, ii3 = -1;
    for (int i = 0; i < n; ++i) {
      const int ll = pivots[i];
      double s0 = b[ll*4], s1 = b[ll*4+1], s2 = b[ll*4+2], s3 = b[ll*4+3];
      wasm_v128_store(b+ll*4, wasm_v128_load(b+i*4));
      wasm_v128_store(b+ll*4+2, wasm_v128_load(b+i*4+2));
      if (ii0 == -1 && s0 != 0.0) ii0 = i;
      if (ii1 == -1 && s1 != 0.0) ii1 = i;
      if (ii2 == -1 && s2 != 0.0) ii2 = i;
      if (ii3 == -1 && s3 != 0.0) ii3 = i;
      int start = ii0 < 0 ? i : ii0;
      if ((ii1 < 0 ? i : ii1) > start) start = ii1 < 0 ? i : ii1;
      if ((ii2 < 0 ? i : ii2) > start) start = ii2 < 0 ? i : ii2;
      if ((ii3 < 0 ? i : ii3) > start) start = ii3 < 0 ? i : ii3;
      if (ii0 >= 0) for (int j = ii0; j < start; ++j) s0 -= a[i*n+j] * b[j*4];
      if (ii1 >= 0) for (int j = ii1; j < start; ++j) s1 -= a[i*n+j] * b[j*4+1];
      if (ii2 >= 0) for (int j = ii2; j < start; ++j) s2 -= a[i*n+j] * b[j*4+2];
      if (ii3 >= 0) for (int j = ii3; j < start; ++j) s3 -= a[i*n+j] * b[j*4+3];
      v128_t s01 = wasm_f64x2_make(s0, s1), s23 = wasm_f64x2_make(s2, s3);
      for (int j = start; j < i; ++j) {
        const v128_t coefficient = wasm_f64x2_splat(a[i*n+j]);
        s01 = wasm_f64x2_sub(s01, wasm_f64x2_mul(coefficient, wasm_v128_load(b+j*4)));
        s23 = wasm_f64x2_sub(s23, wasm_f64x2_mul(coefficient, wasm_v128_load(b+j*4+2)));
      }
      wasm_v128_store(b+i*4, s01);
      wasm_v128_store(b+i*4+2, s23);
    }
    for (int i = n-1; i >= 0; --i) {
      v128_t s01 = wasm_v128_load(b+i*4), s23 = wasm_v128_load(b+i*4+2);
      for (int j = i+1; j < n; ++j) {
        const v128_t coefficient = wasm_f64x2_splat(a[i*n+j]);
        s01 = wasm_f64x2_sub(s01, wasm_f64x2_mul(coefficient, wasm_v128_load(b+j*4)));
        s23 = wasm_f64x2_sub(s23, wasm_f64x2_mul(coefficient, wasm_v128_load(b+j*4+2)));
      }
      const v128_t diagonal = wasm_f64x2_splat(a[i*n+i]);
      wasm_v128_store(b+i*4, wasm_f64x2_div(s01, diagonal));
      wasm_v128_store(b+i*4+2, wasm_f64x2_div(s23, diagonal));
    }
  }
  for (int i = 0; i < n*count; ++i) if (!finite(columns[i])) return 0;
  return 1;
}
