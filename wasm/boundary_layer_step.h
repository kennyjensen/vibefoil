// Share SETBL scratch storage with the existing strict-double BLSOLV kernel.
// Derived from XFOIL by Mark Drela. GNU General Public License v2 or later.
extern "C" int blsolv_flat(int, int, int, double, double, double,
                           double*, const double*, double*, double*, const double*);

V1 flat_tensor(V3 tensor, V1 flat, u32 middle, u32 stride) {
  if (tensor.length != 4 || flat.length != 4 * middle * stride) unsupported();
  for (u32 k = 0; k < 4; k++) {
    if (tensor[k].length != middle) unsupported();
    for (u32 j = 0; j < middle; j++) {
      auto row = tensor[k][j];
      if (row.length != stride || row.pointer != flat.pointer + 8 * (k * middle + j) * stride) unsupported();
    }
  }
  return flat;
}
bool overlaps(V1 a, V1 b) {
  return a.pointer < b.pointer + 8 * b.length && b.pointer < a.pointer + 8 * a.length;
}
void solve_boundary_layer_matrix(Context ctx) {
  u32 n = size_value(ctx.get(scalar_NSYS));
  if (n < 32 || n > 1024) unsupported();
  u32 stride = n + 1;
  V1 matrices[] = {
    flat_tensor(ctx.arr<3>(array_VA), ctx.arr<1>(array_VAF), 3, stride),
    flat_tensor(ctx.arr<3>(array_VB), ctx.arr<1>(array_VBF), 3, stride),
    flat_tensor(ctx.arr<3>(array_VM), ctx.arr<1>(array_VMF), stride, stride),
    flat_tensor(ctx.arr<3>(array_VDEL), ctx.arr<1>(array_VDELF), 3, stride),
  };
  for (u32 i = 0; i < 4; i++) for (u32 j = i + 1; j < 4; j++) {
    if (overlaps(matrices[i], matrices[j])) unsupported();
  }
  auto independent = [&](V1 a) {
    for (auto matrix : matrices) if (overlaps(a, matrix)) unsupported();
  };
  V1 s = ctx.arr<1>(array_S), iblte = ctx.arr<1>(array_IBLTE);
  independent(s); independent(iblte);
  u32 count = size_value(ctx.get(scalar_N));
  if (count < 1 || count >= s.length) unsupported();
  V1 edge = ctx.arr<2>(array_ISYS)[iblte[1]];
  independent(edge);
  u32 ivte1 = size_value(edge[1]);
  if (ivte1 < 1 || ivte1 > n) unsupported();
  u32 ivz = 0;
  auto vz = make_vec<1>(12);
  if (ivte1 < n) {
    V1 wake = ctx.arr<2>(array_ISYS)[iblte[2] + 1];
    independent(wake);
    ivz = size_value(wake[2]);
    if (ivz < 1 || ivz > n) unsupported();
    for (u32 k = 1; k <= 3; k++) {
      V1 row = ctx.arr<2>(array_VZ)[k];
      independent(row);
      vz[3 * k + 1] = row[1]; vz[3 * k + 2] = row[2];
    }
  }
  double vaccel = ctx.get(scalar_VACCEL);
  if (!finite(vaccel) || !finite(s[1]) || !finite(s[count])) unsupported();
  if (!blsolv_flat(n, ivte1, ivz, vaccel, s[1], s[count],
      (double*)matrices[0].pointer, (double*)matrices[1].pointer,
      (double*)matrices[2].pointer, (double*)matrices[3].pointer, (double*)vz.pointer)) unsupported();
}
