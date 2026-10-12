// Optional acceleration for the flat, double-precision BLSOLV path.
// The JavaScript caller remains the reference and fallback implementation.
import { blsolvBase64 } from './blsolv_wasm_data.js';

const set = Float64Array.prototype.set;

function nativeView(value, prototype) {
  return ArrayBuffer.isView(value) && Object.getPrototypeOf(value) === prototype
    && value.buffer instanceof ArrayBuffer;
}

// One scratch memory per solver worker. No context or caller-owned views are
// retained. Returning false always leaves the caller's arrays untouched.
export function createWasmBlsolv() {
  let instance;
  try {
    const decoded = atob(blsolvBase64);
    const bytes = Uint8Array.from(decoded, char => char.charCodeAt(0));
    instance = new WebAssembly.Instance(new WebAssembly.Module(bytes));
  } catch {
    // WebAssembly/SIMD may be unavailable or blocked by the host's CSP.
    return null;
  }
  const { memory, blsolv_flat: solve, __heap_base: heapBase } = instance.exports;
  let size = -1;
  let offsets;
  let views;

  return function trySolve(ctx) {
    const n = ctx.NSYS;
    if (!Number.isInteger(n) || n < 1 || n > 1024) return false;
    const small = 12 * (n + 1);
    const mass = 4 * (n + 1) * (n + 1);
    const arrays = [ctx.VAF, ctx.VBF, ctx.VMF, ctx.VDELF];
    if (arrays.some((value, i) => !nativeView(value, Float64Array.prototype)
      || value.length !== (i === 2 ? mass : small))) return false;
    const buffers = arrays.map(value => value.buffer);
    if (new Set(buffers).size !== 4) return false;

    // Snapshotting metadata is valid only if elimination cannot modify it.
    // Shared memory and unusual storage layouts use the original routine.
    const independent = value => (nativeView(value, Float64Array.prototype)
      || nativeView(value, Int32Array.prototype)) && !buffers.includes(value.buffer);
    const { S, N, IBLTE, ISYS, VZ, VACCEL } = ctx;
    if (!independent(S) || !Number.isInteger(N) || N < 1 || N >= S.length
      || !Number.isFinite(S[1]) || !Number.isFinite(S[N])
      || !Number.isFinite(VACCEL) || !independent(IBLTE)
      || !Array.isArray(ISYS) || !Array.isArray(VZ)) return false;
    const edgeRow = ISYS[IBLTE[1]];
    if (!independent(edgeRow)) return false;
    const ivte1 = edgeRow[1];
    if (!Number.isInteger(ivte1) || ivte1 < 1 || ivte1 > n) return false;
    let ivz = 0;
    if (ivte1 < n) {
      const wakeRow = ISYS[IBLTE[2] + 1];
      if (!independent(wakeRow)) return false;
      ivz = wakeRow[2];
      if (!Number.isInteger(ivz) || ivz < 1 || ivz > n) return false;
      for (let k = 1; k <= 3; k += 1) {
        if (!independent(VZ[k]) || VZ[k].length < 3) return false;
      }
    }

    try {
      if (size !== n) {
        let offset = (heapBase.value + 7) & ~7;
        const lengths = [small, small, mass, small, 12];
        const nextOffsets = [];
        for (const length of lengths) {
          nextOffsets.push(offset);
          offset += length * 8;
        }
        if (offset > memory.buffer.byteLength) {
          memory.grow(Math.ceil((offset - memory.buffer.byteLength) / 65536));
        }
        views = lengths.map((length, i) => new Float64Array(memory.buffer, nextOffsets[i], length));
        offsets = nextOffsets;
        size = n;
      }
      for (let i = 0; i < 4; i += 1) set.call(views[i], arrays[i]);
      if (ivte1 < n) {
        for (let k = 1; k <= 3; k += 1) {
          views[4][3 * k + 1] = VZ[k][1];
          views[4][3 * k + 2] = VZ[k][2];
        }
      }
      // A singular pivot or any nonfinite output discards the scratch result.
      // Running JS from the original input preserves warnings, partial writes,
      // and NaN payloads as well as the return value on these uncommon paths.
      if (!solve(n, ivte1, ivz, VACCEL, S[1], S[N], ...offsets)) return false;
    } catch {
      // A failed growth/instantiation must not leave stale detached views.
      size = -1;
      return false;
    }
    set.call(arrays[0], views[0]);
    set.call(arrays[2], views[2]);
    set.call(arrays[3], views[3]);
    return true;
  };
}

let solver;
export function tryWasmBlsolv(ctx) {
  // Small systems do too little work to justify copying into scratch memory.
  // Inviscid calculations never instantiate the module.
  if (!(ctx.NSYS >= 32)) return false;
  if (solver === undefined) solver = createWasmBlsolv();
  return solver ? solver(ctx) : false;
}
