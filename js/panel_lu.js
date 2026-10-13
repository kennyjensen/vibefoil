import { ludcmp } from './xsolve.js';
import { tryWasmPanelLu } from './panel_wasm.js';

// Retain one panel matrix's factorization. Compare the assembled matrix itself:
// geometry names and angles alone cannot establish that it is safe to reuse.
export function createCachedLudcmp() {
  let size = 0;
  let valid = false;
  let input;
  let factors;
  let pivots;

  // Return whether the factorization was reused (callers need not use this).
  return function cachedLudcmp(n, a, indx) {
    // The worker uses independent Float64 rows and Int32 pivots. Keep other
    // storage types, overlapping buffers, and the original overflow handling
    // on the unmodified path, including their rounding/aliasing behavior.
    let eligible = Number.isInteger(n) && n > 0 && n <= 515
      && indx instanceof Int32Array && indx.length >= n
      && indx.buffer instanceof ArrayBuffer;
    const buffers = new Set();
    if (eligible) {
      buffers.add(indx.buffer);
      for (let i = 0; i < n; i += 1) {
        const row = a[i];
        if (!(row instanceof Float64Array) || row.length < n
          || !(row.buffer instanceof ArrayBuffer) || buffers.has(row.buffer)) {
          eligible = false;
          break;
        }
        buffers.add(row.buffer);
      }
    }
    if (!eligible) {
      valid = false;
      ludcmp(n, a, indx);
      return false;
    }

    let matches = valid && n === size;
    for (let i = 0; matches && i < n; i += 1) {
      const row = a[i];
      const offset = i * n;
      for (let j = 0; j < n; j += 1) {
        if (!Object.is(row[j], input[offset + j])) {
          matches = false;
          break;
        }
      }
    }
    if (matches) {
      // LUDCMP exchanges whole row objects. Replay those exchanges so unused
      // columns and references to individual rows retain the same behavior.
      for (let j = 0; j < n; j += 1) {
        const pivot = pivots[j];
        if (j !== pivot) {
          const row = a[pivot];
          a[pivot] = a[j];
          a[j] = row;
        }
        indx[j] = pivot;
      }
      for (let i = 0; i < n; i += 1) {
        a[i].set(factors.subarray(i * n, (i + 1) * n));
      }
      return true;
    }

    valid = false;
    if (size !== n) {
      size = n;
      input = new Float64Array(n * n);
      factors = new Float64Array(n * n);
      pivots = new Int32Array(n);
    }
    let finite = true;
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) {
        const value = a[i][j];
        input[i * n + j] = value;
        if (!Number.isFinite(value)) finite = false;
      }
    }
    if (!tryWasmPanelLu(n, a, indx)) ludcmp(n, a, indx);
    if (!finite) return false;
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) {
        const value = a[i][j];
        if (!Number.isFinite(value)) return false;
        factors[i * n + j] = value;
      }
    }
    pivots.set(indx.subarray(0, n));
    valid = true;
    return false;
  };
}
