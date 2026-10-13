// GPL-2.0-or-later. Optional panel linear algebra; JS remains the fallback.
// Kernels only touch scratch memory. Rejecting a call changes no input bytes.
import { panelBase64 } from './panel_wasm_data.js';

const set = Float64Array.prototype.set;
const subarray = Float64Array.prototype.subarray;
const ownData = (object, key, writable = false) => {
  const descriptor = Object.getOwnPropertyDescriptor(object, key);
  if (!descriptor || !('value' in descriptor) || writable && !descriptor.writable) throw REJECT;
  return descriptor.value;
};
const REJECT = new Error('Unsupported panel storage');
function native(value, prototype) {
  if (!ArrayBuffer.isView(value) || Object.getPrototypeOf(value) !== prototype) return false;
  for (const key of ['buffer', 'length', 'byteOffset', 'byteLength']) if (Object.hasOwn(value, key)) return false;
  return Object.getPrototypeOf(value.buffer) === ArrayBuffer.prototype && !value.buffer.resizable;
}
function rows(matrix, n, width, buffers, writable = false) {
  if (!Array.isArray(matrix) || Object.getPrototypeOf(matrix) !== Array.prototype || matrix.length < n) throw REJECT;
  const result = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = ownData(matrix, String(i), writable);
    if (!native(row, Float64Array.prototype) || row.length < width || buffers.has(row.buffer)) throw REJECT;
    buffers.add(row.buffer);
    result[i] = row;
  }
  return result;
}

export function createWasmPanel() {
  let instance;
  try {
    const bytes = Uint8Array.from(atob(panelBase64), c => c.charCodeAt(0));
    instance = new WebAssembly.Instance(new WebAssembly.Module(bytes));
  } catch { return null; }
  const { memory, __heap_base: heapBase, panel_lu: lu, panel_solve: solve } = instance.exports;
  const stats = {factorizations: 0, solves: 0, fallback: 0};
  function prepare(n, a, indx, extra, writable) {
    if (!Number.isInteger(n) || n < 1 || n > 515
      || !native(indx, Int32Array.prototype) || indx.length < n) throw REJECT;
    const buffers = new Set([indx.buffer]);
    const matrix = rows(a, n, n, buffers, writable);
    const offset = (heapBase.value + 7) & ~7;
    const pivotOffset = offset + n*n*8;
    const extraOffset = (pivotOffset + n*4 + 7) & ~7;
    const required = extraOffset + extra*8;
    if (required > memory.buffer.byteLength) memory.grow(Math.ceil((required - memory.buffer.byteLength)/65536));
    const packed = new Float64Array(memory.buffer, offset, n*n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const value = matrix[i][j];
        if (!Number.isFinite(value)) throw REJECT;
        packed[i*n+j] = value;
      }
    }
    return {matrix, buffers, offset, pivotOffset, extraOffset, packed,
      pivots: new Int32Array(memory.buffer, pivotOffset, n)};
  }
  return {
    stats,
    factorize(n, a, indx) {
      let work;
      try {
        work = prepare(n, a, indx, n, true);
        if (!lu(n, work.offset, work.pivotOffset, work.extraOffset)) throw REJECT;
      } catch { stats.fallback++;return false; }
      // Replay swaps of whole JS row objects, including their unused padding.
      for (let j = 0; j < n; j++) {
        const pivot = work.pivots[j];
        if (j !== pivot) { const row = a[pivot];a[pivot] = a[j];a[j] = row; }
        indx[j] = pivot;
      }
      for (let i = 0; i < n; i++) set.call(a[i], subarray.call(work.packed, i*n, (i+1)*n));
      stats.factorizations++;
      return true;
    },
    // Read a range of RHS columns from a row matrix; return owned, column-major
    // solutions. The caller chooses which result rows to publish.
    solveColumns(n, a, indx, b, first, count) {
      try {
        if (!Number.isInteger(first) || first < 0 || !Number.isInteger(count) || count < 1 || count > 1024
          || first + count > 2048) throw REJECT;
        const paddedCount = (count + 3) & ~3;
        const work = prepare(n, a, indx, n*paddedCount, false);
        const rhs = rows(b, n, first+count, work.buffers);
        const columns = new Float64Array(memory.buffer, work.extraOffset, n*paddedCount);
        columns.fill(0);
        for (let i = 0; i < n; i++) {
          const pivot = indx[i];
          if (pivot < 0 || pivot >= n) throw REJECT;
          work.pivots[i] = pivot;
          for (let k = 0; k < count; k++) {
            const value = rhs[i][first+k];
            if (!Number.isFinite(value)) throw REJECT;
            columns[((k >> 2)*n+i)*4+(k & 3)] = value;
          }
        }
        if (!solve(n, work.offset, work.pivotOffset, work.extraOffset, paddedCount)) throw REJECT;
        const result = new Float64Array(n*count);
        for (let k = 0; k < count; k++) {
          for (let i = 0; i < n; i++) result[k*n+i] = columns[((k >> 2)*n+i)*4+(k & 3)];
        }
        stats.solves++;
        return result;
      } catch { stats.fallback++;return null; }
    },
  };
}

let accelerator;
function getAccelerator(n) {
  if (!(n >= 32)) return null;
  if (accelerator === undefined) accelerator = createWasmPanel();
  return accelerator;
}
export function tryWasmPanelLu(n, a, indx) {
  return getAccelerator(n)?.factorize(n, a, indx) ?? false;
}
export function tryWasmPanelColumns(n, a, indx, b, first, count) {
  return getAccelerator(n)?.solveColumns(n, a, indx, b, first, count) ?? null;
}
export function getPanelWasmStats() {
  return accelerator ? {...accelerator.stats} : {factorizations: 0, solves: 0, fallback: 0};
}
