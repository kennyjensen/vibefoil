// Preserve the row API and Float64 values while using one owned buffer per
// snapshot matrix. The solver's original rows remain independent of the copy.
export function cloneFloat64Rows(rows) {
  if (!rows) return null;
  let length = 0;
  for (const row of rows) {
    if (row instanceof Float64Array) length += row.length;
  }
  const data = new Float64Array(length);
  const copy = new Array(rows.length);
  let offset = 0;
  for (let i = 0; i < rows.length; i += 1) {
    const row = rows[i];
    if (row instanceof Float64Array) {
      copy[i] = data.subarray(offset, offset + row.length);
      copy[i].set(row);
      offset += row.length;
    } else {
      copy[i] = row;
    }
  }
  return copy;
}

// Only these fields own disposable output buffers. Other fields can reference
// live solver arrays (notably coeffs.ACRIT), which must never be detached.
export function resultTransfers(result) {
  const buffers = new Set();
  const add = value => {
    if (ArrayBuffer.isView(value) && value.buffer instanceof ArrayBuffer) {
      buffers.add(value.buffer);
    }
  };
  add(result.xb);
  add(result.yb);
  add(result.streamlines?.grid);
  for (const snapshot of [result.reuseState?.panel, result.reuseState?.bl]) {
    if (!snapshot) continue;
    for (const value of Object.values(snapshot)) {
      if (Array.isArray(value)) {
        // cloneFloat64Rows leaves other row types untouched.
        for (const row of value) {
          if (row instanceof Float64Array) add(row);
        }
      } else {
        add(value);
      }
    }
  }
  return [...buffers];
}
