// Shared array helpers for XFOIL port data structures.

function createMatrix(rows, cols) {
  const mat = new Array(rows);
  for (let i = 0; i < rows; i += 1) {
    mat[i] = new Float64Array(cols);
  }
  return mat;
}

function createMatrix1(rows, cols) {
  const mat = new Array(rows + 1);
  for (let i = 0; i <= rows; i += 1) {
    mat[i] = new Float64Array(cols + 1);
  }
  return mat;
}

function createTensor3(d1, d2, d3) {
  const arr = new Array(d1 + 1);
  for (let i = 0; i <= d1; i += 1) {
    arr[i] = new Array(d2 + 1);
    for (let j = 0; j <= d2; j += 1) {
      arr[i][j] = new Float64Array(d3 + 1);
    }
  }
  return arr;
}

function createTensor3Flat(d1, d2, d3) {
  const stride3 = d3 + 1;
  const stride2 = (d2 + 1) * stride3;
  const flat = new Float64Array((d1 + 1) * stride2);
  const view = new Array(d1 + 1);
  for (let i = 0; i <= d1; i += 1) {
    view[i] = new Array(d2 + 1);
    for (let j = 0; j <= d2; j += 1) {
      const offset = i * stride2 + j * stride3;
      view[i][j] = flat.subarray(offset, offset + stride3);
    }
  }
  return { flat, view };
}

export { createMatrix, createMatrix1, createTensor3, createTensor3Flat };
