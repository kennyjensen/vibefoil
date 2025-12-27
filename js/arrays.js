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

export { createMatrix, createMatrix1, createTensor3 };
