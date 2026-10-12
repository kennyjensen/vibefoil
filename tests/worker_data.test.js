import assert from 'node:assert/strict';
import { cloneFloat64Rows, resultTransfers } from '../js/worker_data.js';

// Ragged rows and special Float64 values survive packing without aliasing the
// original matrix or neighboring rows. Non-Float64 entries retain old behavior.
const liveIndices = Int32Array.of(2, 3);
const input = [Float64Array.of(1, -0, NaN), null, Float64Array.of(Infinity, -3),
  new Float64Array(0), undefined, [7], liveIndices];
const copy = cloneFloat64Rows(input);
assert.deepStrictEqual(copy, input);
assert.notEqual(copy, input);
assert.equal(copy[0].buffer, copy[2].buffer);
assert.notEqual(copy[0].buffer, input[0].buffer);
copy[0][0] = 42;
assert.equal(input[0][0], 1);
assert.equal(copy[2][0], Infinity);
input[2][1] = 9;
assert.equal(copy[2][1], -3);
assert.equal(cloneFloat64Rows(null), null);
assert.deepEqual(cloneFloat64Rows([]), []);

// Transfer real buffers, checking both the received values and ownership on
// the sender. Sharing a packed buffer must not duplicate it in the transfer list.
const liveCrit = Float64Array.of(9, 9, 9);
const message = {
  xb: Float64Array.of(1, 0, 1), yb: Float64Array.of(0, 0.1, 0),
  streamlines: { grid: Float64Array.of(NaN, -0, 0.3) },
  coeffs: { ACRIT: liveCrit },
  reuseState: {
    panel: { AIJ: copy, X: Float64Array.of(1, 0, 1), AIJPIV: Int32Array.of(1, 2, 0) },
    bl: { ACRIT: new Float64Array(liveCrit), TFORCE: [false, true, false] },
  },
};
const expected = structuredClone(message);
const transfers = resultTransfers(message);
assert.equal(transfers.length, 7);
assert.equal(new Set(transfers).size, transfers.length);
assert.ok(!transfers.includes(liveCrit.buffer));
assert.ok(!transfers.includes(liveIndices.buffer));
const received = structuredClone(message, { transfer: transfers });
assert.deepStrictEqual(received, expected);
assert.ok(transfers.every(buffer => buffer.byteLength === 0));
assert.deepStrictEqual(liveCrit, Float64Array.of(9, 9, 9));
assert.deepStrictEqual(liveIndices, Int32Array.of(2, 3));
assert.equal(received.reuseState.panel.AIJ[0].buffer, received.reuseState.panel.AIJ[2].buffer);
assert.deepEqual(resultTransfers({ ok: false }), []);
assert.equal(resultTransfers({ xb: Float64Array.of(1), yb: Float64Array.of(2) }).length, 2);
console.log('worker_data.test.js: OK (exact values and buffer ownership)');
