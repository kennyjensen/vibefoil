// GPL-2.0-or-later. Optional SETBL/marching/closure and BLSOLV acceleration. All computation
// takes place in scratch memory; failure leaves the original JS context intact.
import { boundaryLayerBase64 } from './boundary_layer_wasm_data.js';
import { layout, logFormats } from './boundary_layer_layout.js';
import * as equations from './xblsys.js';

const arrayEntries = Object.entries(layout.arrays);
const sourceKeys = ['panel', 'system', 'sign', 'mass', 'sensitivity'];
const booleanFields = new Set(layout.booleanFields);
const integerArrays = new Set(['NBL', 'IBLTE', 'ITRAN', 'IV_TO_IS', 'IV_TO_IBL']);
// SETBL, its callees and BLSOLV only read these geometry/influence arrays.
// A buffer shared with any writable field is still copied back in full.
const readOnlyArrays = new Set(['DIJ', 'X', 'Y', 'S', 'XP', 'YP']);
const nativeSet = Float64Array.prototype.set;
const integerSet = Int32Array.prototype.set;
const wordSet = Uint32Array.prototype.set;
const lookupGetter = Object.prototype.__lookupGetter__;
const MAX_BYTES = 64 * 1024 * 1024;
const REJECT = new Error('Unsupported boundary-layer storage');
function finiteArray(values) {
  for (let i = 0; i < values.length; i++) if (!Number.isFinite(values[i])) return false;
  return true;
}
// Matrix rows have adjacent pointer/length descriptors. Pack those immutable
// words once per layout so the bridge can use bulk copies on later calls.
function descriptorBlocks(descriptors) {
  const ordered = [...descriptors].sort((a, b) => a.offset - b.offset);
  const blocks = [];
  for (let first = 0; first < ordered.length;) {
    let end = first + 1;
    while (end < ordered.length && ordered[end].offset === ordered[end - 1].offset + 8) end++;
    const values = new Uint32Array((end - first) * 2);
    for (let i = first; i < end; i++) {
      values[(i - first) * 2] = ordered[i].value[0];
      values[(i - first) * 2 + 1] = ordered[i].value[1];
    }
    blocks.push({offset: ordered[first].offset / 4, values});
    first = end;
  }
  return blocks;
}
function sameDescriptors(words, blocks) {
  for (const {offset, values} of blocks) {
    for (let i = 0; i < values.length; i++) if (words[offset + i] !== values[i]) return false;
  }
  return true;
}
const ownData = (object, key) => {
  const descriptor = Object.getOwnPropertyDescriptor(object, key);
  if (descriptor && (!('value' in descriptor) || !descriptor.writable)) throw REJECT;
  return descriptor?.value;
};
const plain = value => value && Object.getPrototypeOf(value) === Object.prototype
  && Object.isExtensible(value);
const native = (value, prototype) => ArrayBuffer.isView(value)
  && Object.getPrototypeOf(value) === prototype && value.buffer instanceof ArrayBuffer
  && !value.buffer.resizable;

export function createWasmBoundaryLayer({onFallback, onTiming} = {}) {
  let logs = [];
  let instance;
  try {
    const binary = Uint8Array.from(atob(boundaryLayerBase64), c => c.charCodeAt(0));
    instance = new WebAssembly.Instance(new WebAssembly.Module(binary), {
      math: { pow: Math.pow, exp: Math.exp, log: Math.log, log10: Math.log10,
        sin: Math.sin, cos: Math.cos, atan: Math.atan, tanh: Math.tanh },
      diagnostics: { log(id, a, b, c) {
        if (!logFormats[id] || logs.length >= 100000) throw REJECT;
        logs.push([id, [a, b, c]]);
      } },
    });
  } catch { return null; }
  const {memory, boundary_layer: assemble, __heap_base: heapBase} = instance.exports;
  const stats = {calls: 0, accelerated: 0, fallback: 0};
  const plans = new WeakMap();
  let residentPlan;
  let running = false;

  function tryAssemble(ctx, solve = false) {
    const timings = onTiming ? [performance.now()] : null;
    stats.calls++;
    if (running) { stats.fallback++; return false; }
    running = true;
    let commits;
    try {
      if (!plain(ctx) || ownData(ctx, '_ENSURED') !== true || ownData(ctx, 'NCOM') !== 73) throw REJECT;
      if (solve) for (const name of ['VAF', 'VBF', 'VMF', 'VDELF']) {
        if (!native(ownData(ctx, name), Float64Array.prototype)) throw REJECT;
      }
      const nsys = ownData(ctx, 'NSYS');
      if (!Number.isInteger(nsys) || nsys < 1 || nsys > 1024) throw REJECT;
      for (const name of layout.callbacks) {
        if (name === 'ensureViscousArrays') continue;
        if (ownData(ctx, name) !== equations[name]) throw REJECT;
      }
      // These helpers are no-ops in the kernel only when their JS preparation
      // contract is already satisfied. Other contexts retain the original path.
      for (const key of ['VA', 'VB', 'VM', 'VDEL', 'DSSR', 'MDEF', 'PDEF', 'KDEF']) {
        if (!ownData(ctx, key)) throw REJECT;
      }

      let cursor = (heapBase.value + 7) & ~7;
      let allocations = [];
      const allocate = bytes => {
        if (!Number.isSafeInteger(bytes) || bytes < 0 || cursor + bytes > MAX_BYTES) throw REJECT;
        const offset = cursor;
        cursor += Math.max(8, (bytes + 7) & ~7);
        return offset;
      };
      const scalarsOffset = allocate(layout.scalars.length * 16);
      const arraysOffset = allocate(arrayEntries.length * 8);
      const sourcesOffset = allocate(layout.sources.length * 4);
      const buffers = new Map();
      let objects = new Map();
      let descriptors = [];
      const links = [];
      const shapes = [];
      const linked = (object, key) => {
        const value = ownData(object, key);
        links.push({object, key, value, row: Array.isArray(object)});
        return value;
      };

      function pack(value, rank, kind = 'float') {
        if (value == null) return [0, 0];
        const seen = objects.get(value);
        if (seen) { if (seen.rank !== rank || seen.kind !== kind) throw REJECT; return seen.descriptor; }
        const prototype = Object.getPrototypeOf(value);
        const isArray = Array.isArray(value);
        if (isArray ? prototype !== Array.prototype || !Object.isExtensible(value)
          : rank !== 1 || !native(value, Float64Array.prototype) && !native(value, Int32Array.prototype)) throw REJECT;
        const length = value.length;
        if (!Number.isInteger(length) || length < 0 || length > 20000000) throw REJECT;
        shapes.push({value, length, prototype, isArray});
        let pointer;
        if (rank === 1 && native(value, Float64Array.prototype)) {
          let base = buffers.get(value.buffer);
          if (base === undefined) {
            if (value.buffer.byteLength % 8) throw REJECT;
            base = allocate(value.buffer.byteLength);
            buffers.set(value.buffer, base);
            const data = new Float64Array(value.buffer);
            if (!finiteArray(data)) throw REJECT;
            allocations.push({offset: base, data, original: data});
          }
          if (base < 0) throw REJECT;
          pointer = base + value.byteOffset;
        } else if (rank === 1 && native(value, Int32Array.prototype)) {
          // Integer storage would round intermediate equation results, not just
          // the final copy. Only genuine index arrays may use it.
          if (kind !== 'int') throw REJECT;
          // Integer views are converted to exact doubles; overlapping integer
          // views or integer/float aliases are deliberately unsupported.
          if (buffers.has(value.buffer)) throw REJECT;
          pointer = allocate(length * 8);
          buffers.set(value.buffer, -1);
          allocations.push({offset: pointer, data: value, original: value});
        } else {
          if (!Array.isArray(value) || Object.getPrototypeOf(value) !== Array.prototype
            || !Object.isExtensible(value)) throw REJECT;
          pointer = allocate(length * 8);
          const record = {rank, kind, descriptor: [pointer, length]};
          objects.set(value, record);
          if (rank === 1) {
            const data = new Float64Array(length);
            for (let i = 0; i < length; i++) {
              const item = ownData(value, String(i));
              if (typeof item !== (kind === 'bool' ? 'boolean' : 'number')) throw REJECT;
              if (!Number.isFinite(Number(item))) throw REJECT;
              data[i] = Number(item);
            }
            allocations.push({offset: pointer, data, original: value, kind});
          } else {
            for (let i = 0; i < length; i++) descriptors.push({offset: pointer + i * 8,
              value: pack(linked(value, String(i)), rank - 1, kind)});
          }
        }
        const descriptor = [pointer, length];
        objects.set(value, {rank, kind, descriptor});
        return descriptor;
      }

      let plan = plans.get(ctx);
      if (plan && (!plan.links.every(link => {
        if (!link.row || !lookupGetter) return ownData(link.object, link.key) === link.value;
        // Reading a cached matrix row needs no freshly allocated descriptor.
        // Reject getters before invoking them; any actual replacement is checked
        // for writability below before publishing the kernel's changes.
        if (lookupGetter.call(link.object, link.key)) throw REJECT;
        return link.object[link.key] === link.value;
      })
        || !plan.shapes.every(({value, length, prototype, isArray}) => value.length === length
          && (isArray ? Object.getPrototypeOf(value) === prototype && Object.isExtensible(value)
            : Object.getPrototypeOf(value) === prototype)))) plan = null;
      const originalArrays = plan?.originalArrays ?? arrayEntries.map(([name, rank], i) => {
        const value = linked(ctx, name);
        const kind = name === 'TFORCE' ? 'bool' : integerArrays.has(name) ? 'int' : 'float';
        descriptors.push({offset: arraysOffset + i * 8, value: pack(value, rank, kind)});
        return value;
      });
      const originalSources = plan?.originalSources ?? layout.sources.map((name, i) => {
        const value = linked(ctx, name);
        if (value != null && !plain(value)) throw REJECT;
        const offset = value ? allocate(48) : 0;
        if (value) {
          for (let j = 0; j < sourceKeys.length; j++) descriptors.push({offset: offset + j * 8,
            value: pack(linked(value, sourceKeys[j]), 1, j < 2 ? 'int' : 'float')});
        }
        return {value, offset};
      });
      if (plan) {
        cursor = plan.cursor;
        allocations = plan.allocations;
        descriptors = plan.descriptors;
        objects = plan.objects;
      } else {
        const written = new Set();
        const markWritten = (value, rank) => {
          if (value == null) return;
          if (rank === 1) written.add(ArrayBuffer.isView(value) ? value.buffer : value);
          else for (const row of value) markWritten(row, rank - 1);
        };
        for (let i = 0; i < arrayEntries.length; i++) {
          if (!readOnlyArrays.has(arrayEntries[i][0])) markWritten(originalArrays[i], arrayEntries[i][1]);
        }
        for (const {value} of originalSources) if (value) {
          for (const key of sourceKeys) markWritten(value[key], 1);
        }
        for (const allocation of allocations) {
          allocation.writable = written.has(ArrayBuffer.isView(allocation.original)
            ? allocation.original.buffer : allocation.original);
        }
        plan = {cursor, allocations, descriptors, objects, links, shapes, originalArrays, originalSources,
          blocks: descriptorBlocks(descriptors)};
        plans.set(ctx, plan);
      }
      // Local spline arrays and newly needed scratch arrays live after inputs.
      const arenaStart = cursor;
      const required = Math.min(MAX_BYTES, cursor + Math.max(1024 * 1024, cursor - heapBase.value));
      if (required > memory.buffer.byteLength) memory.grow(Math.ceil((required - memory.buffer.byteLength) / 65536));
      const doubles = new Float64Array(memory.buffer);
      const words = new Uint32Array(memory.buffer);
      for (const {offset, data, original, kind} of allocations) {
        if (Array.isArray(original)) for (let i = 0; i < data.length; i++) {
          const value = ownData(original, String(i));
          if (typeof value !== (kind === 'bool' ? 'boolean' : 'number')) throw REJECT;
          data[i] = Number(value);
        }
        // Int32 storage is finite by construction; keep the large Float64 scan
        // specialized instead of mixing typed-array element kinds in its loop.
        if (!(data instanceof Int32Array) && !finiteArray(data)) throw REJECT;
        nativeSet.call(doubles, data, offset / 8);
      }
      // Reuse descriptors only after a successful call proved that the kernel
      // left them intact. Different contexts, rejected calls and replaced rows
      // must republish their own metadata, even when offsets happen to match.
      if (residentPlan !== plan) {
        for (const {offset, values} of plan.blocks) wordSet.call(words, values, offset);
      }
      residentPlan = undefined;
      for (let i = 0; i < layout.scalars.length; i++) {
        const value = ownData(ctx, layout.scalars[i]);
        if (value != null && typeof value !== (booleanFields.has(layout.scalars[i]) ? 'boolean' : 'number')) throw REJECT;
        if (value != null && !Number.isFinite(Number(value))) throw REJECT;
        doubles[(scalarsOffset + i * 16) / 8] = value == null ? NaN : Number(value);
        words[(scalarsOffset + i * 16 + 8) / 4] = value == null ? 0 : 1;
        words[(scalarsOffset + i * 16 + 12) / 4] = 0;
      }
      for (let i = 0; i < originalSources.length; i++) {
        const {value, offset} = originalSources[i];
        words[sourcesOffset / 4 + i] = offset;
        if (offset) {
          if (!plain(value)) throw REJECT;
          const count = ownData(value, 'count') ?? 0;
          if (!Number.isFinite(count)) throw REJECT;
          doubles[(offset + 40) / 8] = count;
        }
      }
      logs = [];
      if (timings) timings.push(performance.now());
      const kernel = solve ? instance.exports.boundary_layer_step : assemble;
      kernel(scalarsOffset, arraysOffset, sourcesOffset, arenaStart, memory.buffer.byteLength);
      if (timings) timings.push(performance.now());

      // Validate and prepare every write before publishing any kernel result.
      commits = [];
      const sameLayout = sameDescriptors(words, plan.blocks)
        && originalSources.every(({offset}, i) => words[sourcesOffset / 4 + i] === offset);
      // Copy each backing buffer once. In particular VM and VMF are views of
      // the same matrix: validating/copying every row again is unnecessary.
      for (const {offset, data: input, original, writable} of allocations) {
        if (!writable) continue;
        const data = new Float64Array(memory.buffer, offset, input.length);
        if (!finiteArray(data)) throw REJECT;
        if (original instanceof Int32Array) {
          for (let i = 0; i < data.length; i++) if (!Number.isInteger(data[i]) || data[i] < -2147483648 || data[i] > 2147483647) throw REJECT;
          commits.push(() => integerSet.call(original, data));
        } else if (original instanceof Float64Array) {
          commits.push(() => nativeSet.call(original, data));
        } else {
          for (let i = 0; i < data.length; i++) if (typeof original[i] === 'boolean' && data[i] !== 0 && data[i] !== 1) throw REJECT;
          commits.push(() => { for (let i = 0; i < data.length; i++) original[i] = typeof original[i] === 'boolean' ? !!data[i] : data[i]; });
        }
      }
      const unpacked = new Map();
      const validRange = (pointer, length, stride = 8) => {
        if (!pointer || pointer % 8 || pointer < scalarsOffset || pointer + length * stride > memory.buffer.byteLength) throw REJECT;
      };
      function unpack(pointer, length, rank, original, kind = 'float') {
        if (!pointer) { if (length) throw REJECT; return original; }
        validRange(pointer, length);
        const packed = original && objects.get(original);
        if (rank === 1 && packed?.descriptor[0] === pointer && packed.descriptor[1] === length) return original;
        const key = `${pointer}:${length}:${rank}`;
        if (unpacked.has(key)) {
          const seen = unpacked.get(key);
          if (original && seen.original && original !== seen.original) {
            // Distinct typed views may share bytes; retain both view identities.
          } else return seen.result;
        }
        if (rank === 1) {
          const data = new Float64Array(memory.buffer, pointer, length);
          if (!finiteArray(data)) throw REJECT;
          const integer = original instanceof Int32Array || kind === 'int';
          const booleans = kind === 'bool';
          if (integer && !data.every(v => Number.isInteger(v) && v >= -2147483648 && v <= 2147483647)) throw REJECT;
          if (booleans && !data.every(v => v === 0 || v === 1)) throw REJECT;
          let result = original;
          if (!result || result.length !== length) result = booleans ? new Array(length) : integer ? new Int32Array(length) : new Float64Array(length);
          commits.push(() => {
            if (Array.isArray(result)) for (let i = 0; i < length; i++) result[i] = booleans ? !!data[i] : data[i];
            else if (result instanceof Int32Array) integerSet.call(result, data);
            else nativeSet.call(result, data);
          });
          unpacked.set(key, {original, result});
          return result;
        }
        let result = original;
        if (!result || result.length !== length) result = new Array(length);
        unpacked.set(key, {original, result});
        for (let i = 0; i < length; i++) {
          const child = unpack(words[pointer / 4 + i * 2], words[pointer / 4 + i * 2 + 1], rank - 1, result[i], kind);
          if (child !== result[i] || !Object.hasOwn(result, i)) {
            ownData(result, String(i));
            commits.push(() => { result[i] = child ?? null; });
          }
        }
        return result;
      }
      for (let i = 0; i < layout.scalars.length; i++) {
        const offset = scalarsOffset + i * 16;
        if (!words[(offset + 12) / 4]) continue;
        const name = layout.scalars[i], value = doubles[offset / 8];
        if (!Number.isFinite(value) || (booleanFields.has(name) && value !== 0 && value !== 1)) throw REJECT;
        commits.push(() => { ctx[name] = booleanFields.has(name) ? !!value : value; });
      }
      for (let i = 0; !sameLayout && i < arrayEntries.length; i++) {
        const [name, rank] = arrayEntries[i], offset = arraysOffset / 4 + i * 2;
        const value = unpack(words[offset], words[offset + 1], rank, originalArrays[i], name === 'TFORCE' ? 'bool' : integerArrays.has(name) ? 'int' : 'float');
        if (value !== originalArrays[i]) commits.push(() => { ctx[name] = value; });
      }
      for (let i = 0; i < originalSources.length; i++) {
        const pointer = words[sourcesOffset / 4 + i];
        if (!pointer) continue;
        validRange(pointer, 6);
        const original = originalSources[i].value;
        const value = original ?? {};
        for (let j = 0; !sameLayout && j < sourceKeys.length; j++) {
          const key = sourceKeys[j], offset = pointer / 4 + j * 2;
          const data = unpack(words[offset], words[offset + 1], 1, value[key], key === 'panel' || key === 'system' ? 'int' : 'float');
          if (data !== value[key]) commits.push(() => { value[key] = data; });
        }
        if (i === 1) {
          const count = doubles[(pointer + 40) / 8];
          if (!Number.isInteger(count) || count < 0 || count > ctx.NSYS) throw REJECT;
          commits.push(() => { value.count = count; });
        }
        if (!original) commits.push(() => { ctx[layout.sources[i]] = value; });
      }
      if (sameLayout) residentPlan = plan;
    } catch (error) {
      residentPlan = undefined;
      stats.fallback++;
      logs = [];
      running = false;
      onFallback?.(error);
      return false;
    }
    if (timings) timings.push(performance.now());
    for (const commit of commits) commit();
    if (timings) { timings.push(performance.now()); onTiming(timings.slice(1).map((t,i)=>t-timings[i])); }
    stats.accelerated++;
    running = false;
    for (const [id, values] of logs) console.log(logFormats[id](values));
    logs = [];
    return true;
  }
  tryAssemble.stats = stats;
  return tryAssemble;
}

let accelerator;
export function tryWasmBoundaryLayer(ctx, solve = false) {
  if (accelerator === undefined) accelerator = createWasmBoundaryLayer();
  return accelerator ? accelerator(ctx, solve) : false;
}

export function getBoundaryLayerWasmStats() {
  return accelerator ? {...accelerator.stats} : {calls: 0, accelerated: 0, fallback: 0};
}
