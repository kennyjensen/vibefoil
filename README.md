# Vibefoil

Vibefoil is a browser-based port of XFOIL. Learn more about XFOIL at https://web.mit.edu/drela/Public/web/xfoil/. The live site is https://vibefoil.com.

## Features

- Numerically faithful port to JavaScript
- Runs locally in the browser
- Interactive UI
- Supports normal and inverse design

## Run locally

With Python 3 installed, run `npm start` (or `python3 serve.py`) from this
directory, then open http://localhost:8080. Airfoil links such as
http://localhost:8080/naca0012 work through the same custom `404.html` page used
on GitHub Pages. To choose another port, run `python3 serve.py --port 8081`.

Use this command in place of `python3 -m http.server`, which serves its own
generic 404 page and cannot open airfoil paths directly.

## Numerical verification

The solver tracks the checked-in XFOIL 6.996 source. See the
[numerical audit and test instructions](python/tests/README.md#xfoil-numerical-audit)
for corrections, Fortran/native/WASM comparisons, precision settings, tolerances,
and the limits of feature equivalence.

## Performance checks

Run `npm run benchmark` to time the browser worker's calculation path in Node.
It measures inviscid and viscous NACA 2412 cases, a 10-degree flap, and a
three-angle sweep with solution reuse, including streamline generation. It uses
240 panels, the existing 60-by-32 streamline grid, and unchanged solver settings.
Viscous cases use Mach 0.1, Reynolds number 3 million, Ncrit 9, and a 30-iteration
limit; the sweep runs at 2, 4, and 6 degrees.
The benchmark excludes browser rendering and worker message transfer.

For comparisons, run `npm run benchmark:compare -- --baseline /path/to/checkout`.
It alternates two isolated workers, reverses their execution order each round,
and checks complete output equality after every pair, including warm-up rounds.
Where Node supports it, the report also includes CPU time used by the solver
thread, separately from elapsed time, to help identify scheduling noise.
On machines with different CPU core types, pin the entire command to one
available core with `taskset -c <core>` to reduce timing noise.

On 2026-10-11, Node 24.16.0 on arm64, pinned to CPU 3, gave these median paired
reductions over nine measured rounds after four warm-up rounds. The baseline
already included source/streamline caching and contiguous matrix elimination;
this comparison measures the additional matrix-assembly changes.

| Workload | Median paired reduction in elapsed time |
| --- | ---: |
| Inviscid | approximately 0% |
| Viscous | 0.8% |
| Flapped airfoil | 6.2% |
| Three-angle viscous sweep | 5.3% |

The inviscid path does not use boundary-layer assembly. Small differences are
within timing variation; timings depend on hardware, runtime, and system load.

The changes gather boundary-layer source data once per evaluation and prepare
panel constants once per streamline grid, calculating only the streamfunction
needed by the display. The viscous matrix solver temporarily transposes its
coefficient planes in place so its inner loops read contiguous memory, restoring
the public layout even on failure. This JavaScript path needs no extra matrix
allocation; the optional WebAssembly path described below uses scratch memory.
Assembly also captures source indices/signs once per call and reads each row's
constant coefficients outside the loop over source stations.
It now assembles all three equation rows in one pass, sharing source reads and
the repeated leading-edge term. Tightly sized tensors use the existing flat
buffer directly; nested or oversized storage retains its row views. Compared
with the preceding optimizations, seven paired rounds after three warm-up rounds
on Node 24.16.0/arm64, pinned to CPU 3, measured 1.4% less time for the viscous
case, 1.1% for the flapped case, and 2.1% for the sweep; inviscid timing was
unchanged. These are small incremental gains and vary with system load.
The arithmetic within each coefficient, solver tolerances, and precision are
preserved.
All returned worker data matched the baseline exactly in these same-runtime
comparisons. A separate Chromium worker comparison also matched 237,185 numeric
values exactly across NACA 2413, NACA 23012, and database SD7003 cases, including
solution reuse and inviscid mode.
All JavaScript checks and eight assembly/viscous numerical regression tests passed,
including the complete operating-point comparisons with double-precision XFOIL.

A later elimination change shares each source-row read across independent
destination rows when the second and third updates are active. It preserves
each coefficient's multiply/subtract expression and the strict acceleration
cutoffs. Compared with the already-optimized worker above, seven paired rounds
after three warm-up rounds on the same runtime/core measured 4.2% less time for
the viscous case and 4.2% for the sweep; inviscid and flap timings showed no
meaningful improvement. Complete worker outputs matched exactly on every pair.

The contiguous mass-row updates now use small helpers that process four entries
per loop, with scalar tails for the remaining entries. Each coefficient keeps
the same multiply/subtract order, and acceleration cutoffs, matrix layouts, and
precision are unchanged. An isolated comparison of six captured solver states
measured 8.2% less matrix-solver thread CPU time on Node 24.16.0/arm64 (21 paired
samples after ten warm-up samples, CPU 3). A Chromium 145 check using a fresh
worker for every solve measured 19.6% less total time in matrix elimination
during each worker's first airfoil calculation (nine paired samples, NACA 2412
angles cycling from 0 to 4 degrees, CPU 3). Complete Node solves showed no
reliable overall speed change; the matrix timings measure only that stage.
All 151 direct matrix comparisons, the complete JavaScript suite, and seven
matrix/viscous numerical reference tests passed, including double-precision
XFOIL operating points.
Browser comparisons matched over four million numeric values, including saved
state, changed Ncrit, and inverse design. Rendered canvases and Cp/BL exports
also matched through mode changes and resizing.

Larger flat boundary-layer systems (32–1024 stations) now use a small standalone
WebAssembly kernel when the browser supports it. The kernel retains Float64
arithmetic, operation order, pivot checks, and acceleration cutoffs. It compiles
lazily on the first viscous solve; inviscid work does not instantiate it. Separate
scratch memory keeps the original inputs intact until a finite solution succeeds.
Unsupported runtimes, unusual/aliased storage, singular blocks, or nonfinite
results run the original JavaScript routine from those intact inputs, preserving
its failure behavior. Scratch capacity is bounded to 64 MiB per worker and reused
between calls (about 2.4 MiB for the default 270-station system).

Against the preceding JavaScript implementation, Chromium 145 measured 46.0%
less time in matrix elimination and 6.4% less time in complete warm viscous worker
calls (11 paired samples after four warm-up pairs, NACA 2412 at 0–4 degrees,
solution reuse disabled, CPU 3). A separate first-solve comparison created, ran,
and terminated each worker sequentially: nine pairs measured 55.8% less matrix
time and 7.9% less complete solver time, including kernel compilation but excluding
worker/module startup and rendering. Inviscid timings were essentially unchanged.
Node 24.16.0/arm64 measured 6.7% less solver-thread CPU time for the fresh viscous
case, 5.7% for the flapped case, and 9.2% for the reused sweep (nine paired samples,
four warm-up pairs, CPU 3). Complete worker outputs matched exactly in every pair;
these are incremental, hardware-dependent measurements, not cumulative gains.

The full JavaScript suite and nine matrix/viscous numerical reference tests passed,
including double-precision XFOIL operating points. Browser checks matched over
four million numeric values across reused/restored state and inverse design;
canvas images and Cp/BL exports also matched through mode changes and resizing.

The kernel source is `wasm/blsolv_kernel.cpp`. Rebuild the generated embedded
module with `python3 wasm/build_blsolv.py`, using Clang and LLD with wasm32 support
(tested with LLVM 19). Use `--compiler /path/to/clang++ --linker /path/to/wasm-ld`
if they are outside PATH. The script disables fast math and fused multiply-add;
SIMD operates on independent coefficients. Static hosting needs no compilation,
additional binary request, or Emscripten runtime. Tests check the source fingerprint
and numerical results so edits cannot silently leave an outdated generated asset.

The WebAssembly row kernels now give the compiler separate, nonoverlapping row
spans instead of manually unrolling indexed updates through one matrix pointer.
The source and destination rows are distinct by construction; this lets the
compiler vectorize independent coefficients while preserving each coefficient's
multiply/subtract order. The strict arithmetic flags and JavaScript fallback are
unchanged. The compiled binary shrank from 13,830 to 11,715 bytes (15.3%).
Compared with the first WebAssembly version, six captured 270-station systems
took 22.5% less elapsed time in an isolated Node benchmark (21 paired samples,
ten warm-up pairs, CPU 3). Complete Node worker comparisons measured 4.9% less
solver-thread CPU time for fresh viscous and flapped cases and 8.0% for the reused
sweep (nine pairs, four warm-up pairs). Chromium 145 measured 26.1% less matrix
time in 11 pairs after four warm-up pairs, with exact complete outputs. Browser
total-solve and first-solve timings varied too much to claim another overall gain.
The existing 231 matrix comparisons, 67 storage/failure checks, full JavaScript
suite, and nine numerical reference tests passed with this kernel. Browser
comparisons again matched over four million numeric values, rendered canvases,
and Cp/BL exports exactly.

Boundary-layer assembly and transition calculations avoid three redundant
73-value state copies. Ordinary assembly uses the snapshot already published by
`BLVAR`; midpoint evaluation leaves those station-2 values intact. Transition
checks and transition differentials let `BLKIN` publish its completed snapshot
without first copying intermediate values that it would overwrite. Saved-state
restoration, numerical expressions, and solver tolerances are unchanged.
Direct comparisons against the preceding implementation matched the entire
routine state in 1,200 cases, including repeated calls, transition boundaries,
alternate array storage, unused entries, and nonfinite inputs.
Nine paired Node 24.16.0/arm64 worker comparisons after four warm-up rounds,
pinned to CPU 3, measured 4.2% less solver-thread CPU time for the fresh viscous
case, 3.7% for the flapped case, and 5.4% for the three-angle sweep. Inviscid
timing was essentially unchanged. All complete outputs matched exactly.
Browser timings were too variable to claim a reliable overall speedup.
All 132 snapshot/transition checks, the full JavaScript suite, and nine
boundary-layer/viscous numerical reference tests passed, including the
double-precision XFOIL operating points. Browser checks matched over four
million numeric values across fresh and reused solves, restored state, changed
Ncrit, and inverse design. Canvas images and Cp/BL exports also matched through
mode changes and resizing.

The remaining boundary-layer snapshots use direct assignments for the fixed
73-field COMMON layout, avoiding repeated property-name lookups. The copy order,
array conversions, and numerical calculations are unchanged. A regression test
reads the original XFOIL `XBL.INC` declarations and checks every offset, including
alternate array types, padding, partial storage, overlapping snapshots, special
numeric values, and observable access order (860 exact comparisons).
Five paired Node 24.16.0/arm64 comparisons after three warm-up rounds measured
2.6–3.2% less solver-thread CPU time for viscous, flapped, and sweep cases.
Chromium 145 measured 18.0% less boundary-layer assembly time and 8.8% less total
viscous solve time in 11 paired rounds after four warm-up rounds, using NACA 2412
at 0–4 degrees with solution reuse disabled. Both runs were pinned to CPU 3;
complete results matched exactly. Browser timings contained
outliers, so gains will vary. Nine additional pairs with a new worker for each
calculation showed no reliable change in overall first-calculation time.
All 1,200 entire-state comparisons against the previous implementation, the
JavaScript suite, and nine boundary-layer/viscous numerical reference tests
passed. Browser checks matched over four million numeric values, canvas images,
and Cp/BL exports, including reused/restored state and inverse design.

Built-in predictor/kinematics pairs now publish one completed 73-value snapshot
instead of copying both the intermediate and completed state. This applies to
both marching routines and system assembly. Standalone `BLPRV` still publishes
its snapshot; custom predictor or kinematics callbacks keep the original call
sequence, receivers, and snapshot visibility. All numerical expressions remain
shared with the standalone routines. Another 909 regression comparisons check
entire station state, repeated calls, array conversions/padding, special values,
method calls, callback replacement, and the halved snapshot-store count.
Chromium 145 measured 9.9% less time in `MRCHDU` in 11 paired rounds after four
warm-up rounds (NACA 2412 at 0–4 degrees, reuse disabled, CPU 3). Complete outputs
matched exactly. Whole-solve Node and browser gains were too small or variable
to claim a reliable improvement. First-calculation browser timings also varied
substantially, including a check that created and ran each worker separately
to avoid simultaneous startup interfering with the measurements.
The full JavaScript suite and nine numerical reference tests passed, including
double-precision XFOIL operating points. Direct before/after station comparisons
and over four million browser numeric values matched exactly. Canvas images and
Cp/BL exports matched through mode changes, reuse, inverse design, and resizing.

The small Gaussian solves used during boundary-layer marching now work directly
on the existing 1-based Float64 rows. This removes five temporary typed arrays
and an outer row array per call while preserving pivot selection and every
changed matrix/RHS value. Nine paired worker rounds after four warm-up rounds
showed small reductions (1.4% viscous, 0.8% flapped, 0.1% sweep); these are close
to timing noise. All complete outputs matched exactly.
The final implementation passed eight marching/viscous regression tests,
including double-precision XFOIL operating points. Chromium again matched all
237,185 compared numeric values, rendered canvases, and Cp/BL downloads exactly.

Streamline evaluation also shares each endpoint's computed logarithm with the
next adjacent panel. Skipped panels break that reuse, and each evaluation point
starts independently. This nearly halves endpoint logarithm calls while keeping
the original expressions, summation order, and grid resolution.
Compared with the preceding optimizations, Chromium 145 measured 8–9% less time
building streamline grids (about 31 ms to 28 ms). Complete inviscid worker calls
took 5.5% less time; viscous timing was essentially unchanged. These are median
paired reductions from 11 samples after four warm-up samples, NACA 2412 angles
cycling from 0 to 4 degrees with solution reuse disabled, pinned to CPU 3.
Complete outputs matched exactly. Node's isolated evaluator measured a 14.7%
reduction, but whole-worker Node timings varied too much for a reliable claim.
The full JavaScript suite passed, and browser canvas images and Cp/BL exports
matched exactly through mode changes and resizing.

The full panel influence routine also reuses shared-endpoint logarithms,
including geometry/source derivatives and reflected geometry. Self-influence
guards remain unchanged, skipped panels break reuse, and the reflected loop
starts with its own distances. A 31-panel example uses 32 endpoint logarithms
instead of 62; source-midpoint logarithms are unchanged. Independent comparisons
against the preceding implementation matched all returned numeric values and
the entire context in 8,943 cases, including derivative arrays and unused entries.
An isolated Node 24.16.0/arm64 benchmark measured 4–7% less elapsed time for panel
evaluation (21 paired samples after ten warm-up samples, three passes over 240
surface nodes per sample, CPU 3). Complete solve differences were around 1% or
less, and browser timings were too variable to claim a reliable overall gain.
The full JavaScript suite and 19 panel/viscous numerical reference tests passed,
including double-precision XFOIL operating points. Browser comparisons matched
over four million numeric values across fresh/reused solves, restoration, changed
Ncrit, and inverse design. Canvas images and Cp/BL exports also matched through
mode changes and resizing.

Repeated panel setup can reuse one LU factorization after checking every entry
of the assembled matrix with `Object.is`. Any changed entry runs the original
factorization. Reuse preserves row-object swaps, pivots, and unused padding;
other storage layouts and nonfinite matrices use the original routine. Matrix
assembly, right-hand sides, and boundary-layer initialization still run normally.
The bounded cache adds about 0.93 MB for the default 240-panel geometry.
Chromium 145 measured median repeated inviscid solve time of 70.4 ms before and
45.9 ms after (32.5% median paired reduction, 11 pairs after four warm-up pairs,
NACA 2412 angles cycling from 0 to 4 degrees, pinned to CPU 3). This times the
complete worker handler, including sending its result, with solution reuse
disabled; it excludes main-thread rendering. Viscous timings were too variable
to claim a reliable gain. All worker outputs matched exactly, including fresh
solves, restored snapshots, and inverse design. Canvas images and Cp/BL exports
also matched, and 22 numerical reference tests passed.
An isolated Node 24.16.0/arm64 comparison (11 pairs after four warm-up rounds,
CPU 3) measured 17.7% less time for inviscid solves and 2.3% for viscous solves.
The geometry-changing flap case was 1.2% slower and the reused sweep was nearly
unchanged, so this cache's main benefit is repeated fresh panel setup.

Wake-coupling back substitution processes four independent right-hand sides
together, sharing matrix reads while preserving each column's pivot sequence,
leading-zero handling, and arithmetic order. Four Float64 scratch vectors serve
all columns, replacing 270 temporary vectors for the default 240-panel geometry
and 30 wake points.
Chromium 145 measured 21.5% less time in this coupling calculation (26.2 ms to
21.2 ms), and a 4.7% median paired reduction in complete fresh viscous solves
(11 pairs after four warm-up pairs, NACA 2412 angles cycling from 0 to 4 degrees,
CPU 3). Reused sweeps showed no reliable overall timing improvement. Nine Node
24.16.0/arm64 pairs after four warm-up rounds on CPU 3 measured 2.4% less solver
thread CPU time for fresh viscous solves and 1.3% for the flapped case; inviscid
and reused-sweep differences were small. Complete worker outputs matched exactly,
including all 5,370,864 numeric values compared in the browser sweep.
The full JavaScript suite and 22 numerical reference tests passed, including
double-precision XFOIL operating points. Rendered canvases and Cp/BL exports
also matched through mode changes and resizing.

Boundary-layer updates store scratch velocities in contiguous vectors for each
side. For the default 240-panel geometry and 30 wake points, this constructs six
temporary Float64 arrays per update instead of 546. Numeric scratch storage
falls from 16,912 to 12,560 bytes; unused side-zero entries are no longer stored.
Every call still starts with fresh, zeroed arrays and uses the same arithmetic.
A separate before/after check matched all state fields through 289 updates,
including changed active regions, repeated calls, holes in the surface mapping,
signed zeros, and nonfinite inputs; allocation counting confirmed the reduction.
Complete worker outputs and browser iteration counts also matched exactly.
Whole-solver timing differences were small; no reliable end-to-end speedup was
established for this change.
The full JavaScript suite and six viscous/reference tests passed, including
double-precision XFOIL operating points. Browser plots and exports matched
exactly through mode changes and resizing.

Startup avoids evaluating numerical code on the main thread: the solver worker
loads those modules. Compared with the original startup, Chromium evaluates 10
fewer main-thread scripts/modules (about 307 KB less source), including the small
cache-tracking and contour modules added below. Initial and changed-angle
results, rendered canvases, and downloaded Cp/BL data matched the baseline exactly.

Streamline redraws retain one grid's contour coordinates and check the grid,
domain, and contour levels for changes, including in-place edits. Rendering
still applies the current rotation, zoom, pan, and airfoil clipping on every
draw. Grid resolution, the 14 contour levels, interpolation, and segment order
are unchanged. Chromium 145 measured median contour-rendering time of 0.5 ms
before and 0.2 ms after, with the complete redraw function falling from 2.8 ms
to 2.1 ms (16 paired zoom samples after four warm-up samples, pinned to CPU 3).
These timings cover synchronous drawing work, not the entire displayed frame.
All four canvases matched pixel for pixel in 26 screen states, including zoom,
pan, node dragging, edited geometry, and resizing; Cp/BL exports also matched.

Reusable snapshots store each matrix in one Float64 buffer with separate row
views. The worker transfers owned output buffers and retains live solver arrays.
For a NACA 2412 result with saved state, this reduced the message's buffer count
from 5,651 to 56 while preserving all 1,753,464 bytes of numeric storage. Chromium
145 measured median worker-side `postMessage` time of 3.4 ms before and 1.2 ms
after (11 paired samples after four warm-up rounds, angles cycling from 0 to 4
degrees). This measures sending the result, not total solve or render time.
The browser comparison matched 3,766,551 numeric values exactly across repeated
solves, cold snapshot restoration, changed Ncrit, downloads, and inverse design.
Rendered canvases and exports also matched, including after resizing.

Run `npm run benchmark:messages -- --baseline /path/to/checkout` for the
equivalent comparison using real Node worker messages. It reports serialization
time, buffer counts, byte counts, and checks complete output equality on each
paired request. Node and browser timings can differ.

Requests also avoid copying a saved solution back to a worker that explicitly
acknowledged matching live panel and boundary-layer caches. The page retains its
full backup for restoration; worker replacement, errors, and inverse design
invalidate the acknowledgment. Inviscid solves and advanced solves with reuse
disabled omit snapshots because those paths do not use them.
For repeated NACA 2412 solves, this removes 52 buffers totaling 1,734,240 bytes
from each request. Chromium 145 measured median page-side `postMessage` time of
5.0 ms before and 1.8 ms after (11 paired samples after four warm-up rounds,
angles cycling from 0 to 4 degrees, pinned to CPU 3). This measures sending a
request, not the complete calculation. All 3,766,551 compared numeric values
matched exactly, including restored state, changed Ncrit, downloads, and inverse
design. The only added response field is opt-in cache-readiness metadata.
Use `npm run benchmark:messages -- --direction request --baseline /path/to/checkout`
to measure this request path with real Node worker messages.

Use `node benchmarks/solver.mjs --root /path/to/checkout --runs 5` to benchmark
another revision. Add `--save /tmp/solver-results.bin` to save the full outputs
with Node's V8 serializer for exact comparisons; pass
`--compare /tmp/solver-results.bin` on a later run to check those outputs.
`npm test` also compares the streamline evaluator with the original panel routine
at 10,248 points, including sharp and blunt edges, rotated foils, duplicate points,
reflected geometry, extreme coordinates, and repeated evaluations after skipped
panels. A work-count check verifies the reduction in logarithm calls.
Contour tests compare exact segment coordinates/order in 70 cases, including
all 16 marching-square configurations, grid/domain/level changes, nonfinite
values, signed zeros, and cache reuse across camera changes.
It also checks the flat matrix solver against the nested-array reference in 231
cases, including different system sizes, all eight acceleration-threshold masks,
every four-entry loop remainder, repeated calls, and singular blocks; all
coefficient arrays must match exactly. The Fortran BLSOLV comparison also runs
the flat storage path and checks every changed coefficient against the nested
reference. The same cases exercise the WebAssembly kernel explicitly, including
208 successful accelerated solves. Another 67 storage/failure/byte-parity cases
cover shifted views, memory growth, signed zeros, subnormal values, NaN payloads,
metadata aliasing, and complete worker results with WebAssembly disabled.
Another 132 checks cover boundary-layer snapshot consistency, stale input
snapshots, midpoint/similarity updates, transition restoration, unused entries,
and the reduced number of snapshot stores.
The fixed COMMON copies have 860 comparisons against the original Fortran
field declarations, including every offset, conversions, padding, overlapping
storage, partial arrays, and property access order.
Another 367 exact comparisons check the direct Gaussian solve against the prior
copying wrapper, including pivot ties, singular systems, nonfinite values, and
untouched padding.
LU caching adds 361 exact matrix/back-substitution comparisons and 72 entire
panel-state comparisons, including changed geometry, signed zeros, reflected
geometry, row identities, padding, and fallback storage types.
Another 796 comparisons check batched back substitution against independent
single-column solves, covering pivots, different leading-zero lengths, singular
systems, nonfinite values, repeated calls, leftover columns, and untouched padding.
Another 254 exact comparisons check combined mass-matrix assembly against the
previous separate-row expressions, covering flat and nested storage, repeated
assembly, nonfinite values, and oversized buffers with inactive entries.
Worker-data tests check Float64 values and buffer ownership. Worker lifecycle
tests compare ordinary copying with actual buffer transfers through ten solve,
restore, download, and inverse-design operations, so detaching live solver state
causes a regression failure.
Request-reuse tests add 20 exact worker comparisons, including retained backups
after failed/unsaved solves, inviscid mode, inverse design, cold restoration, and
an incompatible snapshot that returns a result without populating the live
cache. Protocol checks also cover stale responses and queued inverse requests.

## Airfoil links

Open an airfoil directly by appending its name to the site URL:

- `https://vibefoil.com/naca2412` — NACA 4-digit section.
- `https://vibefoil.com/naca23012` — supported NACA 5-digit section.
- `https://vibefoil.com/naca63-212` or `/naca63212` — NACA 6-series section.
- `https://vibefoil.com/naca64a210` — NACA 6A-series section.
- `https://vibefoil.com/sd7003` — airfoil from the UIUC database.
- `https://vibefoil.com/naca2412.dat` — explicitly use database coordinates.

Names are case-insensitive, and a trailing slash is optional. NACA designations
use the existing 4-digit, 5-digit (210–250 camber lines), and 6-series generators;
other names are fetched from the database. Failed lookups show an error in the
Database panel so the name can be corrected and fetched again.

GitHub Pages must publish the root `404.html` alongside `index.html`. It redirects
airfoil paths through the entry page, which restores the clean URL without an
extra browser-history entry. No server rewrite configuration or build step is
required. The local server above provides this fallback too. With a plain static
server, use `/?_airfoil_path=/naca2412` for local testing instead.
