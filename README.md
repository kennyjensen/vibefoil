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
