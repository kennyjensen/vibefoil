# Vibefoil

Vibefoil is a browser-based port of XFOIL. Learn more about XFOIL at https://web.mit.edu/drela/Public/web/xfoil/. The live site is https://vibefoil.com.

## Features

- Numerically faithful port to JavaScript
- Runs locally in the browser
- Interactive UI
- Supports normal and inverse design

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
required. Plain local static servers need a custom-404 fallback to exercise these
links; otherwise use `/?_airfoil_path=/naca2412` for local testing.
