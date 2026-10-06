import assert from 'node:assert/strict';
import {
  airfoilNameFromPath,
  airfoilRedirectURL,
  restoredAirfoilURL,
  nacaSettingsForName,
} from '../js/airfoil_route.js';

for (const [path, name] of [
  ['/naca2412', 'naca2412'],
  ['/NACA0012/', 'naca0012'],
  ['/sd7003.dat', 'sd7003.dat'],
  ['/fx63-137', 'fx63-137'],
  ['/naca001034a08cli0.2', 'naca001034a08cli0.2'],
  ['/%6Eaca2412', 'naca2412'],
]) {
  assert.equal(airfoilNameFromPath(path), name);
}
for (const path of [
  '/', '/index.html', '/404.html', '/favicon.ico', '/missing.js',
  '/js/missing.js', '/naca2412/extra', '//example.com', '/%2Fsd7003',
  '/%E0%A4%A', '/..', '/<script>',
]) {
  assert.equal(airfoilNameFromPath(path), null, path);
}
assert.equal(airfoilNameFromPath('/vibefoil/s1223/', '/vibefoil/'), 's1223');
assert.equal(airfoilNameFromPath('/another/s1223', '/vibefoil/'), null);

// The Pages detour retains the path, query values and fragment, including
// encoded delimiters, with no redirect loop on the actual entry page.
for (const [address, basePath] of [
  ['https://vibefoil.com/naca0012', '/'],
  ['https://vibefoil.com/SD7003.dat/?label=a%26b%3Dc&label=two#plot', '/'],
  ['https://kennyjensen.github.io/vibefoil/naca63-212/?ref=shared#plot', '/vibefoil/'],
]) {
  const original = new URL(address);
  const redirect = airfoilRedirectURL(original, basePath);
  assert.equal(redirect.pathname, basePath);
  assert.equal(redirect.origin, original.origin);
  assert.equal(airfoilRedirectURL(redirect, basePath), null);
  const restored = restoredAirfoilURL(redirect, basePath);
  assert.equal(restored.pathname, original.pathname);
  assert.deepEqual([...restored.searchParams], [...original.searchParams]);
  assert.equal(restored.hash, original.hash);
  assert.equal(restoredAirfoilURL(restored, basePath), null);
}
for (const path of ['//example.com', '/js/app.js', '/naca2412/extra', '/%2Fexample.com']) {
  const url = new URL('https://vibefoil.com/');
  url.searchParams.set('_airfoil_path', path);
  assert.equal(restoredAirfoilURL(url), null);
}

assert.deepEqual(nacaSettingsForName('naca2412'), { mode: '4', m: 2, p: 4, t: 12 });
assert.deepEqual(nacaSettingsForName('NACA0012'), { mode: '4', m: 0, p: 0, t: 12 });
assert.deepEqual(nacaSettingsForName('naca0040'), { mode: '4', m: 0, p: 0, t: 40 });
for (const series5 of ['210', '220', '230', '240', '250']) {
  assert.deepEqual(nacaSettingsForName(`naca${series5}15`), { mode: '5', series5, t5: 15 });
}
for (const profile6 of ['63', '64', '65', '66', '67', '63A', '64A', '65A']) {
  for (const separator of ['', '-']) {
    assert.deepEqual(nacaSettingsForName(`naca${profile6}${separator}212`), {
      mode: '6', profile6, cl6: 0.2, t6: 12,
    });
  }
}
assert.deepEqual(nacaSettingsForName('naca64a010'), { mode: '6', profile6: '64A', cl6: 0, t6: 10 });
// Database-only variants must not be approximated as a different NACA section.
for (const name of ['sd7003', 'naca2412.dat', 'naca633218', 'naca001234', 'naca23112', 'naca0000', 'naca68-212']) {
  assert.equal(nacaSettingsForName(name), null, name);
}

console.log('airfoil_route.test.js: OK');
