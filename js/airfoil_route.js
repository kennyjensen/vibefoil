const ROUTE_PARAM = '_airfoil_path';

// Only a single airfoil name is a route; missing assets and nested paths stay 404s.
export function airfoilNameFromPath(pathname, basePath = '/') {
  if (!pathname.startsWith(basePath)) return null;
  let name;
  try {
    name = decodeURIComponent(pathname.slice(basePath.length).replace(/\/$/, ''));
  } catch {
    return null;
  }
  if (!/^[a-z0-9][a-z0-9._-]*$/i.test(name)) return null;
  if (/^(index\.html|404\.html)$/i.test(name)) return null;
  if (/\.(html?|js|css|json|ico|png|svg|wasm)$/i.test(name)) return null;
  return name.toLowerCase();
}

// GitHub Pages serves 404.html for clean URLs. Visit the real entry page once,
// then restore the requested path after its assets have a stable base URL.
export function airfoilRedirectURL(url, basePath = '/') {
  if (!airfoilNameFromPath(url.pathname, basePath)) return null;
  const target = new URL(url);
  target.searchParams.set(ROUTE_PARAM, target.pathname);
  target.pathname = basePath;
  return target;
}

export function restoredAirfoilURL(url, basePath = '/') {
  const path = url.searchParams.get(ROUTE_PARAM);
  if (!path || !airfoilNameFromPath(path, basePath)) return null;
  const target = new URL(url);
  target.pathname = path;
  target.searchParams.delete(ROUTE_PARAM);
  return target;
}

// Match the generator families supported by the UI. Other designations fall
// through to the database, including modified/subscripted NACA sections.
export function nacaSettingsForName(name) {
  // An explicit .dat suffix requests the database's coordinates.
  const match = /^naca(?:-|_)?(.+)$/i.exec(name);
  if (!match) return null;
  const code = match[1];
  let parts = /^(\d)(\d)(\d{2})$/.exec(code);
  if (parts && Number(parts[3]) > 0) {
    return { mode: '4', m: Number(parts[1]), p: Number(parts[2]), t: Number(parts[3]) };
  }
  parts = /^(2[1-5]0)(\d{2})$/.exec(code);
  if (parts && Number(parts[2]) > 0) {
    return { mode: '5', series5: parts[1], t5: Number(parts[2]) };
  }
  parts = /^(6[3-7]|6[3-5]a)-?(\d)(\d{2})$/i.exec(code);
  if (parts && Number(parts[3]) > 0) {
    return { mode: '6', profile6: parts[1].toUpperCase(), cl6: Number(parts[2]) / 10, t6: Number(parts[3]) };
  }
  return null;
}
