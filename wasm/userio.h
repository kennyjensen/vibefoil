// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_USERIO_H
#define WASM_USERIO_H

#include <string>
#include <utility>

std::pair<std::string, int> strip_string(const std::string &input);

#endif  // WASM_USERIO_H
