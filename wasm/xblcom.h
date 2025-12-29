// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XBLCOM_H
#define WASM_XBLCOM_H

#include <string>
#include <vector>

struct XBlState;

extern const int NCOM;
extern const std::vector<std::string> COM1_NAMES;
extern const std::vector<std::string> COM2_NAMES;

void sync_com_to_vars(XBlState &bl, int which);
void sync_vars_to_com(XBlState &bl, int which);

#endif  // WASM_XBLCOM_H
