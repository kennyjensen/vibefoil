// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#include "userio.h"

std::pair<std::string, int> strip_string(const std::string &input) {
    const int nlen = static_cast<int>(input.size());

    int k2 = 0;
    for (int idx = nlen; idx >= 1; --idx) {
        if (input[idx - 1] != ' ') {
            k2 = idx;
            break;
        }
    }

    int k1 = 1;
    for (int idx = 1; idx <= k2; ++idx) {
        if (input[idx - 1] != ' ') {
            k1 = idx;
            break;
        }
    }

    const int ns = k2 - k1 + 1;
    if (ns == 0) {
        return {input, 0};
    }

    std::string stripped = input.substr(static_cast<size_t>(k1 - 1), static_cast<size_t>(k2 - k1 + 1));
    stripped.append(static_cast<size_t>(nlen - static_cast<int>(stripped.size())), ' ');
    return {stripped, ns};
}
