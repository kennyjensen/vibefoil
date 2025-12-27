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
