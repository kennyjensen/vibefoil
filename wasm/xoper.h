#ifndef WASM_XOPER_H
#define WASM_XOPER_H

struct XFoilState;
struct XBlState;

void mhinge(XFoilState &ctx);
void viscal(XFoilState &ctx, XBlState &bl, int niter1);

#endif  // WASM_XOPER_H
