// Generated from current JavaScript by generate_boundary_layer.mjs.
// Derived from XFOIL by Mark Drela. GNU General Public License v2 or later.
// https://web.mit.edu/drela/Public/web/xfoil/
// GPL-2.0-or-later. Minimal checked array ABI for the generated BL equations.
// No libc, allocation runtime, fast math, or host-dependent libm is linked.
using u32 = unsigned;
static_assert(sizeof(double) == 8 && sizeof(u32) == 4);
[[noreturn]] void unsupported() { __builtin_trap(); }
double nan_value() { return __builtin_nan(""); }
bool finite(double v) { return __builtin_isfinite(v); }
bool truth(double v) { return v != 0.0 && v == v; }
double js_min(double a,double b) { if(a!=a||b!=b)return nan_value(); if(a==0&&b==0)return __builtin_signbit(a)?a:b; return a<b?a:b; }
double js_max(double a,double b) { if(a!=a||b!=b)return nan_value(); if(a==0&&b==0)return __builtin_signbit(a)?b:a; return a>b?a:b; }
extern "C" {
// These imports are the deterministic host Math functions. Declaring them pure
// allows reuse of identical calls without changing arithmetic or rounding.
__attribute__((const, import_module("math"),import_name("pow"))) double host_pow(double,double);
__attribute__((const, import_module("math"),import_name("exp"))) double host_exp(double);
__attribute__((const, import_module("math"),import_name("log"))) double host_log(double);
__attribute__((const, import_module("math"),import_name("log10"))) double host_log10(double);
__attribute__((const, import_module("math"),import_name("sin"))) double host_sin(double);
__attribute__((const, import_module("math"),import_name("cos"))) double host_cos(double);
__attribute__((const, import_module("math"),import_name("atan"))) double host_atan(double);
__attribute__((const, import_module("math"),import_name("tanh"))) double host_tanh(double);
__attribute__((import_module("diagnostics"),import_name("log"))) void record_log(u32,double,double,double);
}
u32 arena_cursor,arena_limit;
void arena_begin(u32 begin,u32 end) { arena_cursor=begin; arena_limit=end; }
u32 allocate(u32 bytes) {
  if(bytes>arena_limit-arena_cursor)unsupported();
  u32 start=arena_cursor; arena_cursor+=(bytes+7)&~7u;
  if(arena_cursor>arena_limit)unsupported();
  for(u32 i=start;i<arena_cursor;i++) *((unsigned char*)i)=0;
  return start;
}
u32 size_value(double n) { if(!finite(n)||n<0||n>20000000||n!=(u32)n)unsupported();return(u32)n; }
// An existing array's length already bounds its indices. Check this range
// before converting, which also rejects NaN/infinity without another scan of
// their exponent bits. Retain the integer check: fractional JS properties
// must never silently address a neighboring numeric element.
u32 array_index(double n,u32 length) {
  if(!(n>=0 && n<length))unsupported();
  u32 i=(u32)n;
  if(n!=i)unsupported();
  return i;
}
void advance_index(u32& i) { if(i==0xffffffffu)unsupported();++i; }
struct Descriptor { u32 pointer,length; };
template<int R> struct Vec;
template<> struct Vec<1> {
  u32 pointer,length;
  explicit operator bool()const{return pointer!=0;}
  double& at(u32 i)const { if(!pointer||i>=length)unsupported();return ((double*)pointer)[i]; }
  double& operator[](double n)const { if(!pointer)unsupported();return ((double*)pointer)[array_index(n,length)]; }
  double get_or(double n,double fallback)const { return finite(n)&&n>=0&&n<length?(*this)[n]:fallback; }
  void fill(double value,double begin=0,double end=-1)const { u32 a=size_value(begin),b=end<0?length:size_value(end);if(a>length||b>length)unsupported();for(u32 i=a;i<b;i++)(*this)[i]=value; }
  Vec<1> subarray(double begin,double end)const {u32 a=size_value(begin),b=size_value(end);if(a>b||b>length)unsupported();return {pointer+8*a,b-a};}
};
template<int R> struct Vec {
  u32 pointer,length;
  explicit operator bool()const{return pointer!=0;}
  Vec<R-1>& at(u32 i)const {if(!pointer||i>=length)unsupported();return ((Vec<R-1>*)pointer)[i];}
  Vec<R-1>& operator[](double n)const {if(!pointer)unsupported();return ((Vec<R-1>*)pointer)[array_index(n,length)];}
};
using V1=Vec<1>;using V2=Vec<2>;using V3=Vec<3>;
template<int R> Vec<R> make_vec(double n) {u32 size=size_value(n);return {allocate(size*8),size};}
V2 make_matrix(double rows,double columns) {auto m=make_vec<2>(rows);for(u32 i=0;i<m.length;i++)m[i]=make_vec<1>(columns);return m;}
template<class... T> V2 make_rows(T... values) { V1 rows[]{values...};auto m=make_vec<2>(sizeof...(values));u32 i=0;for(auto row:rows)m[i++]=row;return m;}
struct Sources {u32 pointer;explicit operator bool()const{return pointer!=0;}V1& arr(u32 i)const{if(!pointer||i>4)unsupported();return((V1*)pointer)[i];}double& count()const{if(!pointer)unsupported();return*(double*)(pointer+40);} };
Sources make_sources(){return {allocate(48)};}
struct Slot { double value;u32 present,dirty; };
struct Context {
  Slot* scalar;Descriptor* array;Sources* source;
  double get(u32 index)const{return scalar[index].value;}
  double put(u32 index,double value)const{scalar[index]={value,1,1};return value;}
  bool has(u32 index)const{return scalar[index].present!=0;}
  template<int R> Vec<R>& arr(u32 index)const{return *(Vec<R>*)(array+index);}
  Sources& sources(u32 index)const{return source[index];}
};
// Match spline.js's NMAX guard; out-of-range accesses additionally trap below
// the JS bridge so they can be retried with the original error behavior.
void check_n(double n) {if(n>2500)unsupported();}

struct Record {
double ax = nan_value();
double axA1 = nan_value();
double axA2 = nan_value();
double axHk = nan_value();
double axHk1 = nan_value();
double axHk2 = nan_value();
double axRt = nan_value();
double axRt1 = nan_value();
double axRt2 = nan_value();
double axT1 = nan_value();
double axT2 = nan_value();
double axTh = nan_value();
double cf = nan_value();
double cfHk = nan_value();
double cfMsq = nan_value();
double cfRt = nan_value();
double di = nan_value();
double diCf = nan_value();
double diHk = nan_value();
double diHs = nan_value();
double diRt = nan_value();
double diSt = nan_value();
double diUs = nan_value();
double dsi = nan_value();
double hc = nan_value();
double hcHk = nan_value();
double hcMsq = nan_value();
double hk = nan_value();
double hkH = nan_value();
double hkMsq = nan_value();
double hs = nan_value();
double hsHk = nan_value();
double hsMsq = nan_value();
double hsRt = nan_value();
double length = nan_value();
double mCls = nan_value();
double rCls = nan_value();
double thi = nan_value();
};
void fn_syncComToVars(Context ctx, double which);
void fn_syncVarsToCom(Context ctx, double which);
void fn_copyCom(Context ctx, double from, double to);
void fn_trchek(Context ctxIn);
Record fn_axset(double hk1, double t1, double rt1, double a1, double hk2, double t2, double rt2, double a2, double acrit, double idampv);
void fn_trchek2(Context ctx);
void fn_blsys(Context ctxIn);
void fn_tesys(double cte, double tte, double dte, Context ctxIn);
void fn_blprv(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctxIn);
void fn_blprvKin(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctxIn);
void fn_predictBlStation(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctx);
void fn_blkin(Context ctxIn);
void fn_blvar(double ityp, Context ctxIn);
void fn_blmid(double ityp, Context ctxIn);
void fn_trdif(Context ctxIn);
void fn_bldif(double ityp, Context ctxIn);
Record fn_dampl(double hk, double th, double rt);
Record fn_dampl2(double hk, double th, double rt);
Record fn_hkin(double h, double msq);
Record fn_dil(double hk, double rt);
Record fn_dilw(double hk, double rt);
Record fn_hsl(double hk, double rt, double msq);
Record fn_cfl(double hk, double rt, double msq);
Record fn_dit(double hs, double us, double cf, double st);
Record fn_hst(double hk, double rt, double msq);
Record fn_cft(double hk, double rt, double msq, double cffac = 1.0);
Record fn_hct(double hk, double msq);
void fn_prepareBlStation(Context ctx, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei);
double fn_dslim(Context ctx, double dstr, double thet, double uedg, double msq, double hklim);
Record fn_mrcl(Context ctx, double cls);
void fn_comset(Context ctx);
void fn_iblsys(Context ctx);
void fn_mrchue(Context ctx);
void fn_mrchdu(Context ctx);
void fn_assembleMassRows(Context ctx, double iv, V1 d1M, V1 u1M, V1 d2M, V1 u2M, V1 ule1M, V1 ule2M, double xiUle1, double xiUle2);
void fn_setbl(Context ctx);
void fn_xifset(Context ctx, double is);
Sources fn_velocitySources(Context ctx, double corrections);
void fn_ueset(Context ctx);
void fn_gauss1(double n, V2 z, V1 r);
void fn_trisol(V1 a, V1 b, V1 c, V1 d, double kk, double dOffset = 0.0);
void fn_splindOffset(V1 x, V1 xs, V1 s, double n, double xs1, double xs2, double offset);
void fn_splind(V1 x, V1 xs, V1 s, double n, double xs1, double xs2);
double fn_seval(double ss, V1 x, V1 xs, V1 s, double n);
double fn_deval(double ss, V1 x, V1 xs, V1 s, double n);
double fn_sinvrt(double si, double xi, V1 x, V1 xs, V1 s, double n);
// js/xblsys.js:124 (syncComToVars)
void fn_syncComToVars(Context ctx, double which) {
((void)0);
if (truth((which == 1.0))) {
V1 data = ctx.arr<1>(3);
ctx.put(228, data[1.0]);
ctx.put(200, data[2.0]);
ctx.put(188, data[3.0]);
ctx.put(52, data[4.0]);
ctx.put(180, data[5.0]);
ctx.put(1, data[6.0]);
ctx.put(202, data[7.0]);
ctx.put(201, data[8.0]);
ctx.put(81, data[9.0]);
ctx.put(91, data[10.0]);
ctx.put(93, data[11.0]);
ctx.put(92, data[12.0]);
ctx.put(139, data[13.0]);
ctx.put(141, data[14.0]);
ctx.put(140, data[15.0]);
ctx.put(155, data[16.0]);
ctx.put(157, data[17.0]);
ctx.put(156, data[18.0]);
ctx.put(218, data[19.0]);
ctx.put(221, data[20.0]);
ctx.put(219, data[21.0]);
ctx.put(220, data[22.0]);
ctx.put(108, data[23.0]);
ctx.put(112, data[24.0]);
ctx.put(111, data[25.0]);
ctx.put(109, data[26.0]);
ctx.put(110, data[27.0]);
ctx.put(118, data[28.0]);
ctx.put(123, data[29.0]);
ctx.put(122, data[30.0]);
ctx.put(119, data[31.0]);
ctx.put(120, data[32.0]);
ctx.put(121, data[33.0]);
ctx.put(98, data[34.0]);
ctx.put(102, data[35.0]);
ctx.put(101, data[36.0]);
ctx.put(99, data[37.0]);
ctx.put(100, data[38.0]);
ctx.put(170, data[39.0]);
ctx.put(174, data[40.0]);
ctx.put(173, data[41.0]);
ctx.put(171, data[42.0]);
ctx.put(172, data[43.0]);
ctx.put(5, data[44.0]);
ctx.put(10, data[45.0]);
ctx.put(9, data[46.0]);
ctx.put(6, data[47.0]);
ctx.put(7, data[48.0]);
ctx.put(8, data[49.0]);
ctx.put(64, data[50.0]);
ctx.put(70, data[51.0]);
ctx.put(69, data[52.0]);
ctx.put(65, data[53.0]);
ctx.put(68, data[54.0]);
ctx.put(66, data[55.0]);
ctx.put(67, data[56.0]);
ctx.put(206, data[57.0]);
ctx.put(211, data[58.0]);
ctx.put(210, data[59.0]);
ctx.put(207, data[60.0]);
ctx.put(208, data[61.0]);
ctx.put(209, data[62.0]);
ctx.put(36, data[63.0]);
ctx.put(41, data[64.0]);
ctx.put(40, data[65.0]);
ctx.put(37, data[66.0]);
ctx.put(38, data[67.0]);
ctx.put(39, data[68.0]);
ctx.put(54, data[69.0]);
ctx.put(58, data[70.0]);
ctx.put(57, data[71.0]);
ctx.put(55, data[72.0]);
ctx.put(56, data[73.0]);
} else {
V1 data = ctx.arr<1>(4);
ctx.put(229, data[1.0]);
ctx.put(203, data[2.0]);
ctx.put(189, data[3.0]);
ctx.put(53, data[4.0]);
ctx.put(181, data[5.0]);
ctx.put(2, data[6.0]);
ctx.put(205, data[7.0]);
ctx.put(204, data[8.0]);
ctx.put(82, data[9.0]);
ctx.put(94, data[10.0]);
ctx.put(96, data[11.0]);
ctx.put(95, data[12.0]);
ctx.put(142, data[13.0]);
ctx.put(144, data[14.0]);
ctx.put(143, data[15.0]);
ctx.put(158, data[16.0]);
ctx.put(160, data[17.0]);
ctx.put(159, data[18.0]);
ctx.put(222, data[19.0]);
ctx.put(225, data[20.0]);
ctx.put(223, data[21.0]);
ctx.put(224, data[22.0]);
ctx.put(113, data[23.0]);
ctx.put(117, data[24.0]);
ctx.put(116, data[25.0]);
ctx.put(114, data[26.0]);
ctx.put(115, data[27.0]);
ctx.put(124, data[28.0]);
ctx.put(129, data[29.0]);
ctx.put(128, data[30.0]);
ctx.put(125, data[31.0]);
ctx.put(126, data[32.0]);
ctx.put(127, data[33.0]);
ctx.put(103, data[34.0]);
ctx.put(107, data[35.0]);
ctx.put(106, data[36.0]);
ctx.put(104, data[37.0]);
ctx.put(105, data[38.0]);
ctx.put(175, data[39.0]);
ctx.put(179, data[40.0]);
ctx.put(178, data[41.0]);
ctx.put(176, data[42.0]);
ctx.put(177, data[43.0]);
ctx.put(11, data[44.0]);
ctx.put(19, data[45.0]);
ctx.put(18, data[46.0]);
ctx.put(12, data[47.0]);
ctx.put(15, data[48.0]);
ctx.put(16, data[49.0]);
ctx.put(71, data[50.0]);
ctx.put(77, data[51.0]);
ctx.put(76, data[52.0]);
ctx.put(72, data[53.0]);
ctx.put(75, data[54.0]);
ctx.put(73, data[55.0]);
ctx.put(74, data[56.0]);
ctx.put(212, data[57.0]);
ctx.put(217, data[58.0]);
ctx.put(216, data[59.0]);
ctx.put(213, data[60.0]);
ctx.put(214, data[61.0]);
ctx.put(215, data[62.0]);
ctx.put(42, data[63.0]);
ctx.put(47, data[64.0]);
ctx.put(46, data[65.0]);
ctx.put(43, data[66.0]);
ctx.put(44, data[67.0]);
ctx.put(45, data[68.0]);
ctx.put(59, data[69.0]);
ctx.put(63, data[70.0]);
ctx.put(62, data[71.0]);
ctx.put(60, data[72.0]);
ctx.put(61, data[73.0]);
}
}

// js/xblsys.js:280 (syncVarsToCom)
void fn_syncVarsToCom(Context ctx, double which) {
((void)0);
if (truth((which == 1.0))) {
V1 data = ctx.arr<1>(3);
(data[1.0] = ctx.get(228));
(data[2.0] = ctx.get(200));
(data[3.0] = ctx.get(188));
(data[4.0] = ctx.get(52));
(data[5.0] = ctx.get(180));
(data[6.0] = ctx.get(1));
(data[7.0] = ctx.get(202));
(data[8.0] = ctx.get(201));
(data[9.0] = ctx.get(81));
(data[10.0] = ctx.get(91));
(data[11.0] = ctx.get(93));
(data[12.0] = ctx.get(92));
(data[13.0] = ctx.get(139));
(data[14.0] = ctx.get(141));
(data[15.0] = ctx.get(140));
(data[16.0] = ctx.get(155));
(data[17.0] = ctx.get(157));
(data[18.0] = ctx.get(156));
(data[19.0] = ctx.get(218));
(data[20.0] = ctx.get(221));
(data[21.0] = ctx.get(219));
(data[22.0] = ctx.get(220));
(data[23.0] = ctx.get(108));
(data[24.0] = ctx.get(112));
(data[25.0] = ctx.get(111));
(data[26.0] = ctx.get(109));
(data[27.0] = ctx.get(110));
(data[28.0] = ctx.get(118));
(data[29.0] = ctx.get(123));
(data[30.0] = ctx.get(122));
(data[31.0] = ctx.get(119));
(data[32.0] = ctx.get(120));
(data[33.0] = ctx.get(121));
(data[34.0] = ctx.get(98));
(data[35.0] = ctx.get(102));
(data[36.0] = ctx.get(101));
(data[37.0] = ctx.get(99));
(data[38.0] = ctx.get(100));
(data[39.0] = ctx.get(170));
(data[40.0] = ctx.get(174));
(data[41.0] = ctx.get(173));
(data[42.0] = ctx.get(171));
(data[43.0] = ctx.get(172));
(data[44.0] = ctx.get(5));
(data[45.0] = ctx.get(10));
(data[46.0] = ctx.get(9));
(data[47.0] = ctx.get(6));
(data[48.0] = ctx.get(7));
(data[49.0] = ctx.get(8));
(data[50.0] = ctx.get(64));
(data[51.0] = ctx.get(70));
(data[52.0] = ctx.get(69));
(data[53.0] = ctx.get(65));
(data[54.0] = ctx.get(68));
(data[55.0] = ctx.get(66));
(data[56.0] = ctx.get(67));
(data[57.0] = ctx.get(206));
(data[58.0] = ctx.get(211));
(data[59.0] = ctx.get(210));
(data[60.0] = ctx.get(207));
(data[61.0] = ctx.get(208));
(data[62.0] = ctx.get(209));
(data[63.0] = ctx.get(36));
(data[64.0] = ctx.get(41));
(data[65.0] = ctx.get(40));
(data[66.0] = ctx.get(37));
(data[67.0] = ctx.get(38));
(data[68.0] = ctx.get(39));
(data[69.0] = ctx.get(54));
(data[70.0] = ctx.get(58));
(data[71.0] = ctx.get(57));
(data[72.0] = ctx.get(55));
(data[73.0] = ctx.get(56));
} else {
V1 data = ctx.arr<1>(4);
(data[1.0] = ctx.get(229));
(data[2.0] = ctx.get(203));
(data[3.0] = ctx.get(189));
(data[4.0] = ctx.get(53));
(data[5.0] = ctx.get(181));
(data[6.0] = ctx.get(2));
(data[7.0] = ctx.get(205));
(data[8.0] = ctx.get(204));
(data[9.0] = ctx.get(82));
(data[10.0] = ctx.get(94));
(data[11.0] = ctx.get(96));
(data[12.0] = ctx.get(95));
(data[13.0] = ctx.get(142));
(data[14.0] = ctx.get(144));
(data[15.0] = ctx.get(143));
(data[16.0] = ctx.get(158));
(data[17.0] = ctx.get(160));
(data[18.0] = ctx.get(159));
(data[19.0] = ctx.get(222));
(data[20.0] = ctx.get(225));
(data[21.0] = ctx.get(223));
(data[22.0] = ctx.get(224));
(data[23.0] = ctx.get(113));
(data[24.0] = ctx.get(117));
(data[25.0] = ctx.get(116));
(data[26.0] = ctx.get(114));
(data[27.0] = ctx.get(115));
(data[28.0] = ctx.get(124));
(data[29.0] = ctx.get(129));
(data[30.0] = ctx.get(128));
(data[31.0] = ctx.get(125));
(data[32.0] = ctx.get(126));
(data[33.0] = ctx.get(127));
(data[34.0] = ctx.get(103));
(data[35.0] = ctx.get(107));
(data[36.0] = ctx.get(106));
(data[37.0] = ctx.get(104));
(data[38.0] = ctx.get(105));
(data[39.0] = ctx.get(175));
(data[40.0] = ctx.get(179));
(data[41.0] = ctx.get(178));
(data[42.0] = ctx.get(176));
(data[43.0] = ctx.get(177));
(data[44.0] = ctx.get(11));
(data[45.0] = ctx.get(19));
(data[46.0] = ctx.get(18));
(data[47.0] = ctx.get(12));
(data[48.0] = ctx.get(15));
(data[49.0] = ctx.get(16));
(data[50.0] = ctx.get(71));
(data[51.0] = ctx.get(77));
(data[52.0] = ctx.get(76));
(data[53.0] = ctx.get(72));
(data[54.0] = ctx.get(75));
(data[55.0] = ctx.get(73));
(data[56.0] = ctx.get(74));
(data[57.0] = ctx.get(212));
(data[58.0] = ctx.get(217));
(data[59.0] = ctx.get(216));
(data[60.0] = ctx.get(213));
(data[61.0] = ctx.get(214));
(data[62.0] = ctx.get(215));
(data[63.0] = ctx.get(42));
(data[64.0] = ctx.get(47));
(data[65.0] = ctx.get(46));
(data[66.0] = ctx.get(43));
(data[67.0] = ctx.get(44));
(data[68.0] = ctx.get(45));
(data[69.0] = ctx.get(59));
(data[70.0] = ctx.get(63));
(data[71.0] = ctx.get(62));
(data[72.0] = ctx.get(60));
(data[73.0] = ctx.get(61));
}
}

// js/xblsys.js:436 (copyCom)
void fn_copyCom(Context ctx, double from, double to) {
((void)0);
V1 src = (truth((from == 1.0)) ? ctx.arr<1>(3) : ctx.arr<1>(4));
V1 dst = (truth((to == 1.0)) ? ctx.arr<1>(3) : ctx.arr<1>(4));
for (u32 i = 1u; truth((double(i) <= 73.0)); advance_index(i)) {
(dst.at(i) = src.at(i));
}
fn_syncComToVars(ctx, to);
}

// js/xblsys.js:447 (trchek)
void fn_trchek(Context ctxIn) {
Context ctx = ctxIn;
fn_trchek2(ctx);
}

// js/xblsys.js:453 (axset)
Record fn_axset(double hk1, double t1, double rt1, double a1, double hk2, double t2, double rt2, double a2, double acrit, double idampv) {
double ax1 = nan_value();
double ax1Hk1 = nan_value();
double ax1T1 = nan_value();
double ax1Rt1 = nan_value();
double ax2 = nan_value();
double ax2Hk2 = nan_value();
double ax2T2 = nan_value();
double ax2Rt2 = nan_value();
if (truth((idampv == 0.0))) {
([&]() { auto value = fn_dampl(hk1, t1, rt1); ax1 = value.ax; ax1Hk1 = value.axHk; ax1T1 = value.axTh; ax1Rt1 = value.axRt; return value; }());
([&]() { auto value = fn_dampl(hk2, t2, rt2); ax2 = value.ax; ax2Hk2 = value.axHk; ax2T2 = value.axTh; ax2Rt2 = value.axRt; return value; }());
} else {
([&]() { auto value = fn_dampl2(hk1, t1, rt1); ax1 = value.ax; ax1Hk1 = value.axHk; ax1T1 = value.axTh; ax1Rt1 = value.axRt; return value; }());
([&]() { auto value = fn_dampl2(hk2, t2, rt2); ax2 = value.ax; ax2Hk2 = value.axHk; ax2T2 = value.axTh; ax2Rt2 = value.axRt; return value; }());
}
double axsq = (0.5 * (host_pow(ax1, 2.0) + host_pow(ax2, 2.0)));
double axa = nan_value();
double axaAx1 = nan_value();
double axaAx2 = nan_value();
if (truth((axsq <= 0.0))) {
(axa = 0.0);
(axaAx1 = 0.0);
(axaAx2 = 0.0);
} else {
(axa = __builtin_sqrt(axsq));
(axaAx1 = ((0.5 * ax1) / axa));
(axaAx2 = ((0.5 * ax2) / axa));
}
double arg = js_min((20.0 * (acrit - (0.5 * (a1 + a2)))), 20.0);
double exn = nan_value();
double exnA1 = nan_value();
double exnA2 = nan_value();
if (truth((arg <= 0.0))) {
(exn = 1.0);
(exnA1 = 0.0);
(exnA2 = 0.0);
} else {
(exn = host_exp((-arg)));
(exnA1 = ((20.0 * 0.5) * exn));
(exnA2 = ((20.0 * 0.5) * exn));
}
double dax = ((exn * 0.002) / (t1 + t2));
double daxA1 = ((exnA1 * 0.002) / (t1 + t2));
double daxA2 = ((exnA2 * 0.002) / (t1 + t2));
double daxT1 = ((-dax) / (t1 + t2));
double daxT2 = ((-dax) / (t1 + t2));
double ax = (axa + dax);
double axHk1 = (axaAx1 * ax1Hk1);
double axT1 = ((axaAx1 * ax1T1) + daxT1);
double axRt1 = (axaAx1 * ax1Rt1);
double axA1 = daxA1;
double axHk2 = (axaAx2 * ax2Hk2);
double axT2 = ((axaAx2 * ax2T2) + daxT2);
double axRt2 = (axaAx2 * ax2Rt2);
double axA2 = daxA2;
return ([&]() { Record value{}; value.ax = ax; value.axHk1 = axHk1; value.axT1 = axT1; value.axRt1 = axRt1; value.axA1 = axA1; value.axHk2 = axHk2; value.axT2 = axT2; value.axRt2 = axRt2; value.axA2 = axA2; return value; }());
}

// js/xblsys.js:529 (trchek2)
void fn_trchek2(Context ctx) {
((void)0);
double daeps = 5.0e-5;
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(ctx.arr<1>(2).at(icom) = ctx.arr<1>(4).at(icom));
}
Record axOut = fn_axset(ctx.get(108), ctx.get(188), ctx.get(170), ctx.get(1), ctx.get(113), ctx.get(189), ctx.get(175), ctx.get(2), ctx.get(0), ctx.get(134));
double ax = axOut.ax;
ctx.put(2, (ctx.get(1) + (ax * (ctx.get(229) - ctx.get(228)))));
double wf2 = 0.0;
double wf2A1 = 0.0;
double wf2A2 = 0.0;
double wf2X1 = 0.0;
double wf2X2 = 0.0;
double wf2Xf = 0.0;
double wf1 = 0.0;
double wf1A1 = 0.0;
double wf1A2 = 0.0;
double wf1X1 = 0.0;
double wf1X2 = 0.0;
double wf1Xf = 0.0;
double xt = 0.0;
double tt = 0.0;
double dt = 0.0;
double ut = 0.0;
double xtA2 = 0.0;
double ttA2 = 0.0;
double dtA2 = 0.0;
double utA2 = 0.0;
double hkt = 0.0;
double hktTt = 0.0;
double hktDt = 0.0;
double hktUt = 0.0;
double hktMs = 0.0;
double rtt = 0.0;
double rttTt = 0.0;
double rttUt = 0.0;
double rttMs = 0.0;
double rttRe = 0.0;
double amplt = 0.0;
double ampltA2 = 0.0;
for (u32 itam = 1u; truth((double(itam) <= 30.0)); advance_index(itam)) {
double sfa = nan_value();
double sfaA1 = nan_value();
double sfaA2 = nan_value();
if (truth((ctx.get(2) <= ctx.get(0)))) {
(amplt = ctx.get(2));
(ampltA2 = 1.0);
(sfa = 1.0);
(sfaA1 = 0.0);
(sfaA2 = 0.0);
} else {
(amplt = ctx.get(0));
(ampltA2 = 0.0);
(sfa = ((amplt - ctx.get(1)) / (ctx.get(2) - ctx.get(1))));
(sfaA1 = ((sfa - 1.0) / (ctx.get(2) - ctx.get(1))));
(sfaA2 = ((-sfa) / (ctx.get(2) - ctx.get(1))));
}
double sfx = nan_value();
double sfxX1 = nan_value();
double sfxX2 = nan_value();
double sfxXf = nan_value();
double xifEps = (1.0e-12 * js_max(1.0, __builtin_fabs(ctx.get(229))));
double xif = ((truth((ctx.get(230) >= (ctx.get(229) - xifEps))) && truth((ctx.get(230) <= (ctx.get(229) + xifEps)))) ? ctx.get(229) : ctx.get(230));
if (truth((xif < ctx.get(229)))) {
(sfx = ((xif - ctx.get(228)) / (ctx.get(229) - ctx.get(228))));
(sfxX1 = ((sfx - 1.0) / (ctx.get(229) - ctx.get(228))));
(sfxX2 = ((-sfx) / (ctx.get(229) - ctx.get(228))));
(sfxXf = (1.0 / (ctx.get(229) - ctx.get(228))));
} else {
(sfx = 1.0);
(sfxX1 = 0.0);
(sfxX2 = 0.0);
(sfxXf = 0.0);
}
if (truth((sfa < sfx))) {
(wf2 = sfa);
(wf2A1 = sfaA1);
(wf2A2 = sfaA2);
(wf2X1 = 0.0);
(wf2X2 = 0.0);
(wf2Xf = 0.0);
} else {
(wf2 = sfx);
(wf2A1 = 0.0);
(wf2A2 = 0.0);
(wf2X1 = sfxX1);
(wf2X2 = sfxX2);
(wf2Xf = sfxXf);
}
(wf1 = (1.0 - wf2));
(wf1A1 = (-wf2A1));
(wf1A2 = (-wf2A2));
(wf1X1 = (-wf2X1));
(wf1X2 = (-wf2X2));
(wf1Xf = (-wf2Xf));
(xt = ((ctx.get(228) * wf1) + (ctx.get(229) * wf2)));
(tt = ((ctx.get(188) * wf1) + (ctx.get(189) * wf2)));
(dt = ((ctx.get(52) * wf1) + (ctx.get(53) * wf2)));
(ut = ((ctx.get(200) * wf1) + (ctx.get(203) * wf2)));
(xtA2 = ((ctx.get(228) * wf1A2) + (ctx.get(229) * wf2A2)));
(ttA2 = ((ctx.get(188) * wf1A2) + (ctx.get(189) * wf2A2)));
(dtA2 = ((ctx.get(52) * wf1A2) + (ctx.get(53) * wf2A2)));
(utA2 = ((ctx.get(200) * wf1A2) + (ctx.get(203) * wf2A2)));
ctx.put(229, xt);
ctx.put(189, tt);
ctx.put(53, dt);
ctx.put(203, ut);
fn_blkin(ctx);
(hkt = ctx.get(113));
(hktTt = ctx.get(116));
(hktDt = ctx.get(114));
(hktUt = ctx.get(117));
(hktMs = ctx.get(115));
(rtt = ctx.get(175));
(rttTt = ctx.get(178));
(rttUt = ctx.get(179));
(rttMs = ctx.get(176));
(rttRe = ctx.get(177));
double amsave = ctx.get(2);
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(ctx.arr<1>(4).at(icom) = ctx.arr<1>(2).at(icom));
}
fn_syncComToVars(ctx, 2.0);
ctx.put(2, amsave);
(axOut = fn_axset(ctx.get(108), ctx.get(188), ctx.get(170), ctx.get(1), hkt, tt, rtt, amplt, ctx.get(0), ctx.get(134)));
(ax = axOut.ax);
if (truth((ax <= 0.0))) break;
double axA2 = (((((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttA2) + ((axOut.axHk2 * hktDt) * dtA2)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utA2)) + (axOut.axA2 * ampltA2));
double res = ((ctx.get(2) - ctx.get(1)) - (ax * (ctx.get(229) - ctx.get(228))));
double resA2 = (1.0 - (axA2 * (ctx.get(229) - ctx.get(228))));
double da2 = ((-res) / resA2);
double rlx = 1.0;
double dxt = (xtA2 * da2);
if (truth(((rlx * __builtin_fabs((dxt / (ctx.get(229) - ctx.get(228))))) > 0.05))) {
(rlx = (0.05 * __builtin_fabs(((ctx.get(229) - ctx.get(228)) / dxt))));
}
if (truth(((rlx * __builtin_fabs(da2)) > 1.0))) {
(rlx = (1.0 * __builtin_fabs((1.0 / da2))));
}
if (truth((__builtin_fabs(da2) < daeps))) break;
if (((truth((ctx.get(2) > ctx.get(0))) && truth(((ctx.get(2) + (rlx * da2)) < ctx.get(0)))) || (truth((ctx.get(2) < ctx.get(0))) && truth(((ctx.get(2) + (rlx * da2)) > ctx.get(0)))))) {
ctx.put(2, ctx.get(0));
} else {
ctx.put(2, (ctx.get(2) + (rlx * da2)));
}
}
ctx.put(232, xt);
ctx.put(235, xtA2);
ctx.put(197, (ctx.get(2) >= ctx.get(0)));
double xifEpsFinal = (1.0e-12 * js_max(1.0, __builtin_fabs(ctx.get(229))));
double xifFinal = ((truth((ctx.get(230) >= (ctx.get(229) - xifEpsFinal))) && truth((ctx.get(230) <= (ctx.get(229) + xifEpsFinal)))) ? ctx.get(229) : ctx.get(230));
ctx.put(196, ([&]() -> double { double value = (xifFinal > ctx.get(228)); return truth(value) ? (xifFinal <= ctx.get(229)) : value; }()));
ctx.put(195, ([&]() -> double { double value = ctx.get(196); return truth(value) ? value : ctx.get(197); }()));
if (truth((!truth(ctx.get(195))))) return;
if ((truth(ctx.get(197)) && truth(ctx.get(196)))) {
ctx.put(196, (xifFinal < xt));
ctx.put(197, (xifFinal >= xt));
}
if (truth(ctx.get(196))) {
ctx.put(232, xifFinal);
ctx.put(234, 0.0);
ctx.put(244, 0.0);
ctx.put(240, 0.0);
ctx.put(236, 0.0);
ctx.put(242, 0.0);
ctx.put(245, 0.0);
ctx.put(241, 0.0);
ctx.put(237, 0.0);
ctx.put(243, 0.0);
ctx.put(238, 0.0);
ctx.put(239, 0.0);
ctx.put(246, 1.0);
return;
}
ctx.put(232, ((ctx.get(228) * wf1) + (ctx.get(229) * wf2)));
ctx.put(235, xtA2);
(tt = ((ctx.get(188) * wf1) + (ctx.get(189) * wf2)));
(dt = ((ctx.get(52) * wf1) + (ctx.get(53) * wf2)));
(ut = ((ctx.get(200) * wf1) + (ctx.get(203) * wf2)));
ctx.put(244, wf1);
double ttT1 = wf1;
double dtD1 = wf1;
double utU1 = wf1;
ctx.put(245, wf2);
double ttT2 = wf2;
double dtD2 = wf2;
double utU2 = wf2;
ctx.put(234, ((ctx.get(228) * wf1A1) + (ctx.get(229) * wf2A1)));
double ttA1 = ((ctx.get(188) * wf1A1) + (ctx.get(189) * wf2A1));
double dtA1 = ((ctx.get(52) * wf1A1) + (ctx.get(53) * wf2A1));
double utA1 = ((ctx.get(200) * wf1A1) + (ctx.get(203) * wf2A1));
ctx.put(244, (((ctx.get(228) * wf1X1) + (ctx.get(229) * wf2X1)) + ctx.get(244)));
double ttX1 = ((ctx.get(188) * wf1X1) + (ctx.get(189) * wf2X1));
double dtX1 = ((ctx.get(52) * wf1X1) + (ctx.get(53) * wf2X1));
double utX1 = ((ctx.get(200) * wf1X1) + (ctx.get(203) * wf2X1));
ctx.put(245, (((ctx.get(228) * wf1X2) + (ctx.get(229) * wf2X2)) + ctx.get(245)));
double ttX2 = ((ctx.get(188) * wf1X2) + (ctx.get(189) * wf2X2));
double dtX2 = ((ctx.get(52) * wf1X2) + (ctx.get(53) * wf2X2));
double utX2 = ((ctx.get(200) * wf1X2) + (ctx.get(203) * wf2X2));
ctx.put(246, ((ctx.get(228) * wf1Xf) + (ctx.get(229) * wf2Xf)));
double ttXf = ((ctx.get(188) * wf1Xf) + (ctx.get(189) * wf2Xf));
double dtXf = ((ctx.get(52) * wf1Xf) + (ctx.get(53) * wf2Xf));
double utXf = ((ctx.get(200) * wf1Xf) + (ctx.get(203) * wf2Xf));
double axT1 = ((((axOut.axHk1 * ctx.get(111)) + axOut.axT1) + (axOut.axRt1 * ctx.get(173))) + ((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttT1));
double axD1 = ((axOut.axHk1 * ctx.get(109)) + ((axOut.axHk2 * hktDt) * dtD1));
double axU1 = (((axOut.axHk1 * ctx.get(112)) + (axOut.axRt1 * ctx.get(174))) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utU1));
double axA1 = (((axOut.axA1 + ((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttA1)) + ((axOut.axHk2 * hktDt) * dtA1)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utA1));
double axX1 = ((((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttX1) + ((axOut.axHk2 * hktDt) * dtX1)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utX1));
double axT2 = ((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttT2);
double axD2 = ((axOut.axHk2 * hktDt) * dtD2);
double axU2 = (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utU2);
double axA2 = ((((axOut.axA2 * ampltA2) + ((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttA2)) + ((axOut.axHk2 * hktDt) * dtA2)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utA2));
double axX2 = ((((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttX2) + ((axOut.axHk2 * hktDt) * dtX2)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utX2));
double axXf = ((((((axOut.axHk2 * hktTt) + axOut.axT2) + (axOut.axRt2 * rttTt)) * ttXf) + ((axOut.axHk2 * hktDt) * dtXf)) + (((axOut.axHk2 * hktUt) + (axOut.axRt2 * rttUt)) * utXf));
double axMs = ((((axOut.axHk2 * hktMs) + (axOut.axRt2 * rttMs)) + (axOut.axHk1 * ctx.get(110))) + (axOut.axRt1 * ctx.get(171)));
double axRe = ((axOut.axRt2 * rttRe) + (axOut.axRt1 * ctx.get(172)));
double zAx = (-(ctx.get(229) - ctx.get(228)));
double zA1 = ((zAx * axA1) - 1.0);
double zT1 = (zAx * axT1);
double zD1 = (zAx * axD1);
double zU1 = (zAx * axU1);
double zX1 = ((zAx * axX1) + ax);
double zA2 = ((zAx * axA2) + 1.0);
double zT2 = (zAx * axT2);
double zD2 = (zAx * axD2);
double zU2 = (zAx * axU2);
double zX2 = ((zAx * axX2) - ax);
double zXf = (zAx * axXf);
double zMs = (zAx * axMs);
double zRe = (zAx * axRe);
ctx.put(234, (ctx.get(234) - ((ctx.get(235) / zA2) * zA1)));
ctx.put(240, ((-(ctx.get(235) / zA2)) * zT1));
ctx.put(236, ((-(ctx.get(235) / zA2)) * zD1));
ctx.put(242, ((-(ctx.get(235) / zA2)) * zU1));
ctx.put(244, (ctx.get(244) - ((ctx.get(235) / zA2) * zX1)));
ctx.put(241, ((-(ctx.get(235) / zA2)) * zT2));
ctx.put(237, ((-(ctx.get(235) / zA2)) * zD2));
ctx.put(243, ((-(ctx.get(235) / zA2)) * zU2));
ctx.put(245, (ctx.get(245) - ((ctx.get(235) / zA2) * zX2)));
ctx.put(238, ((-(ctx.get(235) / zA2)) * zMs));
ctx.put(239, ((-(ctx.get(235) / zA2)) * zRe));
ctx.put(246, 0.0);
}

// js/xblsys.js:847 (blsys)
void fn_blsys(Context ctxIn) {
Context ctx = ctxIn;
((void)0);
double simi = ctx.get(183);
double tran = ctx.get(195);
double turb = ctx.get(199);
double wake = ctx.get(227);
V1 COM1 = ctx.arr<1>(3);
V1 COM2 = ctx.arr<1>(4);
V1 VS1F = ctx.arr<1>(43);
V1 VS2F = ctx.arr<1>(45);
double vsStride = 6.0;
V1 VSREZ = ctx.arr<1>(48);
V1 VSM = ctx.arr<1>(46);
double u1Uei = ctx.get(202);
double u2Uei = ctx.get(205);
double u1Ms = ctx.get(201);
double u2Ms = ctx.get(204);
if (truth(wake)) {
fn_blvar(3.0, ctx);
fn_blmid(3.0, ctx);
} else if ((truth(turb) || truth(tran))) {
fn_blvar(2.0, ctx);
fn_blmid(2.0, ctx);
} else {
fn_blvar(1.0, ctx);
fn_blmid(1.0, ctx);
}
if (truth(simi)) {
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(COM1.at(icom) = COM2.at(icom));
}
fn_syncComToVars(ctx, 1.0);
}
if (truth(tran)) {
fn_trdif(ctx);
} else if (truth(simi)) {
fn_bldif(0.0, ctx);
} else if (truth((!truth(turb)))) {
fn_bldif(1.0, ctx);
} else if (truth(wake)) {
fn_bldif(3.0, ctx);
} else {
fn_bldif(2.0, ctx);
}
if (truth(simi)) {
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
double base = (double(k) * vsStride);
for (u32 l = 1u; truth((double(l) <= 5.0)); advance_index(l)) {
double idx = (base + double(l));
(VS2F[idx] = (VS1F[idx] + VS2F[idx]));
(VS1F[idx] = 0.0);
}
}
}
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
double base = (double(k) * vsStride);
double resU1 = VS1F[(base + 4.0)];
double resU2 = VS2F[(base + 4.0)];
double resMs = VSM.at(k);
(VS1F[(base + 4.0)] = (resU1 * u1Uei));
(VS2F[(base + 4.0)] = (resU2 * u2Uei));
(VSM.at(k) = (((resU1 * u1Ms) + (resU2 * u2Ms)) + resMs));
}
}

// js/xblsys.js:931 (tesys)
void fn_tesys(double cte, double tte, double dte, Context ctxIn) {
Context ctx = ctxIn;
((void)0);
V1 VS1F = ctx.arr<1>(43);
V1 VS2F = ctx.arr<1>(45);
double vsStride = 6.0;
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
(ctx.arr<1>(48).at(k) = 0.0);
(ctx.arr<1>(46).at(k) = 0.0);
(ctx.arr<1>(47).at(k) = 0.0);
(ctx.arr<1>(49).at(k) = 0.0);
double base = (double(k) * vsStride);
for (u32 l = 1u; truth((double(l) <= 5.0)); advance_index(l)) {
double idx = (base + double(l));
(VS1F[idx] = 0.0);
(VS2F[idx] = 0.0);
}
}
fn_blvar(3.0, ctx);
(VS1F[((1.0 * vsStride) + 1.0)] = (-1.0));
(VS2F[((1.0 * vsStride) + 1.0)] = 1.0);
(ctx.arr<1>(48)[1.0] = (cte - ctx.get(181)));
(VS1F[((2.0 * vsStride) + 2.0)] = (-1.0));
(VS2F[((2.0 * vsStride) + 2.0)] = 1.0);
(ctx.arr<1>(48)[2.0] = (tte - ctx.get(189)));
(VS1F[((3.0 * vsStride) + 3.0)] = (-1.0));
(VS2F[((3.0 * vsStride) + 3.0)] = 1.0);
(ctx.arr<1>(48)[3.0] = ((dte - ctx.get(53)) - ctx.get(82)));
}

// js/xblsys.js:970 (blprv)
void fn_blprv(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctxIn) {
Context ctx = ctxIn;
fn_predictBlStation(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_syncVarsToCom(ctx, 2.0);
}

// js/xblsys.js:977 (blprvKin)
void fn_blprvKin(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctxIn) {
Context ctx = ctxIn;
fn_predictBlStation(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_blkin(ctx);
}

// js/xblsys.js:983 (predictBlStation)
void fn_predictBlStation(double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei, Context ctx) {
((void)0);
ctx.put(229, xsi);
ctx.put(2, ami);
ctx.put(181, cti);
ctx.put(189, thi);
ctx.put(53, (dsi - dswaki));
ctx.put(82, dswaki);
ctx.put(203, ((uei * (1.0 - ctx.get(190))) / (1.0 - (ctx.get(190) * host_pow((uei / ctx.get(153)), 2.0)))));
ctx.put(205, ((1.0 + (ctx.get(190) * ((((2.0 * ctx.get(203)) * uei) / host_pow(ctx.get(153), 2.0)) - 1.0))) / (1.0 - (ctx.get(190) * host_pow((uei / ctx.get(153)), 2.0)))));
ctx.put(204, ((((ctx.get(203) * host_pow((uei / ctx.get(153)), 2.0)) - uei) * ctx.get(191)) / (1.0 - (ctx.get(190) * host_pow((uei / ctx.get(153)), 2.0)))));
}

// js/xblsys.js:1001 (blkin)
void fn_blkin(Context ctxIn) {
Context ctx = ctxIn;
((void)0);
double u2 = ctx.get(203);
double gm1bl = ctx.get(90);
double hstinv = ctx.get(130);
double hstinvMs = ctx.get(131);
double rstbl = ctx.get(168);
double rstblMs = ctx.get(169);
ctx.put(142, (((u2 * u2) * hstinv) / (gm1bl * (1.0 - (((0.5 * u2) * u2) * hstinv)))));
double tr2 = (1.0 + ((0.5 * gm1bl) * ctx.get(142)));
ctx.put(144, (((2.0 * ctx.get(142)) * tr2) / u2));
ctx.put(143, ((((u2 * u2) * tr2) / (gm1bl * (1.0 - (((0.5 * u2) * u2) * hstinv)))) * hstinvMs));
ctx.put(158, (rstbl * host_pow(tr2, ((-1.0) / gm1bl))));
ctx.put(160, ((((-ctx.get(158)) / tr2) * 0.5) * ctx.get(144)));
ctx.put(159, (((((-ctx.get(158)) / tr2) * 0.5) * ctx.get(143)) + (rstblMs * host_pow(tr2, ((-1.0) / gm1bl)))));
ctx.put(94, (ctx.get(53) / ctx.get(189)));
ctx.put(95, (1.0 / ctx.get(189)));
ctx.put(96, ((-ctx.get(94)) / ctx.get(189)));
double herat = (1.0 - (((0.5 * u2) * u2) * hstinv));
double heU2 = ((-u2) * hstinv);
double heMs = ((((-0.5) * u2) * u2) * hstinvMs);
ctx.put(222, (((__builtin_sqrt(host_pow(herat, 3.0)) * (1.0 + ctx.get(132))) / (herat + ctx.get(132))) / ctx.get(165)));
double v2He = (ctx.get(222) * ((1.5 / herat) - (1.0 / (herat + ctx.get(132)))));
ctx.put(225, (v2He * heU2));
ctx.put(223, ((((-ctx.get(222)) / ctx.get(165)) * ctx.get(166)) + (v2He * heMs)));
ctx.put(224, (((-ctx.get(222)) / ctx.get(165)) * ctx.get(167)));
Record hkinOut = fn_hkin(ctx.get(94), ctx.get(142));
ctx.put(113, hkinOut.hk);
double hk2H2 = hkinOut.hkH;
double hk2M2 = hkinOut.hkMsq;
ctx.put(117, (hk2M2 * ctx.get(144)));
ctx.put(116, (hk2H2 * ctx.get(96)));
ctx.put(114, (hk2H2 * ctx.get(95)));
ctx.put(115, (hk2M2 * ctx.get(143)));
ctx.put(175, (((ctx.get(158) * ctx.get(203)) * ctx.get(189)) / ctx.get(222)));
ctx.put(179, (ctx.get(175) * (((1.0 / ctx.get(203)) + (ctx.get(160) / ctx.get(158))) - (ctx.get(225) / ctx.get(222)))));
ctx.put(178, (ctx.get(175) / ctx.get(189)));
ctx.put(176, (ctx.get(175) * ((ctx.get(159) / ctx.get(158)) - (ctx.get(223) / ctx.get(222)))));
ctx.put(177, (ctx.get(175) * ((-ctx.get(224)) / ctx.get(222))));
fn_syncVarsToCom(ctx, 2.0);
}

// js/xblsys.js:1059 (blvar)
void fn_blvar(double ityp, Context ctxIn) {
Context ctx = ctxIn;
((void)0);
double gbcon = ctx.get(88);
double ctcon = ctx.get(48);
double gccon = ctx.get(89);
double h2 = ctx.get(94);
double h2U2 = ctx.get(97);
double h2T2 = ctx.get(96);
double h2D2 = ctx.get(95);
double m2 = ctx.get(142);
double m2U2 = ctx.get(144);
double m2Ms = ctx.get(143);
double rt2 = ctx.get(175);
double hk2U2 = ctx.get(117);
double hk2T2 = ctx.get(116);
double hk2D2 = ctx.get(114);
double hk2Ms = ctx.get(115);
double rt2U2 = ctx.get(179);
double rt2T2 = ctx.get(178);
double rt2Ms = ctx.get(176);
double rt2Re = ctx.get(177);
if (truth((ityp == 3.0))) ctx.put(113, js_max(ctx.get(113), 1.00005));
if (truth((ityp != 3.0))) ctx.put(113, js_max(ctx.get(113), 1.05000));
double hk2 = ctx.get(113);
Record hcOut = fn_hct(hk2, m2);
ctx.put(103, hcOut.hc);
double hc2Hk2 = hcOut.hcHk;
double hc2M2 = hcOut.hcMsq;
ctx.put(107, ((hc2Hk2 * hk2U2) + (hc2M2 * m2U2)));
ctx.put(106, (hc2Hk2 * hk2T2));
ctx.put(104, (hc2Hk2 * hk2D2));
ctx.put(105, ((hc2Hk2 * hk2Ms) + (hc2M2 * m2Ms)));
Record hsOut = Record{};
if (truth((ityp == 1.0))) {
(hsOut = fn_hsl(hk2, rt2, m2));
} else {
(hsOut = fn_hst(hk2, rt2, m2));
}
ctx.put(124, hsOut.hs);
double hs2Hk2 = hsOut.hsHk;
double hs2Rt2 = hsOut.hsRt;
double hs2M2 = hsOut.hsMsq;
ctx.put(129, (((hs2Hk2 * hk2U2) + (hs2Rt2 * rt2U2)) + (hs2M2 * m2U2)));
ctx.put(128, ((hs2Hk2 * hk2T2) + (hs2Rt2 * rt2T2)));
ctx.put(125, (hs2Hk2 * hk2D2));
ctx.put(126, (((hs2Hk2 * hk2Ms) + (hs2Rt2 * rt2Ms)) + (hs2M2 * m2Ms)));
ctx.put(127, (hs2Rt2 * rt2Re));
ctx.put(212, ((0.5 * ctx.get(124)) * (1.0 - ((hk2 - 1.0) / (gbcon * h2)))));
double us2Hs2 = (0.5 * (1.0 - ((hk2 - 1.0) / (gbcon * h2))));
double us2Hk2 = ((0.5 * ctx.get(124)) * ((-1.0) / (gbcon * h2)));
double us2H2 = (((0.5 * ctx.get(124)) * (hk2 - 1.0)) / (gbcon * host_pow(h2, 2.0)));
ctx.put(217, ((us2Hs2 * ctx.get(129)) + (us2Hk2 * hk2U2)));
ctx.put(216, (((us2Hs2 * ctx.get(128)) + (us2Hk2 * hk2T2)) + (us2H2 * h2T2)));
ctx.put(213, (((us2Hs2 * ctx.get(125)) + (us2Hk2 * hk2D2)) + (us2H2 * h2D2)));
ctx.put(214, ((us2Hs2 * ctx.get(126)) + (us2Hk2 * hk2Ms)));
ctx.put(215, (us2Hs2 * ctx.get(127)));
if ((truth((ityp <= 2.0)) && truth((ctx.get(212) > 0.95)))) {
ctx.put(212, 0.98);
ctx.put(217, 0.0);
ctx.put(216, 0.0);
ctx.put(213, 0.0);
ctx.put(214, 0.0);
ctx.put(215, 0.0);
}
if ((truth((ityp == 3.0)) && truth((ctx.get(212) > 0.99995)))) {
ctx.put(212, 0.99995);
ctx.put(217, 0.0);
ctx.put(216, 0.0);
ctx.put(213, 0.0);
ctx.put(214, 0.0);
ctx.put(215, 0.0);
}
double gcc = 0.0;
double hkc = (hk2 - 1.0);
double hkcHk2 = 1.0;
double hkcRt2 = 0.0;
if (truth((ityp == 2.0))) {
(gcc = gccon);
(hkc = ((hk2 - 1.0) - (gcc / ctx.get(175))));
(hkcHk2 = 1.0);
(hkcRt2 = (gcc / host_pow(ctx.get(175), 2.0)));
if (truth((hkc < 0.01))) {
(hkc = 0.01);
(hkcHk2 = 0.0);
(hkcRt2 = 0.0);
}
}
double hkb = (hk2 - 1.0);
double usb = (1.0 - ctx.get(212));
ctx.put(42, __builtin_sqrt(((((ctcon * ctx.get(124)) * hkb) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 2.0)))));
double cq2Hs2 = (((((ctcon * hkb) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 2.0))) * 0.5) / ctx.get(42));
double cq2Us2 = (((((((ctcon * ctx.get(124)) * hkb) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 2.0))) / usb) * 0.5) / ctx.get(42));
double cq2Hk2 = (((((((ctcon * ctx.get(124)) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 2.0))) * 0.5) / ctx.get(42)) - (((((((ctcon * ctx.get(124)) * hkb) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 3.0))) * 2.0) * 0.5) / ctx.get(42))) + ((((((((ctcon * ctx.get(124)) * hkb) * hkc) / ((usb * h2) * host_pow(hk2, 2.0))) * 2.0) * 0.5) / ctx.get(42)) * hkcHk2));
double cq2Rt2 = ((((((((ctcon * ctx.get(124)) * hkb) * hkc) / ((usb * h2) * host_pow(hk2, 2.0))) * 2.0) * 0.5) / ctx.get(42)) * hkcRt2);
double cq2H2 = ((((((((-ctcon) * ctx.get(124)) * hkb) * host_pow(hkc, 2.0)) / ((usb * h2) * host_pow(hk2, 2.0))) / h2) * 0.5) / ctx.get(42));
ctx.put(47, (((cq2Hs2 * ctx.get(129)) + (cq2Us2 * ctx.get(217))) + (cq2Hk2 * hk2U2)));
ctx.put(46, (((cq2Hs2 * ctx.get(128)) + (cq2Us2 * ctx.get(216))) + (cq2Hk2 * hk2T2)));
ctx.put(43, (((cq2Hs2 * ctx.get(125)) + (cq2Us2 * ctx.get(213))) + (cq2Hk2 * hk2D2)));
ctx.put(44, (((cq2Hs2 * ctx.get(126)) + (cq2Us2 * ctx.get(214))) + (cq2Hk2 * hk2Ms)));
ctx.put(45, ((cq2Hs2 * ctx.get(127)) + (cq2Us2 * ctx.get(215))));
ctx.put(47, (ctx.get(47) + (cq2Rt2 * rt2U2)));
ctx.put(46, (ctx.get(46) + ((cq2H2 * h2T2) + (cq2Rt2 * rt2T2))));
ctx.put(43, (ctx.get(43) + (cq2H2 * h2D2)));
ctx.put(44, (ctx.get(44) + (cq2Rt2 * rt2Ms)));
ctx.put(45, (ctx.get(45) + (cq2Rt2 * rt2Re)));
Record cfOut = Record{};
if (truth((ityp == 3.0))) {
ctx.put(11, 0.0);
ctx.put(13, 0.0);
ctx.put(17, 0.0);
ctx.put(14, 0.0);
} else if (truth((ityp == 1.0))) {
(cfOut = fn_cfl(hk2, ctx.get(175), ctx.get(142)));
ctx.put(11, cfOut.cf);
ctx.put(13, cfOut.cfHk);
ctx.put(17, cfOut.cfRt);
ctx.put(14, cfOut.cfMsq);
} else {
(cfOut = fn_cft(hk2, ctx.get(175), ctx.get(142), (ctx.has(20) ? ctx.get(20) : 1.0)));
ctx.put(11, cfOut.cf);
ctx.put(13, cfOut.cfHk);
ctx.put(17, cfOut.cfRt);
ctx.put(14, cfOut.cfMsq);
Record cfL = fn_cfl(hk2, ctx.get(175), ctx.get(142));
if (truth((cfL.cf > ctx.get(11)))) {
ctx.put(11, cfL.cf);
ctx.put(13, cfL.cfHk);
ctx.put(17, cfL.cfRt);
ctx.put(14, cfL.cfMsq);
}
}
ctx.put(19, (((ctx.get(13) * hk2U2) + (ctx.get(17) * rt2U2)) + (ctx.get(14) * ctx.get(144))));
ctx.put(18, ((ctx.get(13) * hk2T2) + (ctx.get(17) * rt2T2)));
ctx.put(12, (ctx.get(13) * hk2D2));
ctx.put(15, (((ctx.get(13) * hk2Ms) + (ctx.get(17) * rt2Ms)) + (ctx.get(14) * ctx.get(143))));
ctx.put(16, (ctx.get(17) * rt2Re));
if (truth((ityp == 1.0))) {
Record diOut = fn_dil(hk2, ctx.get(175));
ctx.put(71, diOut.di);
double di2Hk2 = diOut.diHk;
double di2Rt2 = diOut.diRt;
ctx.put(77, ((di2Hk2 * hk2U2) + (di2Rt2 * rt2U2)));
ctx.put(76, ((di2Hk2 * hk2T2) + (di2Rt2 * rt2T2)));
ctx.put(72, (di2Hk2 * hk2D2));
ctx.put(75, 0.0);
ctx.put(73, ((di2Hk2 * hk2Ms) + (di2Rt2 * rt2Ms)));
ctx.put(74, (di2Rt2 * rt2Re));
} else if (truth((ityp == 2.0))) {
Record cf2t = fn_cft(hk2, ctx.get(175), ctx.get(142), (ctx.has(20) ? ctx.get(20) : 1.0));
double cf2tU2 = (((cf2t.cfHk * hk2U2) + (cf2t.cfRt * rt2U2)) + (cf2t.cfMsq * ctx.get(144)));
double cf2tT2 = ((cf2t.cfHk * hk2T2) + (cf2t.cfRt * rt2T2));
double cf2tD2 = (cf2t.cfHk * hk2D2);
double cf2tMs = (((cf2t.cfHk * hk2Ms) + (cf2t.cfRt * rt2Ms)) + (cf2t.cfMsq * ctx.get(143)));
double cf2tRe = (cf2t.cfRt * rt2Re);
ctx.put(71, ((((0.5 * cf2t.cf) * ctx.get(212)) * 2.0) / ctx.get(124)));
double di2Hs2 = (((-((0.5 * cf2t.cf) * ctx.get(212))) * 2.0) / host_pow(ctx.get(124), 2.0));
double di2Us2 = (((0.5 * cf2t.cf) * 2.0) / ctx.get(124));
double di2Cf2t = (((0.5 * ctx.get(212)) * 2.0) / ctx.get(124));
ctx.put(75, 0.0);
ctx.put(77, (((di2Hs2 * ctx.get(129)) + (di2Us2 * ctx.get(217))) + (di2Cf2t * cf2tU2)));
ctx.put(76, (((di2Hs2 * ctx.get(128)) + (di2Us2 * ctx.get(216))) + (di2Cf2t * cf2tT2)));
ctx.put(72, (((di2Hs2 * ctx.get(125)) + (di2Us2 * ctx.get(213))) + (di2Cf2t * cf2tD2)));
ctx.put(73, (((di2Hs2 * ctx.get(126)) + (di2Us2 * ctx.get(214))) + (di2Cf2t * cf2tMs)));
ctx.put(74, (((di2Hs2 * ctx.get(127)) + (di2Us2 * ctx.get(215))) + (di2Cf2t * cf2tRe)));
double grt = host_log(ctx.get(175));
double hmin = (1.0 + (2.1 / grt));
double hmRt2 = ((-(2.1 / host_pow(grt, 2.0))) / ctx.get(175));
double fl = ((ctx.get(113) - 1.0) / (hmin - 1.0));
double flHk2 = (1.0 / (hmin - 1.0));
double flRt2 = (((-fl) / (hmin - 1.0)) * hmRt2);
double tfl = host_tanh(fl);
double dfac = (0.5 + (0.5 * tfl));
double dfFl = (0.5 * (1.0 - host_pow(tfl, 2.0)));
double dfHk2 = (dfFl * flHk2);
double dfRt2 = (dfFl * flRt2);
ctx.put(75, (ctx.get(75) * dfac));
ctx.put(77, ((ctx.get(77) * dfac) + (ctx.get(71) * ((dfHk2 * hk2U2) + (dfRt2 * rt2U2)))));
ctx.put(76, ((ctx.get(76) * dfac) + (ctx.get(71) * ((dfHk2 * hk2T2) + (dfRt2 * rt2T2)))));
ctx.put(72, ((ctx.get(72) * dfac) + (ctx.get(71) * (dfHk2 * hk2D2))));
ctx.put(73, ((ctx.get(73) * dfac) + (ctx.get(71) * ((dfHk2 * hk2Ms) + (dfRt2 * rt2Ms)))));
ctx.put(74, ((ctx.get(74) * dfac) + (ctx.get(71) * (dfRt2 * rt2Re))));
ctx.put(71, (ctx.get(71) * dfac));
} else {
ctx.put(71, 0.0);
ctx.put(75, 0.0);
ctx.put(77, 0.0);
ctx.put(76, 0.0);
ctx.put(72, 0.0);
ctx.put(73, 0.0);
ctx.put(74, 0.0);
}
if (truth((ityp != 1.0))) {
double dd = (((host_pow(ctx.get(181), 2.0) * (0.995 - ctx.get(212))) * 2.0) / ctx.get(124));
double ddHs2 = ((((-host_pow(ctx.get(181), 2.0)) * (0.995 - ctx.get(212))) * 2.0) / host_pow(ctx.get(124), 2.0));
double ddUs2 = (((-host_pow(ctx.get(181), 2.0)) * 2.0) / ctx.get(124));
double ddS2 = ((((ctx.get(181) * 2.0) * (0.995 - ctx.get(212))) * 2.0) / ctx.get(124));
ctx.put(71, (ctx.get(71) + dd));
ctx.put(75, ddS2);
ctx.put(77, (ctx.get(77) + ((ddHs2 * ctx.get(129)) + (ddUs2 * ctx.get(217)))));
ctx.put(76, (ctx.get(76) + ((ddHs2 * ctx.get(128)) + (ddUs2 * ctx.get(216)))));
ctx.put(72, (ctx.get(72) + ((ddHs2 * ctx.get(125)) + (ddUs2 * ctx.get(213)))));
ctx.put(73, (ctx.get(73) + ((ddHs2 * ctx.get(126)) + (ddUs2 * ctx.get(214)))));
ctx.put(74, (ctx.get(74) + ((ddHs2 * ctx.get(127)) + (ddUs2 * ctx.get(215)))));
(dd = ((((0.15 * host_pow((0.995 - ctx.get(212)), 2.0)) / ctx.get(175)) * 2.0) / ctx.get(124)));
double ddUs2b = ((((((-0.15) * (0.995 - ctx.get(212))) * 2.0) / ctx.get(175)) * 2.0) / ctx.get(124));
double ddHs2b = ((-dd) / ctx.get(124));
double ddRt2 = ((-dd) / ctx.get(175));
ctx.put(71, (ctx.get(71) + dd));
ctx.put(77, (ctx.get(77) + (((ddHs2b * ctx.get(129)) + (ddUs2b * ctx.get(217))) + (ddRt2 * ctx.get(179)))));
ctx.put(76, (ctx.get(76) + (((ddHs2b * ctx.get(128)) + (ddUs2b * ctx.get(216))) + (ddRt2 * ctx.get(178)))));
ctx.put(72, (ctx.get(72) + ((ddHs2b * ctx.get(125)) + (ddUs2b * ctx.get(213)))));
ctx.put(73, (ctx.get(73) + (((ddHs2b * ctx.get(126)) + (ddUs2b * ctx.get(214))) + (ddRt2 * ctx.get(176)))));
ctx.put(74, (ctx.get(74) + (((ddHs2b * ctx.get(127)) + (ddUs2b * ctx.get(215))) + (ddRt2 * ctx.get(177)))));
}
if (truth((ityp == 2.0))) {
Record di2l = fn_dil(hk2, ctx.get(175));
if (truth((di2l.di > ctx.get(71)))) {
ctx.put(71, di2l.di);
ctx.put(75, 0.0);
ctx.put(77, ((di2l.diHk * hk2U2) + (di2l.diRt * rt2U2)));
ctx.put(76, ((di2l.diHk * hk2T2) + (di2l.diRt * rt2T2)));
ctx.put(72, (di2l.diHk * hk2D2));
ctx.put(73, ((di2l.diHk * hk2Ms) + (di2l.diRt * rt2Ms)));
ctx.put(74, (di2l.diRt * rt2Re));
}
}
if (truth((ityp == 3.0))) {
Record di2l = fn_dilw(ctx.get(113), ctx.get(175));
if (truth((di2l.di > ctx.get(71)))) {
ctx.put(71, di2l.di);
ctx.put(75, 0.0);
ctx.put(77, ((di2l.diHk * ctx.get(117)) + (di2l.diRt * ctx.get(179))));
ctx.put(76, ((di2l.diHk * ctx.get(116)) + (di2l.diRt * ctx.get(178))));
ctx.put(72, (di2l.diHk * ctx.get(114)));
ctx.put(73, ((di2l.diHk * ctx.get(115)) + (di2l.diRt * ctx.get(176))));
ctx.put(74, (di2l.diRt * ctx.get(177)));
}
}
if (truth((ityp == 3.0))) {
ctx.put(71, (ctx.get(71) * 2.0));
ctx.put(75, (ctx.get(75) * 2.0));
ctx.put(77, (ctx.get(77) * 2.0));
ctx.put(76, (ctx.get(76) * 2.0));
ctx.put(72, (ctx.get(72) * 2.0));
ctx.put(73, (ctx.get(73) * 2.0));
ctx.put(74, (ctx.get(74) * 2.0));
}
ctx.put(59, (((3.15 + (1.72 / (ctx.get(113) - 1.0))) * ctx.get(189)) + ctx.get(53)));
double de2Hk2 = (((-1.72) / host_pow((ctx.get(113) - 1.0), 2.0)) * ctx.get(189));
ctx.put(63, (de2Hk2 * ctx.get(117)));
ctx.put(62, ((de2Hk2 * ctx.get(116)) + (3.15 + (1.72 / (ctx.get(113) - 1.0)))));
ctx.put(60, ((de2Hk2 * ctx.get(114)) + 1.0));
ctx.put(61, (de2Hk2 * ctx.get(115)));
double hdmax = 12.0;
if (truth((ctx.get(59) > (hdmax * ctx.get(189))))) {
ctx.put(59, (hdmax * ctx.get(189)));
ctx.put(63, 0.0);
ctx.put(62, hdmax);
ctx.put(60, 0.0);
ctx.put(61, 0.0);
}
fn_syncVarsToCom(ctx, 2.0);
}

// js/xblsys.js:1362 (blmid)
void fn_blmid(double ityp, Context ctxIn) {
Context ctx = ctxIn;
((void)0);
if (truth(ctx.get(183))) {
ctx.put(108, ctx.get(113));
ctx.put(111, ctx.get(116));
ctx.put(109, ctx.get(114));
ctx.put(112, ctx.get(117));
ctx.put(110, ctx.get(115));
ctx.put(170, ctx.get(175));
ctx.put(173, ctx.get(178));
ctx.put(174, ctx.get(179));
ctx.put(171, ctx.get(176));
ctx.put(172, ctx.get(177));
ctx.put(139, ctx.get(142));
ctx.put(141, ctx.get(144));
ctx.put(140, ctx.get(143));
}
double hka = (0.5 * (ctx.get(108) + ctx.get(113)));
double rta = (0.5 * (ctx.get(170) + ctx.get(175)));
double ma = (0.5 * (ctx.get(139) + ctx.get(142)));
if (truth((ityp == 3.0))) {
ctx.put(21, 0.0);
ctx.put(24, 0.0);
ctx.put(28, 0.0);
ctx.put(25, 0.0);
ctx.put(26, 0.0);
} else if (truth((ityp == 1.0))) {
Record cfOut = fn_cfl(hka, rta, ma);
ctx.put(21, cfOut.cf);
ctx.put(24, cfOut.cfHk);
ctx.put(28, cfOut.cfRt);
ctx.put(25, cfOut.cfMsq);
} else {
Record cfOut = fn_cft(hka, rta, ma, (ctx.has(20) ? ctx.get(20) : 1.0));
ctx.put(21, cfOut.cf);
ctx.put(24, cfOut.cfHk);
ctx.put(28, cfOut.cfRt);
ctx.put(25, cfOut.cfMsq);
Record cfml = fn_cfl(hka, rta, ma);
if (truth((cfml.cf > ctx.get(21)))) {
ctx.put(21, cfml.cf);
ctx.put(24, cfml.cfHk);
ctx.put(28, cfml.cfRt);
ctx.put(25, cfml.cfMsq);
}
}
ctx.put(31, (0.5 * (((ctx.get(24) * ctx.get(112)) + (ctx.get(25) * ctx.get(141))) + (ctx.get(28) * ctx.get(174)))));
ctx.put(29, (0.5 * ((ctx.get(24) * ctx.get(111)) + (ctx.get(28) * ctx.get(173)))));
ctx.put(22, (0.5 * (ctx.get(24) * ctx.get(109))));
ctx.put(32, (0.5 * (((ctx.get(24) * ctx.get(117)) + (ctx.get(25) * ctx.get(144))) + (ctx.get(28) * ctx.get(179)))));
ctx.put(30, (0.5 * ((ctx.get(24) * ctx.get(116)) + (ctx.get(28) * ctx.get(178)))));
ctx.put(23, (0.5 * (ctx.get(24) * ctx.get(114))));
ctx.put(26, (0.5 * ((((((ctx.get(24) * ctx.get(110)) + (ctx.get(25) * ctx.get(140))) + (ctx.get(28) * ctx.get(171))) + (ctx.get(24) * ctx.get(115))) + (ctx.get(25) * ctx.get(143))) + (ctx.get(28) * ctx.get(176)))));
ctx.put(27, (0.5 * ((ctx.get(28) * ctx.get(172)) + (ctx.get(28) * ctx.get(177)))));
}

// js/xblsys.js:1429 (trdif)
void fn_trdif(Context ctxIn) {
Context ctx = ctxIn;
((void)0);
V1 VS1F = ctx.arr<1>(43);
V1 VS2F = ctx.arr<1>(45);
double vsStride = 6.0;
V2 bl1 = ctx.arr<2>(70);
V2 bl2 = ctx.arr<2>(71);
V2 bt1 = ctx.arr<2>(76);
V2 bt2 = ctx.arr<2>(77);
V1 blrez = ctx.arr<1>(74);
V1 blm = ctx.arr<1>(72);
V1 blr = ctx.arr<1>(73);
V1 blx = ctx.arr<1>(75);
V1 btrez = ctx.arr<1>(80);
V1 btm = ctx.arr<1>(78);
V1 btr = ctx.arr<1>(79);
V1 btx = ctx.arr<1>(81);
if (truth((!bool(bl1)))) {
(bl1 = make_matrix(4.0 + 1, 5.0 + 1));
(bl2 = make_matrix(4.0 + 1, 5.0 + 1));
(bt1 = make_matrix(4.0 + 1, 5.0 + 1));
(bt2 = make_matrix(4.0 + 1, 5.0 + 1));
(blrez = make_vec<1>(5.0));
(blm = make_vec<1>(5.0));
(blr = make_vec<1>(5.0));
(blx = make_vec<1>(5.0));
(btrez = make_vec<1>(5.0));
(btm = make_vec<1>(5.0));
(btr = make_vec<1>(5.0));
(btx = make_vec<1>(5.0));
(ctx.arr<2>(70) = bl1);
(ctx.arr<2>(71) = bl2);
(ctx.arr<2>(76) = bt1);
(ctx.arr<2>(77) = bt2);
(ctx.arr<1>(74) = blrez);
(ctx.arr<1>(72) = blm);
(ctx.arr<1>(73) = blr);
(ctx.arr<1>(75) = blx);
(ctx.arr<1>(80) = btrez);
(ctx.arr<1>(78) = btm);
(ctx.arr<1>(79) = btr);
(ctx.arr<1>(81) = btx);
} else {
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
bl1.at(k).fill(0.0);
bl2.at(k).fill(0.0);
bt1.at(k).fill(0.0);
bt2.at(k).fill(0.0);
}
blrez.fill(0.0);
blm.fill(0.0);
blr.fill(0.0);
blx.fill(0.0);
btrez.fill(0.0);
btm.fill(0.0);
btr.fill(0.0);
btx.fill(0.0);
}
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(ctx.arr<1>(1).at(icom) = ctx.arr<1>(3).at(icom));
(ctx.arr<1>(2).at(icom) = ctx.arr<1>(4).at(icom));
}
double wf2 = ((ctx.get(232) - ctx.get(228)) / (ctx.get(229) - ctx.get(228)));
double wf2Xt = (1.0 / (ctx.get(229) - ctx.get(228)));
double wf2A1 = (wf2Xt * ctx.get(234));
double wf2X1 = ((wf2Xt * ctx.get(244)) + ((wf2 - 1.0) / (ctx.get(229) - ctx.get(228))));
double wf2X2 = ((wf2Xt * ctx.get(245)) - (wf2 / (ctx.get(229) - ctx.get(228))));
double wf2T1 = (wf2Xt * ctx.get(240));
double wf2T2 = (wf2Xt * ctx.get(241));
double wf2D1 = (wf2Xt * ctx.get(236));
double wf2D2 = (wf2Xt * ctx.get(237));
double wf2U1 = (wf2Xt * ctx.get(242));
double wf2U2 = (wf2Xt * ctx.get(243));
double wf2Ms = (wf2Xt * ctx.get(238));
double wf2Re = (wf2Xt * ctx.get(239));
double wf2Xf = (wf2Xt * ctx.get(246));
double wf1 = (1.0 - wf2);
double wf1A1 = (-wf2A1);
double wf1X1 = (-wf2X1);
double wf1X2 = (-wf2X2);
double wf1T1 = (-wf2T1);
double wf1T2 = (-wf2T2);
double wf1D1 = (-wf2D1);
double wf1D2 = (-wf2D2);
double wf1U1 = (-wf2U1);
double wf1U2 = (-wf2U2);
double wf1Ms = (-wf2Ms);
double wf1Re = (-wf2Re);
double wf1Xf = (-wf2Xf);
double tt = ((ctx.get(188) * wf1) + (ctx.get(189) * wf2));
double ttA1 = ((ctx.get(188) * wf1A1) + (ctx.get(189) * wf2A1));
double ttX1 = ((ctx.get(188) * wf1X1) + (ctx.get(189) * wf2X1));
double ttX2 = ((ctx.get(188) * wf1X2) + (ctx.get(189) * wf2X2));
double ttT1 = (((ctx.get(188) * wf1T1) + (ctx.get(189) * wf2T1)) + wf1);
double ttT2 = (((ctx.get(188) * wf1T2) + (ctx.get(189) * wf2T2)) + wf2);
double ttD1 = ((ctx.get(188) * wf1D1) + (ctx.get(189) * wf2D1));
double ttD2 = ((ctx.get(188) * wf1D2) + (ctx.get(189) * wf2D2));
double ttU1 = ((ctx.get(188) * wf1U1) + (ctx.get(189) * wf2U1));
double ttU2 = ((ctx.get(188) * wf1U2) + (ctx.get(189) * wf2U2));
double ttMs = ((ctx.get(188) * wf1Ms) + (ctx.get(189) * wf2Ms));
double ttRe = ((ctx.get(188) * wf1Re) + (ctx.get(189) * wf2Re));
double ttXf = ((ctx.get(188) * wf1Xf) + (ctx.get(189) * wf2Xf));
double dt = ((ctx.get(52) * wf1) + (ctx.get(53) * wf2));
double dtA1 = ((ctx.get(52) * wf1A1) + (ctx.get(53) * wf2A1));
double dtX1 = ((ctx.get(52) * wf1X1) + (ctx.get(53) * wf2X1));
double dtX2 = ((ctx.get(52) * wf1X2) + (ctx.get(53) * wf2X2));
double dtT1 = ((ctx.get(52) * wf1T1) + (ctx.get(53) * wf2T1));
double dtT2 = ((ctx.get(52) * wf1T2) + (ctx.get(53) * wf2T2));
double dtD1 = (((ctx.get(52) * wf1D1) + (ctx.get(53) * wf2D1)) + wf1);
double dtD2 = (((ctx.get(52) * wf1D2) + (ctx.get(53) * wf2D2)) + wf2);
double dtU1 = ((ctx.get(52) * wf1U1) + (ctx.get(53) * wf2U1));
double dtU2 = ((ctx.get(52) * wf1U2) + (ctx.get(53) * wf2U2));
double dtMs = ((ctx.get(52) * wf1Ms) + (ctx.get(53) * wf2Ms));
double dtRe = ((ctx.get(52) * wf1Re) + (ctx.get(53) * wf2Re));
double dtXf = ((ctx.get(52) * wf1Xf) + (ctx.get(53) * wf2Xf));
double ut = ((ctx.get(200) * wf1) + (ctx.get(203) * wf2));
double utA1 = ((ctx.get(200) * wf1A1) + (ctx.get(203) * wf2A1));
double utX1 = ((ctx.get(200) * wf1X1) + (ctx.get(203) * wf2X1));
double utX2 = ((ctx.get(200) * wf1X2) + (ctx.get(203) * wf2X2));
double utT1 = ((ctx.get(200) * wf1T1) + (ctx.get(203) * wf2T1));
double utT2 = ((ctx.get(200) * wf1T2) + (ctx.get(203) * wf2T2));
double utD1 = ((ctx.get(200) * wf1D1) + (ctx.get(203) * wf2D1));
double utD2 = ((ctx.get(200) * wf1D2) + (ctx.get(203) * wf2D2));
double utU1 = (((ctx.get(200) * wf1U1) + (ctx.get(203) * wf2U1)) + wf1);
double utU2 = (((ctx.get(200) * wf1U2) + (ctx.get(203) * wf2U2)) + wf2);
double utMs = ((ctx.get(200) * wf1Ms) + (ctx.get(203) * wf2Ms));
double utRe = ((ctx.get(200) * wf1Re) + (ctx.get(203) * wf2Re));
double utXf = ((ctx.get(200) * wf1Xf) + (ctx.get(203) * wf2Xf));
ctx.put(229, ctx.get(232));
ctx.put(189, tt);
ctx.put(53, dt);
ctx.put(203, ut);
ctx.put(2, ctx.get(0));
ctx.put(181, 0.0);
fn_blkin(ctx);
fn_blvar(1.0, ctx);
fn_blmid(1.0, ctx);
fn_bldif(1.0, ctx);
for (u32 k = 2u; truth((double(k) <= 3.0)); advance_index(k)) {
double base = (double(k) * vsStride);
(blrez.at(k) = ctx.arr<1>(48).at(k));
(blm.at(k) = ((((ctx.arr<1>(46).at(k) + (VS2F[(base + 2.0)] * ttMs)) + (VS2F[(base + 3.0)] * dtMs)) + (VS2F[(base + 4.0)] * utMs)) + (VS2F[(base + 5.0)] * ctx.get(238))));
(blr.at(k) = ((((ctx.arr<1>(47).at(k) + (VS2F[(base + 2.0)] * ttRe)) + (VS2F[(base + 3.0)] * dtRe)) + (VS2F[(base + 4.0)] * utRe)) + (VS2F[(base + 5.0)] * ctx.get(239))));
(blx.at(k) = ((((ctx.arr<1>(49).at(k) + (VS2F[(base + 2.0)] * ttXf)) + (VS2F[(base + 3.0)] * dtXf)) + (VS2F[(base + 4.0)] * utXf)) + (VS2F[(base + 5.0)] * ctx.get(246))));
(bl1.at(k)[1.0] = ((((VS1F[(base + 1.0)] + (VS2F[(base + 2.0)] * ttA1)) + (VS2F[(base + 3.0)] * dtA1)) + (VS2F[(base + 4.0)] * utA1)) + (VS2F[(base + 5.0)] * ctx.get(234))));
(bl1.at(k)[2.0] = ((((VS1F[(base + 2.0)] + (VS2F[(base + 2.0)] * ttT1)) + (VS2F[(base + 3.0)] * dtT1)) + (VS2F[(base + 4.0)] * utT1)) + (VS2F[(base + 5.0)] * ctx.get(240))));
(bl1.at(k)[3.0] = ((((VS1F[(base + 3.0)] + (VS2F[(base + 2.0)] * ttD1)) + (VS2F[(base + 3.0)] * dtD1)) + (VS2F[(base + 4.0)] * utD1)) + (VS2F[(base + 5.0)] * ctx.get(236))));
(bl1.at(k)[4.0] = ((((VS1F[(base + 4.0)] + (VS2F[(base + 2.0)] * ttU1)) + (VS2F[(base + 3.0)] * dtU1)) + (VS2F[(base + 4.0)] * utU1)) + (VS2F[(base + 5.0)] * ctx.get(242))));
(bl1.at(k)[5.0] = ((((VS1F[(base + 5.0)] + (VS2F[(base + 2.0)] * ttX1)) + (VS2F[(base + 3.0)] * dtX1)) + (VS2F[(base + 4.0)] * utX1)) + (VS2F[(base + 5.0)] * ctx.get(244))));
(bl2.at(k)[1.0] = 0.0);
(bl2.at(k)[2.0] = ((((VS2F[(base + 2.0)] * ttT2) + (VS2F[(base + 3.0)] * dtT2)) + (VS2F[(base + 4.0)] * utT2)) + (VS2F[(base + 5.0)] * ctx.get(241))));
(bl2.at(k)[3.0] = ((((VS2F[(base + 2.0)] * ttD2) + (VS2F[(base + 3.0)] * dtD2)) + (VS2F[(base + 4.0)] * utD2)) + (VS2F[(base + 5.0)] * ctx.get(237))));
(bl2.at(k)[4.0] = ((((VS2F[(base + 2.0)] * ttU2) + (VS2F[(base + 3.0)] * dtU2)) + (VS2F[(base + 4.0)] * utU2)) + (VS2F[(base + 5.0)] * ctx.get(243))));
(bl2.at(k)[5.0] = ((((VS2F[(base + 2.0)] * ttX2) + (VS2F[(base + 3.0)] * dtX2)) + (VS2F[(base + 4.0)] * utX2)) + (VS2F[(base + 5.0)] * ctx.get(245))));
}
fn_blvar(2.0, ctx);
double ctr = (ctx.get(51) * host_exp(((-ctx.get(50)) / (ctx.get(113) - 1.0))));
double ctrHk2 = ((ctr * ctx.get(50)) / host_pow((ctx.get(113) - 1.0), 2.0));
double st = (ctr * ctx.get(42));
double stTt = ((ctr * ctx.get(46)) + ((ctx.get(42) * ctrHk2) * ctx.get(116)));
double stDt = ((ctr * ctx.get(43)) + ((ctx.get(42) * ctrHk2) * ctx.get(114)));
double stUt = ((ctr * ctx.get(47)) + ((ctx.get(42) * ctrHk2) * ctx.get(117)));
double stMs = ((ctr * ctx.get(44)) + ((ctx.get(42) * ctrHk2) * ctx.get(115)));
double stRe = (ctr * ctx.get(45));
double stA1 = (((stTt * ttA1) + (stDt * dtA1)) + (stUt * utA1));
double stX1 = (((stTt * ttX1) + (stDt * dtX1)) + (stUt * utX1));
double stX2 = (((stTt * ttX2) + (stDt * dtX2)) + (stUt * utX2));
double stT1 = (((stTt * ttT1) + (stDt * dtT1)) + (stUt * utT1));
double stT2 = (((stTt * ttT2) + (stDt * dtT2)) + (stUt * utT2));
double stD1 = (((stTt * ttD1) + (stDt * dtD1)) + (stUt * utD1));
double stD2 = (((stTt * ttD2) + (stDt * dtD2)) + (stUt * utD2));
double stU1 = (((stTt * ttU1) + (stDt * dtU1)) + (stUt * utU1));
double stU2 = (((stTt * ttU2) + (stDt * dtU2)) + (stUt * utU2));
double stMs2 = ((((stTt * ttMs) + (stDt * dtMs)) + (stUt * utMs)) + stMs);
double stRe2 = ((((stTt * ttRe) + (stDt * dtRe)) + (stUt * utRe)) + stRe);
double stXf = (((stTt * ttXf) + (stDt * dtXf)) + (stUt * utXf));
ctx.put(2, 0.0);
ctx.put(181, st);
fn_blvar(2.0, ctx);
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(ctx.arr<1>(3).at(icom) = ctx.arr<1>(4).at(icom));
(ctx.arr<1>(4).at(icom) = ctx.arr<1>(2).at(icom));
}
fn_syncComToVars(ctx, 1.0);
fn_syncComToVars(ctx, 2.0);
fn_blmid(2.0, ctx);
fn_bldif(2.0, ctx);
for (u32 k = 1u; truth((double(k) <= 3.0)); advance_index(k)) {
double base = (double(k) * vsStride);
(btrez.at(k) = ctx.arr<1>(48).at(k));
(btm.at(k) = (((((ctx.arr<1>(46).at(k) + (VS1F[(base + 1.0)] * stMs2)) + (VS1F[(base + 2.0)] * ttMs)) + (VS1F[(base + 3.0)] * dtMs)) + (VS1F[(base + 4.0)] * utMs)) + (VS1F[(base + 5.0)] * ctx.get(238))));
(btr.at(k) = (((((ctx.arr<1>(47).at(k) + (VS1F[(base + 1.0)] * stRe2)) + (VS1F[(base + 2.0)] * ttRe)) + (VS1F[(base + 3.0)] * dtRe)) + (VS1F[(base + 4.0)] * utRe)) + (VS1F[(base + 5.0)] * ctx.get(239))));
(btx.at(k) = (((((ctx.arr<1>(49).at(k) + (VS1F[(base + 1.0)] * stXf)) + (VS1F[(base + 2.0)] * ttXf)) + (VS1F[(base + 3.0)] * dtXf)) + (VS1F[(base + 4.0)] * utXf)) + (VS1F[(base + 5.0)] * ctx.get(246))));
(bt1.at(k)[1.0] = (((((VS1F[(base + 1.0)] * stA1) + (VS1F[(base + 2.0)] * ttA1)) + (VS1F[(base + 3.0)] * dtA1)) + (VS1F[(base + 4.0)] * utA1)) + (VS1F[(base + 5.0)] * ctx.get(234))));
(bt1.at(k)[2.0] = (((((VS1F[(base + 1.0)] * stT1) + (VS1F[(base + 2.0)] * ttT1)) + (VS1F[(base + 3.0)] * dtT1)) + (VS1F[(base + 4.0)] * utT1)) + (VS1F[(base + 5.0)] * ctx.get(240))));
(bt1.at(k)[3.0] = (((((VS1F[(base + 1.0)] * stD1) + (VS1F[(base + 2.0)] * ttD1)) + (VS1F[(base + 3.0)] * dtD1)) + (VS1F[(base + 4.0)] * utD1)) + (VS1F[(base + 5.0)] * ctx.get(236))));
(bt1.at(k)[4.0] = (((((VS1F[(base + 1.0)] * stU1) + (VS1F[(base + 2.0)] * ttU1)) + (VS1F[(base + 3.0)] * dtU1)) + (VS1F[(base + 4.0)] * utU1)) + (VS1F[(base + 5.0)] * ctx.get(242))));
(bt1.at(k)[5.0] = (((((VS1F[(base + 1.0)] * stX1) + (VS1F[(base + 2.0)] * ttX1)) + (VS1F[(base + 3.0)] * dtX1)) + (VS1F[(base + 4.0)] * utX1)) + (VS1F[(base + 5.0)] * ctx.get(244))));
(bt2.at(k)[1.0] = VS2F[(base + 1.0)]);
(bt2.at(k)[2.0] = (((((VS2F[(base + 2.0)] + (VS1F[(base + 1.0)] * stT2)) + (VS1F[(base + 2.0)] * ttT2)) + (VS1F[(base + 3.0)] * dtT2)) + (VS1F[(base + 4.0)] * utT2)) + (VS1F[(base + 5.0)] * ctx.get(241))));
(bt2.at(k)[3.0] = (((((VS2F[(base + 3.0)] + (VS1F[(base + 1.0)] * stD2)) + (VS1F[(base + 2.0)] * ttD2)) + (VS1F[(base + 3.0)] * dtD2)) + (VS1F[(base + 4.0)] * utD2)) + (VS1F[(base + 5.0)] * ctx.get(237))));
(bt2.at(k)[4.0] = (((((VS2F[(base + 4.0)] + (VS1F[(base + 1.0)] * stU2)) + (VS1F[(base + 2.0)] * ttU2)) + (VS1F[(base + 3.0)] * dtU2)) + (VS1F[(base + 4.0)] * utU2)) + (VS1F[(base + 5.0)] * ctx.get(243))));
(bt2.at(k)[5.0] = (((((VS2F[(base + 5.0)] + (VS1F[(base + 1.0)] * stX2)) + (VS1F[(base + 2.0)] * ttX2)) + (VS1F[(base + 3.0)] * dtX2)) + (VS1F[(base + 4.0)] * utX2)) + (VS1F[(base + 5.0)] * ctx.get(245))));
}
(ctx.arr<1>(48)[1.0] = btrez[1.0]);
(ctx.arr<1>(48)[2.0] = (blrez[2.0] + btrez[2.0]));
(ctx.arr<1>(48)[3.0] = (blrez[3.0] + btrez[3.0]));
(ctx.arr<1>(46)[1.0] = btm[1.0]);
(ctx.arr<1>(46)[2.0] = (blm[2.0] + btm[2.0]));
(ctx.arr<1>(46)[3.0] = (blm[3.0] + btm[3.0]));
(ctx.arr<1>(47)[1.0] = btr[1.0]);
(ctx.arr<1>(47)[2.0] = (blr[2.0] + btr[2.0]));
(ctx.arr<1>(47)[3.0] = (blr[3.0] + btr[3.0]));
(ctx.arr<1>(49)[1.0] = btx[1.0]);
(ctx.arr<1>(49)[2.0] = (blx[2.0] + btx[2.0]));
(ctx.arr<1>(49)[3.0] = (blx[3.0] + btx[3.0]));
for (u32 l = 1u; truth((double(l) <= 5.0)); advance_index(l)) {
(VS1F[((1.0 * vsStride) + double(l))] = bt1[1.0].at(l));
(VS2F[((1.0 * vsStride) + double(l))] = bt2[1.0].at(l));
(VS1F[((2.0 * vsStride) + double(l))] = (bl1[2.0].at(l) + bt1[2.0].at(l)));
(VS2F[((2.0 * vsStride) + double(l))] = (bl2[2.0].at(l) + bt2[2.0].at(l)));
(VS1F[((3.0 * vsStride) + double(l))] = (bl1[3.0].at(l) + bt1[3.0].at(l)));
(VS2F[((3.0 * vsStride) + double(l))] = (bl2[3.0].at(l) + bt2[3.0].at(l)));
}
for (u32 icom = 1u; truth((double(icom) <= 73.0)); advance_index(icom)) {
(ctx.arr<1>(3).at(icom) = ctx.arr<1>(1).at(icom));
}
fn_syncComToVars(ctx, 1.0);
}

// js/xblsys.js:1789 (bldif)
void fn_bldif(double ityp, Context ctxIn) {
Context ctx = ctxIn;
((void)0);
V1 VS1F = ctx.arr<1>(43);
V1 VS2F = ctx.arr<1>(45);
double vsStride = 6.0;
V1 VSREZ = ctx.arr<1>(48);
V1 VSM = ctx.arr<1>(46);
V1 VSR = ctx.arr<1>(47);
V1 VSX = ctx.arr<1>(49);
double hk1 = ctx.get(108);
double hk2 = ctx.get(113);
double s1 = ctx.get(180);
double s2 = ctx.get(181);
double cq1 = ctx.get(36);
double cq2 = ctx.get(42);
double cf1 = ctx.get(5);
double cf2 = ctx.get(11);
double us1 = ctx.get(206);
double us2 = ctx.get(212);
double rt1 = ctx.get(170);
double rt2 = ctx.get(175);
double de1 = ctx.get(54);
double de2 = ctx.get(59);
double d1 = ctx.get(52);
double d2 = ctx.get(53);
double dlcon = ctx.get(78);
double gacon = ctx.get(84);
double gbcon = ctx.get(88);
double gccon = ctx.get(89);
double hk1U1 = ctx.get(112);
double hk1T1 = ctx.get(111);
double hk1D1 = ctx.get(109);
double hk1Ms = ctx.get(110);
double hk2U2 = ctx.get(117);
double hk2T2 = ctx.get(116);
double hk2D2 = ctx.get(114);
double hk2Ms = ctx.get(115);
double rt1T1 = ctx.get(173);
double rt1U1 = ctx.get(174);
double rt1Ms = ctx.get(171);
double rt1Re = ctx.get(172);
double rt2T2 = ctx.get(178);
double rt2U2 = ctx.get(179);
double rt2Ms = ctx.get(176);
double rt2Re = ctx.get(177);
double x1 = ctx.get(228);
double x2 = ctx.get(229);
double t1 = ctx.get(188);
double t2 = ctx.get(189);
double u1 = ctx.get(200);
double u2 = ctx.get(203);
double xlog = nan_value();
double ulog = nan_value();
double tlog = nan_value();
double hlog = nan_value();
double ddlog = nan_value();
if (truth((ityp == 0.0))) {
(xlog = 1.0);
(ulog = ctx.get(4));
(tlog = (0.5 * (1.0 - ctx.get(4))));
(hlog = 0.0);
(ddlog = 0.0);
} else {
(xlog = host_log((x2 / x1)));
(ulog = host_log((u2 / u1)));
(tlog = host_log((t2 / t1)));
(hlog = host_log((ctx.get(124) / ctx.get(118))));
(ddlog = 1.0);
}
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
(VSREZ.at(k) = 0.0);
(VSM.at(k) = 0.0);
(VSR.at(k) = 0.0);
(VSX.at(k) = 0.0);
double base = (double(k) * vsStride);
for (u32 l = 1u; truth((double(l) <= 5.0)); advance_index(l)) {
double idx = (base + double(l));
(VS1F[idx] = 0.0);
(VS2F[idx] = 0.0);
}
}
double hupwt = 1.0;
double hdcon = ((5.0 * hupwt) / host_pow(hk2, 2.0));
double hdHk1 = 0.0;
double hdHk2 = (((-hdcon) * 2.0) / hk2);
if (truth((ityp == 3.0))) {
(hdcon = (hupwt / host_pow(hk2, 2.0)));
(hdHk1 = 0.0);
(hdHk2 = (((-hdcon) * 2.0) / hk2));
}
double arg = __builtin_fabs(((hk2 - 1.0) / (hk1 - 1.0)));
double hl = host_log(arg);
double hlHk1 = ((-1.0) / (hk1 - 1.0));
double hlHk2 = (1.0 / (hk2 - 1.0));
double hlsq = js_min(host_pow(hl, 2.0), 15.0);
double ehh = host_exp(((-hlsq) * hdcon));
double upw = (1.0 - (0.5 * ehh));
double upwHl = ((ehh * hl) * hdcon);
double upwHd = ((0.5 * ehh) * hlsq);
double upwHk1 = ((upwHl * hlHk1) + (upwHd * hdHk1));
double upwHk2 = ((upwHl * hlHk2) + (upwHd * hdHk2));
double upwU1 = (upwHk1 * hk1U1);
double upwT1 = (upwHk1 * hk1T1);
double upwD1 = (upwHk1 * hk1D1);
double upwU2 = (upwHk2 * hk2U2);
double upwT2 = (upwHk2 * hk2T2);
double upwD2 = (upwHk2 * hk2D2);
double upwMs = ((upwHk1 * hk1Ms) + (upwHk2 * hk2Ms));
if (truth((ityp == 0.0))) {
(VS2F[((1.0 * vsStride) + 1.0)] = 1.0);
(VSR[1.0] = 0.0);
(VSREZ[1.0] = (-ctx.get(2)));
} else if (truth((ityp == 1.0))) {
Record axOut = fn_axset(hk1, t1, rt1, ctx.get(1), hk2, t2, rt2, ctx.get(2), ctx.get(0), ctx.get(134));
double ax = axOut.ax;
double rezc = ((ctx.get(2) - ctx.get(1)) - (ax * (x2 - x1)));
double zAx = (-(x2 - x1));
(VS1F[((1.0 * vsStride) + 1.0)] = ((zAx * axOut.axA1) - 1.0));
(VS1F[((1.0 * vsStride) + 2.0)] = (zAx * (((axOut.axHk1 * hk1T1) + axOut.axT1) + (axOut.axRt1 * rt1T1))));
(VS1F[((1.0 * vsStride) + 3.0)] = (zAx * (axOut.axHk1 * hk1D1)));
(VS1F[((1.0 * vsStride) + 4.0)] = (zAx * ((axOut.axHk1 * hk1U1) + (axOut.axRt1 * rt1U1))));
(VS1F[((1.0 * vsStride) + 5.0)] = ax);
(VS2F[((1.0 * vsStride) + 1.0)] = ((zAx * axOut.axA2) + 1.0));
(VS2F[((1.0 * vsStride) + 2.0)] = (zAx * (((axOut.axHk2 * hk2T2) + axOut.axT2) + (axOut.axRt2 * rt2T2))));
(VS2F[((1.0 * vsStride) + 3.0)] = (zAx * (axOut.axHk2 * hk2D2)));
(VS2F[((1.0 * vsStride) + 4.0)] = (zAx * ((axOut.axHk2 * hk2U2) + (axOut.axRt2 * rt2U2))));
(VS2F[((1.0 * vsStride) + 5.0)] = (-ax));
(VSM[1.0] = (zAx * ((((axOut.axHk1 * hk1Ms) + (axOut.axRt1 * rt1Ms)) + (axOut.axHk2 * hk2Ms)) + (axOut.axRt2 * rt2Ms))));
(VSR[1.0] = (zAx * ((axOut.axRt1 * rt1Re) + (axOut.axRt2 * rt2Re))));
(VSX[1.0] = 0.0);
(VSREZ[1.0] = (-rezc));
} else {
double sa = (((1.0 - upw) * s1) + (upw * s2));
double cqa = (((1.0 - upw) * cq1) + (upw * cq2));
double cfa = (((1.0 - upw) * cf1) + (upw * cf2));
double hka = (((1.0 - upw) * hk1) + (upw * hk2));
double usa = (0.5 * (us1 + us2));
double rta = (0.5 * (rt1 + rt2));
double dea = (0.5 * (de1 + de2));
double da = (0.5 * (d1 + d2));
double ald = (truth((ityp == 3.0)) ? dlcon : 1.0);
double gcc = 0.0;
double hkc = (hka - 1.0);
double hkcHka = 1.0;
double hkcRta = 0.0;
if (truth((ityp == 2.0))) {
(gcc = gccon);
(hkc = ((hka - 1.0) - (gcc / rta)));
(hkcHka = 1.0);
(hkcRta = (gcc / host_pow(rta, 2.0)));
if (truth((hkc < 0.01))) {
(hkc = 0.01);
(hkcHka = 0.0);
(hkcRta = 0.0);
}
}
double hr = (hkc / ((gacon * ald) * hka));
double hrHka = ((hkcHka / ((gacon * ald) * hka)) - (hr / hka));
double hrRta = (hkcRta / ((gacon * ald) * hka));
double uq = (((0.5 * cfa) - host_pow(hr, 2.0)) / (gbcon * da));
double uqHka = ((((-2.0) * hr) * hrHka) / (gbcon * da));
double uqRta = ((((-2.0) * hr) * hrRta) / (gbcon * da));
double uqCfa = (0.5 / (gbcon * da));
double uqDa = ((-uq) / da);
double uqUpw = ((uqCfa * (cf2 - cf1)) + (uqHka * (hk2 - hk1)));
double uqT1 = (((1.0 - upw) * ((uqCfa * ctx.get(9)) + (uqHka * ctx.get(111)))) + (uqUpw * upwT1));
double uqD1 = (((1.0 - upw) * ((uqCfa * ctx.get(6)) + (uqHka * ctx.get(109)))) + (uqUpw * upwD1));
double uqU1 = (((1.0 - upw) * ((uqCfa * ctx.get(10)) + (uqHka * ctx.get(112)))) + (uqUpw * upwU1));
double uqT2 = ((upw * ((uqCfa * ctx.get(18)) + (uqHka * ctx.get(116)))) + (uqUpw * upwT2));
double uqD2 = ((upw * ((uqCfa * ctx.get(12)) + (uqHka * ctx.get(114)))) + (uqUpw * upwD2));
double uqU2 = ((upw * ((uqCfa * ctx.get(19)) + (uqHka * ctx.get(117)))) + (uqUpw * upwU2));
double uqMs = ((((1.0 - upw) * ((uqCfa * ctx.get(7)) + (uqHka * ctx.get(110)))) + (uqUpw * upwMs)) + (upw * ((uqCfa * ctx.get(15)) + (uqHka * ctx.get(115)))));
double uqRe = ((((1.0 - upw) * uqCfa) * ctx.get(8)) + ((upw * uqCfa) * ctx.get(16)));
(uqT1 += ((0.5 * uqRta) * ctx.get(173)));
(uqD1 += (0.5 * uqDa));
(uqU1 += ((0.5 * uqRta) * ctx.get(174)));
(uqT2 += ((0.5 * uqRta) * ctx.get(178)));
(uqD2 += (0.5 * uqDa));
(uqU2 += ((0.5 * uqRta) * ctx.get(179)));
(uqMs += (((0.5 * uqRta) * ctx.get(171)) + ((0.5 * uqRta) * ctx.get(176))));
(uqRe += (((0.5 * uqRta) * ctx.get(172)) + ((0.5 * uqRta) * ctx.get(177))));
double scc = ((ctx.get(182) * 1.333) / (1.0 + usa));
double sccUsa = ((-scc) / (1.0 + usa));
double sccUs1 = (sccUsa * 0.5);
double sccUs2 = (sccUsa * 0.5);
double slog = host_log((s2 / s1));
double dxi = (x2 - x1);
double rezc = ((((scc * (cqa - (sa * ald))) * dxi) - ((dea * 2.0) * slog)) + (((dea * 2.0) * ((uq * dxi) - ulog)) * ctx.get(80)));
double zCfa = ((((dea * 2.0) * uqCfa) * dxi) * ctx.get(80));
double zHka = ((((dea * 2.0) * uqHka) * dxi) * ctx.get(80));
double zDa = ((((dea * 2.0) * uqDa) * dxi) * ctx.get(80));
double zSl = ((-dea) * 2.0);
double zUl = (((-dea) * 2.0) * ctx.get(80));
double zDxi = ((scc * (cqa - (sa * ald))) + (((dea * 2.0) * uq) * ctx.get(80)));
double zUsa = ((sccUsa * (cqa - (sa * ald))) * dxi);
double zCqa = (scc * dxi);
double zSa = (((-scc) * dxi) * ald);
double zDea = (2.0 * ((((uq * dxi) - ulog) * ctx.get(80)) - slog));
double zUpw = ((((zCqa * (ctx.get(42) - ctx.get(36))) + (zSa * (ctx.get(181) - ctx.get(180)))) + (zCfa * (ctx.get(11) - ctx.get(5)))) + (zHka * (ctx.get(113) - ctx.get(108))));
double zDe1 = (0.5 * zDea);
double zDe2 = (0.5 * zDea);
double zUs1 = (0.5 * zUsa);
double zUs2 = (0.5 * zUsa);
double zD1 = (0.5 * zDa);
double zD2 = (0.5 * zDa);
double zU1 = ((-zUl) / ctx.get(200));
double zU2 = (zUl / ctx.get(203));
double zX1 = (-zDxi);
double zX2 = zDxi;
double zS1 = (((1.0 - upw) * zSa) - (zSl / ctx.get(180)));
double zS2 = ((upw * zSa) + (zSl / ctx.get(181)));
double zCq1 = ((1.0 - upw) * zCqa);
double zCq2 = (upw * zCqa);
double zCf1 = ((1.0 - upw) * zCfa);
double zCf2 = (upw * zCfa);
double zHk1 = ((1.0 - upw) * zHka);
double zHk2 = (upw * zHka);
(VS1F[((1.0 * vsStride) + 1.0)] = zS1);
(VS1F[((1.0 * vsStride) + 2.0)] = (((zUpw * upwT1) + (zDe1 * ctx.get(57))) + (zUs1 * ctx.get(210))));
(VS1F[((1.0 * vsStride) + 3.0)] = (((zD1 + (zUpw * upwD1)) + (zDe1 * ctx.get(55))) + (zUs1 * ctx.get(207))));
(VS1F[((1.0 * vsStride) + 4.0)] = (((zU1 + (zUpw * upwU1)) + (zDe1 * ctx.get(58))) + (zUs1 * ctx.get(211))));
(VS1F[((1.0 * vsStride) + 5.0)] = zX1);
(VS2F[((1.0 * vsStride) + 1.0)] = zS2);
(VS2F[((1.0 * vsStride) + 2.0)] = (((zUpw * upwT2) + (zDe2 * ctx.get(62))) + (zUs2 * ctx.get(216))));
(VS2F[((1.0 * vsStride) + 3.0)] = (((zD2 + (zUpw * upwD2)) + (zDe2 * ctx.get(60))) + (zUs2 * ctx.get(213))));
(VS2F[((1.0 * vsStride) + 4.0)] = (((zU2 + (zUpw * upwU2)) + (zDe2 * ctx.get(63))) + (zUs2 * ctx.get(217))));
(VS2F[((1.0 * vsStride) + 5.0)] = zX2);
(ctx.arr<1>(46)[1.0] = (((((zUpw * upwMs) + (zDe1 * ctx.get(56))) + (zUs1 * ctx.get(208))) + (zDe2 * ctx.get(61))) + (zUs2 * ctx.get(214))));
(VS1F[((1.0 * vsStride) + 2.0)] += (((zCq1 * ctx.get(40)) + (zCf1 * ctx.get(9))) + (zHk1 * ctx.get(111))));
(VS1F[((1.0 * vsStride) + 3.0)] += (((zCq1 * ctx.get(37)) + (zCf1 * ctx.get(6))) + (zHk1 * ctx.get(109))));
(VS1F[((1.0 * vsStride) + 4.0)] += (((zCq1 * ctx.get(41)) + (zCf1 * ctx.get(10))) + (zHk1 * ctx.get(112))));
(VS2F[((1.0 * vsStride) + 2.0)] += (((zCq2 * ctx.get(46)) + (zCf2 * ctx.get(18))) + (zHk2 * ctx.get(116))));
(VS2F[((1.0 * vsStride) + 3.0)] += (((zCq2 * ctx.get(43)) + (zCf2 * ctx.get(12))) + (zHk2 * ctx.get(114))));
(VS2F[((1.0 * vsStride) + 4.0)] += (((zCq2 * ctx.get(47)) + (zCf2 * ctx.get(19))) + (zHk2 * ctx.get(117))));
(ctx.arr<1>(46)[1.0] += ((((((zCq1 * ctx.get(38)) + (zCf1 * ctx.get(7))) + (zHk1 * ctx.get(110))) + (zCq2 * ctx.get(44))) + (zCf2 * ctx.get(15))) + (zHk2 * ctx.get(115))));
(ctx.arr<1>(47)[1.0] = ((((zCq1 * ctx.get(39)) + (zCf1 * ctx.get(8))) + (zCq2 * ctx.get(45))) + (zCf2 * ctx.get(16))));
(ctx.arr<1>(49)[1.0] = 0.0);
(ctx.arr<1>(48)[1.0] = (-rezc));
}
double ha = (0.5 * (ctx.get(91) + ctx.get(94)));
double ma = (0.5 * (ctx.get(139) + ctx.get(142)));
double xa = (0.5 * (ctx.get(228) + ctx.get(229)));
double ta = (0.5 * (ctx.get(188) + ctx.get(189)));
double hwa = (0.5 * ((ctx.get(81) / ctx.get(188)) + (ctx.get(82) / ctx.get(189))));
double cfx = ((((0.50 * ctx.get(21)) * xa) / ta) + (0.25 * (((ctx.get(5) * ctx.get(228)) / ctx.get(188)) + ((ctx.get(11) * ctx.get(229)) / ctx.get(189)))));
double cfxXa = ((0.50 * ctx.get(21)) / ta);
double cfxTa = ((((-0.50) * ctx.get(21)) * xa) / host_pow(ta, 2.0));
double cfxX1 = (((0.25 * ctx.get(5)) / ctx.get(188)) + (cfxXa * 0.5));
double cfxX2 = (((0.25 * ctx.get(11)) / ctx.get(189)) + (cfxXa * 0.5));
double cfxT1 = (((((-0.25) * ctx.get(5)) * ctx.get(228)) / host_pow(ctx.get(188), 2.0)) + (cfxTa * 0.5));
double cfxT2 = (((((-0.25) * ctx.get(11)) * ctx.get(229)) / host_pow(ctx.get(189), 2.0)) + (cfxTa * 0.5));
double cfxCf1 = ((0.25 * ctx.get(228)) / ctx.get(188));
double cfxCf2 = ((0.25 * ctx.get(229)) / ctx.get(189));
double cfxCfm = ((0.50 * xa) / ta);
double btmp = (((ha + 2.0) - ma) + hwa);
double rezt = ((tlog + (btmp * ulog)) - ((xlog * 0.5) * cfx));
double zCfx = ((-xlog) * 0.5);
double zHa = ulog;
double zHwa = ulog;
double zMa = (-ulog);
double zXl = (((-ddlog) * 0.5) * cfx);
double zUl = (ddlog * btmp);
double zTl = ddlog;
double zCfm = (zCfx * cfxCfm);
double zCf1 = (zCfx * cfxCf1);
double zCf2 = (zCfx * cfxCf2);
double zT1 = ((((-zTl) / ctx.get(188)) + (zCfx * cfxT1)) + ((zHwa * 0.5) * ((-ctx.get(81)) / host_pow(ctx.get(188), 2.0))));
double zT2 = (((zTl / ctx.get(189)) + (zCfx * cfxT2)) + ((zHwa * 0.5) * ((-ctx.get(82)) / host_pow(ctx.get(189), 2.0))));
double zX1 = (((-zXl) / ctx.get(228)) + (zCfx * cfxX1));
double zX2 = ((zXl / ctx.get(229)) + (zCfx * cfxX2));
double zU1 = ((-zUl) / ctx.get(200));
double zU2 = (zUl / ctx.get(203));
(VS1F[((2.0 * vsStride) + 2.0)] = (((((0.5 * zHa) * ctx.get(93)) + (zCfm * ctx.get(29))) + (zCf1 * ctx.get(9))) + zT1));
(VS1F[((2.0 * vsStride) + 3.0)] = ((((0.5 * zHa) * ctx.get(92)) + (zCfm * ctx.get(22))) + (zCf1 * ctx.get(6))));
(VS1F[((2.0 * vsStride) + 4.0)] = (((((0.5 * zMa) * ctx.get(141)) + (zCfm * ctx.get(31))) + (zCf1 * ctx.get(10))) + zU1));
(VS1F[((2.0 * vsStride) + 5.0)] = zX1);
(VS2F[((2.0 * vsStride) + 2.0)] = (((((0.5 * zHa) * ctx.get(96)) + (zCfm * ctx.get(30))) + (zCf2 * ctx.get(18))) + zT2));
(VS2F[((2.0 * vsStride) + 3.0)] = ((((0.5 * zHa) * ctx.get(95)) + (zCfm * ctx.get(23))) + (zCf2 * ctx.get(12))));
(VS2F[((2.0 * vsStride) + 4.0)] = (((((0.5 * zMa) * ctx.get(144)) + (zCfm * ctx.get(32))) + (zCf2 * ctx.get(19))) + zU2));
(VS2F[((2.0 * vsStride) + 5.0)] = zX2);
(ctx.arr<1>(46)[2.0] = ((((((0.5 * zMa) * ctx.get(140)) + (zCfm * ctx.get(26))) + (zCf1 * ctx.get(7))) + ((0.5 * zMa) * ctx.get(143))) + (zCf2 * ctx.get(15))));
(ctx.arr<1>(47)[2.0] = (((zCfm * ctx.get(27)) + (zCf1 * ctx.get(8))) + (zCf2 * ctx.get(16))));
(ctx.arr<1>(49)[2.0] = 0.0);
(ctx.arr<1>(48)[2.0] = (-rezt));
double xot1 = (ctx.get(228) / ctx.get(188));
double xot2 = (ctx.get(229) / ctx.get(189));
double ha2 = (0.5 * (ctx.get(91) + ctx.get(94)));
double hsa = (0.5 * (ctx.get(118) + ctx.get(124)));
double hca = (0.5 * (ctx.get(98) + ctx.get(103)));
double hwa2 = (0.5 * ((ctx.get(81) / ctx.get(188)) + (ctx.get(82) / ctx.get(189))));
double dix = ((((1.0 - upw) * ctx.get(64)) * xot1) + ((upw * ctx.get(71)) * xot2));
double cfx2 = ((((1.0 - upw) * ctx.get(5)) * xot1) + ((upw * ctx.get(11)) * xot2));
double dixUpw = ((ctx.get(71) * xot2) - (ctx.get(64) * xot1));
double cfxUpw = ((ctx.get(11) * xot2) - (ctx.get(5) * xot1));
double btmp2 = (((((2.0 * hca) / hsa) + 1.0) - ha2) - hwa2);
double rezh = ((hlog + (btmp2 * ulog)) + (xlog * ((0.5 * cfx2) - dix)));
double zCfx2 = (xlog * 0.5);
double zDix = (-xlog);
double zHca = ((2.0 * ulog) / hsa);
double zHa2 = (-ulog);
double zHwa2 = (-ulog);
double zXl2 = (ddlog * ((0.5 * cfx2) - dix));
double zUl2 = (ddlog * btmp2);
double zHl = ddlog;
double zUpw2 = ((zCfx2 * cfxUpw) + (zDix * dixUpw));
double zHs1 = ((((-hca) * ulog) / host_pow(hsa, 2.0)) - (zHl / ctx.get(118)));
double zHs2 = ((((-hca) * ulog) / host_pow(hsa, 2.0)) + (zHl / ctx.get(124)));
double zCf1b = (((1.0 - upw) * zCfx2) * xot1);
double zCf2b = ((upw * zCfx2) * xot2);
double zDi1 = (((1.0 - upw) * zDix) * xot1);
double zDi2 = ((upw * zDix) * xot2);
double zT1b = (((1.0 - upw) * ((zCfx2 * ctx.get(5)) + (zDix * ctx.get(64)))) * ((-xot1) / ctx.get(188)));
double zT2b = ((upw * ((zCfx2 * ctx.get(11)) + (zDix * ctx.get(71)))) * ((-xot2) / ctx.get(189)));
double zX1b = ((((1.0 - upw) * ((zCfx2 * ctx.get(5)) + (zDix * ctx.get(64)))) / ctx.get(188)) - (zXl2 / ctx.get(228)));
double zX2b = (((upw * ((zCfx2 * ctx.get(11)) + (zDix * ctx.get(71)))) / ctx.get(189)) + (zXl2 / ctx.get(229)));
double zU1b = ((-zUl2) / ctx.get(200));
double zU2b = (zUl2 / ctx.get(203));
(zT1b += ((zHwa2 * 0.5) * ((-ctx.get(81)) / host_pow(ctx.get(188), 2.0))));
(zT2b += ((zHwa2 * 0.5) * ((-ctx.get(82)) / host_pow(ctx.get(189), 2.0))));
(VS1F[((3.0 * vsStride) + 1.0)] = (zDi1 * ctx.get(68)));
(VS1F[((3.0 * vsStride) + 2.0)] = ((((zHs1 * ctx.get(122)) + (zCf1b * ctx.get(9))) + (zDi1 * ctx.get(69))) + zT1b));
(VS1F[((3.0 * vsStride) + 3.0)] = (((zHs1 * ctx.get(119)) + (zCf1b * ctx.get(6))) + (zDi1 * ctx.get(65))));
(VS1F[((3.0 * vsStride) + 4.0)] = ((((zHs1 * ctx.get(123)) + (zCf1b * ctx.get(10))) + (zDi1 * ctx.get(70))) + zU1b));
(VS1F[((3.0 * vsStride) + 5.0)] = zX1b);
(VS2F[((3.0 * vsStride) + 1.0)] = (zDi2 * ctx.get(75)));
(VS2F[((3.0 * vsStride) + 2.0)] = ((((zHs2 * ctx.get(128)) + (zCf2b * ctx.get(18))) + (zDi2 * ctx.get(76))) + zT2b));
(VS2F[((3.0 * vsStride) + 3.0)] = (((zHs2 * ctx.get(125)) + (zCf2b * ctx.get(12))) + (zDi2 * ctx.get(72))));
(VS2F[((3.0 * vsStride) + 4.0)] = ((((zHs2 * ctx.get(129)) + (zCf2b * ctx.get(19))) + (zDi2 * ctx.get(77))) + zU2b));
(VS2F[((3.0 * vsStride) + 5.0)] = zX2b);
(ctx.arr<1>(46)[3.0] = ((((((zHs1 * ctx.get(120)) + (zCf1b * ctx.get(7))) + (zDi1 * ctx.get(66))) + (zHs2 * ctx.get(126))) + (zCf2b * ctx.get(15))) + (zDi2 * ctx.get(73))));
(ctx.arr<1>(47)[3.0] = ((((((zHs1 * ctx.get(121)) + (zCf1b * ctx.get(8))) + (zDi1 * ctx.get(67))) + (zHs2 * ctx.get(127))) + (zCf2b * ctx.get(16))) + (zDi2 * ctx.get(74))));
(VS1F[((3.0 * vsStride) + 2.0)] += ((0.5 * ((zHca * ctx.get(101)) + (zHa2 * ctx.get(93)))) + (zUpw2 * upwT1)));
(VS1F[((3.0 * vsStride) + 3.0)] += ((0.5 * ((zHca * ctx.get(99)) + (zHa2 * ctx.get(92)))) + (zUpw2 * upwD1)));
(VS1F[((3.0 * vsStride) + 4.0)] += ((0.5 * (zHca * ctx.get(102))) + (zUpw2 * upwU1)));
(VS2F[((3.0 * vsStride) + 2.0)] += ((0.5 * ((zHca * ctx.get(106)) + (zHa2 * ctx.get(96)))) + (zUpw2 * upwT2)));
(VS2F[((3.0 * vsStride) + 3.0)] += ((0.5 * ((zHca * ctx.get(104)) + (zHa2 * ctx.get(95)))) + (zUpw2 * upwD2)));
(VS2F[((3.0 * vsStride) + 4.0)] += ((0.5 * (zHca * ctx.get(107))) + (zUpw2 * upwU2)));
(ctx.arr<1>(46)[3.0] += (((0.5 * (zHca * ctx.get(100))) + (zUpw2 * upwMs)) + (0.5 * (zHca * ctx.get(105)))));
(ctx.arr<1>(49)[3.0] = 0.0);
(ctx.arr<1>(48)[3.0] = (-rezh));
}

// js/xblsys.js:2197 (dampl)
Record fn_dampl(double hk, double th, double rt) {
double dgr = 0.08;
double hmi = (1.0 / (hk - 1.0));
double hmiHk = (-host_pow(hmi, 2.0));
double aa = (2.492 * host_pow(hmi, 0.43));
double aaHk = (((aa / hmi) * 0.43) * hmiHk);
double bb = host_tanh(((14.0 * hmi) - 9.24));
double bbHk = (((1.0 - (bb * bb)) * 14.0) * hmiHk);
double grcrit = (aa + (0.7 * (bb + 1.0)));
double grcHk = (aaHk + (0.7 * bbHk));
double gr = host_log10(rt);
double grRt = (1.0 / (2.3025851 * rt));
if (truth((gr < (grcrit - dgr)))) {
return ([&]() { Record value{}; value.ax = 0.0; value.axHk = 0.0; value.axTh = 0.0; value.axRt = 0.0; return value; }());
}
double rnorm = ((gr - (grcrit - dgr)) / (2.0 * dgr));
double rnHk = ((-grcHk) / (2.0 * dgr));
double rnRt = (grRt / (2.0 * dgr));
double rfac = nan_value();
double rfacHk = nan_value();
double rfacRt = nan_value();
if (truth((rnorm >= 1.0))) {
(rfac = 1.0);
(rfacHk = 0.0);
(rfacRt = 0.0);
} else {
(rfac = ((3.0 * host_pow(rnorm, 2.0)) - (2.0 * host_pow(rnorm, 3.0))));
double rfacRn = ((6.0 * rnorm) - (6.0 * host_pow(rnorm, 2.0)));
(rfacHk = (rfacRn * rnHk));
(rfacRt = (rfacRn * rnRt));
}
double arg = ((3.87 * hmi) - 2.52);
double argHk = (3.87 * hmiHk);
double ex = host_exp((-host_pow(arg, 2.0)));
double exHk = (ex * (((-2.0) * arg) * argHk));
double dadr = ((0.028 * (hk - 1.0)) - (0.0345 * ex));
double dadrHk = (0.028 - (0.0345 * exHk));
double af = ((((-0.05) + (2.7 * hmi)) - (5.5 * host_pow(hmi, 2.0))) + (3.0 * host_pow(hmi, 3.0)));
double afHmi = ((2.7 - (11.0 * hmi)) + (9.0 * host_pow(hmi, 2.0)));
double afHk = (afHmi * hmiHk);
double ax = (((af * dadr) / th) * rfac);
double axHk = (((((afHk * dadr) / th) + ((af * dadrHk) / th)) * rfac) + (((af * dadr) / th) * rfacHk));
double axTh = ((-ax) / th);
double axRt = (((af * dadr) / th) * rfacRt);
return ([&]() { Record value{}; value.ax = ax; value.axHk = axHk; value.axTh = axTh; value.axRt = axRt; return value; }());
}

// js/xblsys.js:2249 (dampl2)
Record fn_dampl2(double hk, double th, double rt) {
double dgr = 0.08;
double hk1 = 3.5;
double hk2 = 4.0;
double hmi = (1.0 / (hk - 1.0));
double hmiHk = (-host_pow(hmi, 2.0));
double aa = (2.492 * host_pow(hmi, 0.43));
double aaHk = (((aa / hmi) * 0.43) * hmiHk);
double bb = host_tanh(((14.0 * hmi) - 9.24));
double bbHk = (((1.0 - (bb * bb)) * 14.0) * hmiHk);
double grc = (aa + (0.7 * (bb + 1.0)));
double grcHk = (aaHk + (0.7 * bbHk));
double gr = host_log10(rt);
double grRt = (1.0 / (2.3025851 * rt));
if (truth((gr < (grc - dgr)))) {
return ([&]() { Record value{}; value.ax = 0.0; value.axHk = 0.0; value.axTh = 0.0; value.axRt = 0.0; return value; }());
}
double rnorm = ((gr - (grc - dgr)) / (2.0 * dgr));
double rnHk = ((-grcHk) / (2.0 * dgr));
double rnRt = (grRt / (2.0 * dgr));
double rfac = nan_value();
double rfacHk = nan_value();
double rfacRt = nan_value();
if (truth((rnorm >= 1.0))) {
(rfac = 1.0);
(rfacHk = 0.0);
(rfacRt = 0.0);
} else {
(rfac = ((3.0 * host_pow(rnorm, 2.0)) - (2.0 * host_pow(rnorm, 3.0))));
double rfacRn = ((6.0 * rnorm) - (6.0 * host_pow(rnorm, 2.0)));
(rfacHk = (rfacRn * rnHk));
(rfacRt = (rfacRn * rnRt));
}
double arg = ((3.87 * hmi) - 2.52);
double argHk = (3.87 * hmiHk);
double ex = host_exp((-host_pow(arg, 2.0)));
double exHk = (ex * (((-2.0) * arg) * argHk));
double dadr = ((0.028 * (hk - 1.0)) - (0.0345 * ex));
double dadrHk = (0.028 - (0.0345 * exHk));
double brg = ((-20.0) * hmi);
double af = (((((-0.05) + (2.7 * hmi)) - (5.5 * host_pow(hmi, 2.0))) + (3.0 * host_pow(hmi, 3.0))) + (0.1 * host_exp(brg)));
double afHmi = (((2.7 - (11.0 * hmi)) + (9.0 * host_pow(hmi, 2.0))) - (2.0 * host_exp(brg)));
double afHk = (afHmi * hmiHk);
double ax = (((af * dadr) / th) * rfac);
double axHk = (((((afHk * dadr) / th) + ((af * dadrHk) / th)) * rfac) + (((af * dadr) / th) * rfacHk));
double axTh = ((-ax) / th);
double axRt = (((af * dadr) / th) * rfacRt);
if (truth((hk < hk1))) {
return ([&]() { Record value{}; value.ax = ax; value.axHk = axHk; value.axTh = axTh; value.axRt = axRt; return value; }());
}
double hnorm = ((hk - hk1) / (hk2 - hk1));
double hnHk = (1.0 / (hk2 - hk1));
double hfac = nan_value();
double hfHk = nan_value();
if (truth((hnorm >= 1.0))) {
(hfac = 1.0);
(hfHk = 0.0);
} else {
(hfac = ((3.0 * host_pow(hnorm, 2.0)) - (2.0 * host_pow(hnorm, 3.0))));
(hfHk = (((6.0 * hnorm) - (6.0 * host_pow(hnorm, 2.0))) * hnHk));
}
double ax1 = ax;
double ax1Hk = axHk;
double ax1Th = axTh;
double ax1Rt = axRt;
double gr0 = (0.30 + (0.35 * host_exp(((-0.15) * (hk - 5.0)))));
double gr0Hk = (((-0.35) * host_exp(((-0.15) * (hk - 5.0)))) * 0.15);
double tnr = host_tanh((1.2 * (gr - gr0)));
double tnrRt = (((1.0 - host_pow(tnr, 2.0)) * 1.2) * grRt);
double tnrHk = (((-(1.0 - host_pow(tnr, 2.0))) * 1.2) * gr0Hk);
double ax2 = (((0.086 * tnr) - (0.25 / host_pow((hk - 1.0), 1.5))) / th);
double ax2Hk = (((0.086 * tnrHk) + ((1.5 * 0.25) / host_pow((hk - 1.0), 2.5))) / th);
double ax2Rt = ((0.086 * tnrRt) / th);
double ax2Th = ((-ax2) / th);
if (truth((ax2 < 0.0))) {
(ax2 = 0.0);
(ax2Hk = 0.0);
(ax2Rt = 0.0);
(ax2Th = 0.0);
}
(ax = ((hfac * ax2) + ((1.0 - hfac) * ax1)));
(axHk = (((hfac * ax2Hk) + ((1.0 - hfac) * ax1Hk)) + (hfHk * (ax2 - ax1))));
(axRt = ((hfac * ax2Rt) + ((1.0 - hfac) * ax1Rt)));
(axTh = ((hfac * ax2Th) + ((1.0 - hfac) * ax1Th)));
return ([&]() { Record value{}; value.ax = ax; value.axHk = axHk; value.axTh = axTh; value.axRt = axRt; return value; }());
}

// js/xblsys.js:2351 (hkin)
Record fn_hkin(double h, double msq) {
double hk = ((h - (0.29 * msq)) / (1.0 + (0.113 * msq)));
double hkH = (1.0 / (1.0 + (0.113 * msq)));
double hkMsq = (((-0.29) - (0.113 * hk)) / (1.0 + (0.113 * msq)));
return ([&]() { Record value{}; value.hk = hk; value.hkH = hkH; value.hkMsq = hkMsq; return value; }());
}

// js/xblsys.js:2359 (dil)
Record fn_dil(double hk, double rt) {
double di = nan_value();
double diHk = nan_value();
if (truth((hk < 4.0))) {
(di = (((0.00205 * host_pow((4.0 - hk), 5.5)) + 0.207) / rt));
(diHk = ((((-0.00205) * 5.5) * host_pow((4.0 - hk), 4.5)) / rt));
} else {
double hkb = (hk - 4.0);
double den = (1.0 + (0.02 * host_pow(hkb, 2.0)));
(di = (((((-0.0016) * host_pow(hkb, 2.0)) / den) + 0.207) / rt));
(diHk = (((((-0.0016) * 2.0) * hkb) * ((1.0 / den) - ((0.02 * host_pow(hkb, 2.0)) / host_pow(den, 2.0)))) / rt));
}
double diRt = ((-di) / rt);
return ([&]() { Record value{}; value.di = di; value.diHk = diHk; value.diRt = diRt; return value; }());
}

// js/xblsys.js:2376 (dilw)
Record fn_dilw(double hk, double rt) {
double msq = 0.0;
Record hslRes = fn_hsl(hk, rt, msq);
double rcd = ((1.10 * host_pow((1.0 - (1.0 / hk)), 2.0)) / hk);
double rcdHk = (((((-1.10) * (1.0 - (1.0 / hk))) * 2.0) / host_pow(hk, 3.0)) - (rcd / hk));
double di = ((2.0 * rcd) / (hslRes.hs * rt));
double diHk = (((2.0 * rcdHk) / (hslRes.hs * rt)) - ((di / hslRes.hs) * hslRes.hsHk));
double diRt = (((-di) / rt) - ((di / hslRes.hs) * hslRes.hsRt));
return ([&]() { Record value{}; value.di = di; value.diHk = diHk; value.diRt = diRt; return value; }());
}

// js/xblsys.js:2388 (hsl)
Record fn_hsl(double hk, double rt, double msq) {
double hs = nan_value();
double hsHk = nan_value();
if (truth((hk < 4.35))) {
double tmp = (hk - 4.35);
(hs = (((((0.0111 * host_pow(tmp, 2.0)) / (hk + 1.0)) - ((0.0278 * host_pow(tmp, 3.0)) / (hk + 1.0))) + 1.528) - (0.0002 * host_pow((tmp * hk), 2.0))));
(hsHk = ((((0.0111 * ((2.0 * tmp) - (host_pow(tmp, 2.0) / (hk + 1.0)))) / (hk + 1.0)) - ((0.0278 * ((3.0 * host_pow(tmp, 2.0)) - (host_pow(tmp, 3.0) / (hk + 1.0)))) / (hk + 1.0))) - ((((0.0002 * 2.0) * tmp) * hk) * (tmp + hk))));
} else {
double hs2 = 0.015;
(hs = (((hs2 * host_pow((hk - 4.35), 2.0)) / hk) + 1.528));
(hsHk = ((((hs2 * 2.0) * (hk - 4.35)) / hk) - ((hs2 * host_pow((hk - 4.35), 2.0)) / host_pow(hk, 2.0))));
}
return ([&]() { Record value{}; value.hs = hs; value.hsHk = hsHk; value.hsRt = 0.0; value.hsMsq = 0.0; return value; }());
}

// js/xblsys.js:2408 (cfl)
Record fn_cfl(double hk, double rt, double msq) {
double cf = nan_value();
double cfHk = nan_value();
if (truth((hk < 5.5))) {
double tmp = (host_pow((5.5 - hk), 3.0) / (hk + 1.0));
(cf = (((0.0727 * tmp) - 0.07) / rt));
(cfHk = ((((((-0.0727) * tmp) * 3.0) / (5.5 - hk)) - ((0.0727 * tmp) / (hk + 1.0))) / rt));
} else {
double tmp = (1.0 - (1.0 / (hk - 4.5)));
(cf = (((0.015 * host_pow(tmp, 2.0)) - 0.07) / rt));
(cfHk = ((((0.015 * tmp) * 2.0) / host_pow((hk - 4.5), 2.0)) / rt));
}
double cfRt = ((-cf) / rt);
double cfMsq = 0.0;
return ([&]() { Record value{}; value.cf = cf; value.cfHk = cfHk; value.cfRt = cfRt; value.cfMsq = cfMsq; return value; }());
}

// js/xblsys.js:2426 (dit)
Record fn_dit(double hs, double us, double cf, double st) {
double di = (((((0.5 * cf) * us) + ((st * st) * (1.0 - us))) * 2.0) / hs);
double diHs = (((-(((0.5 * cf) * us) + ((st * st) * (1.0 - us)))) * 2.0) / host_pow(hs, 2.0));
double diUs = ((((0.5 * cf) - (st * st)) * 2.0) / hs);
double diCf = (((0.5 * us) * 2.0) / hs);
double diSt = ((((2.0 * st) * (1.0 - us)) * 2.0) / hs);
return ([&]() { Record value{}; value.di = di; value.diHs = diHs; value.diUs = diUs; value.diCf = diCf; value.diSt = diSt; return value; }());
}

// js/xblsys.js:2436 (hst)
Record fn_hst(double hk, double rt, double msq) {
double hsmin = 1.5;
double dhsinf = 0.015;
double ho = nan_value();
double hoRt = nan_value();
if (truth((rt > 400.0))) {
(ho = (3.0 + (400.0 / rt)));
(hoRt = ((-400.0) / host_pow(rt, 2.0)));
} else {
(ho = 4.0);
(hoRt = 0.0);
}
double rtz = nan_value();
double rtzRt = nan_value();
if (truth((rt > 200.0))) {
(rtz = rt);
(rtzRt = 1.0);
} else {
(rtz = 200.0);
(rtzRt = 0.0);
}
double hs = nan_value();
double hsHk = nan_value();
double hsRt = nan_value();
if (truth((hk < ho))) {
double hr = ((ho - hk) / (ho - 1.0));
double hrHk = ((-1.0) / (ho - 1.0));
double hrRt = (((1.0 - hr) / (ho - 1.0)) * hoRt);
(hs = (((((((2.0 - hsmin) - (4.0 / rtz)) * host_pow(hr, 2.0)) * 1.5) / (hk + 0.5)) + hsmin) + (4.0 / rtz)));
(hsHk = (((((-((2.0 - hsmin) - (4.0 / rtz))) * host_pow(hr, 2.0)) * 1.5) / host_pow((hk + 0.5), 2.0)) + (((((((2.0 - hsmin) - (4.0 / rtz)) * hr) * 2.0) * 1.5) / (hk + 0.5)) * hrHk)));
(hsRt = ((((((((2.0 - hsmin) - (4.0 / rtz)) * hr) * 2.0) * 1.5) / (hk + 0.5)) * hrRt) + ((((((host_pow(hr, 2.0) * 1.5) / (hk + 0.5)) - 1.0) * 4.0) / host_pow(rtz, 2.0)) * rtzRt)));
} else {
double grt = host_log(rtz);
double hdif = (hk - ho);
double rtmp = ((hk - ho) + (4.0 / grt));
double htmp = (((0.007 * grt) / host_pow(rtmp, 2.0)) + (dhsinf / hk));
double htmpHk = ((((-0.014) * grt) / host_pow(rtmp, 3.0)) - (dhsinf / host_pow(hk, 2.0)));
double htmpRt = (((((-0.014) * grt) / host_pow(rtmp, 3.0)) * ((-hoRt) - (((4.0 / host_pow(grt, 2.0)) / rtz) * rtzRt))) + (((0.007 / host_pow(rtmp, 2.0)) / rtz) * rtzRt));
(hs = (((host_pow(hdif, 2.0) * htmp) + hsmin) + (4.0 / rtz)));
(hsHk = (((hdif * 2.0) * htmp) + (host_pow(hdif, 2.0) * htmpHk)));
(hsRt = (((host_pow(hdif, 2.0) * htmpRt) - ((4.0 / host_pow(rtz, 2.0)) * rtzRt)) + (((hdif * 2.0) * htmp) * (-hoRt))));
}
double fm = (1.0 + (0.014 * msq));
(hs = ((hs + (0.028 * msq)) / fm));
(hsHk /= fm);
(hsRt /= fm);
double hsMsq = ((0.028 / fm) - ((0.014 * hs) / fm));
return ([&]() { Record value{}; value.hs = hs; value.hsHk = hsHk; value.hsRt = hsRt; value.hsMsq = hsMsq; return value; }());
}

// js/xblsys.js:2493 (cft)
Record fn_cft(double hk, double rt, double msq, double cffac) {
double gam = 1.4;
double gm1 = (gam - 1.0);
double fc = __builtin_sqrt((1.0 + ((0.5 * gm1) * msq)));
double grt = host_log((rt / fc));
double grtRt = (1.0 / rt);
double grtMsq = (((-0.25) * gm1) / host_pow(fc, 2.0));
if (truth((grt < 3.0))) {
(grt = 3.0);
(grtRt = 0.0);
(grtMsq = 0.0);
}
double gex = ((-1.74) - (0.31 * hk));
double arg = ((-1.33) * hk);
(arg = js_max((-20.0), arg));
double thk = host_tanh((4.0 - (hk / 0.875)));
double cfo = (((cffac * 0.3) * host_exp(arg)) * host_pow((grt / 2.3026), gex));
double cf = ((cfo + (1.1e-4 * (thk - 1.0))) / fc);
double cfHk = (((((-1.33) * cfo) - ((0.31 * host_log((grt / 2.3026))) * cfo)) - ((1.1e-4 * (1.0 - host_pow(thk, 2.0))) / 0.875)) / fc);
double cfRt = (((gex * cfo) / (fc * grt)) * grtRt);
double cfMsq = ((((gex * cfo) / (fc * grt)) * grtMsq) - (((0.25 * gm1) * cf) / host_pow(fc, 2.0)));
return ([&]() { Record value{}; value.cf = cf; value.cfHk = cfHk; value.cfRt = cfRt; value.cfMsq = cfMsq; return value; }());
}

// js/xblsys.js:2519 (hct)
Record fn_hct(double hk, double msq) {
double hc = (msq * ((0.064 / (hk - 0.8)) + 0.251));
double hcHk = (msq * ((-0.064) / host_pow((hk - 0.8), 2.0)));
double hcMsq = ((0.064 / (hk - 0.8)) + 0.251);
return ([&]() { Record value{}; value.hc = hc; value.hcHk = hcHk; value.hcMsq = hcMsq; return value; }());
}

// js/xbl.js:16 (prepareBlStation)
void fn_prepareBlStation(Context ctx, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei) {
if ((truth((fn_blprv == fn_blprv)) && truth((fn_blkin == fn_blkin)))) {
fn_blprvKin(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
} else {
fn_blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_blkin(ctx);
}
}

// js/xbl.js:26 (dslim)
double fn_dslim(Context ctx, double dstr, double thet, double uedg, double msq, double hklim) {
double h = (dstr / thet);
auto record_1240 = fn_hkin(h, msq);
double hk = record_1240.hk;
double hkH = record_1240.hkH;
double dh = (js_max(0.0, (hklim - hk)) / hkH);
return (dstr + (dh * thet));
}

// js/xbl.js:52 (mrcl)
Record fn_mrcl(Context ctx, double cls) {
double cla = js_max(cls, 1.0e-6);
double minf = nan_value();
double mCls = 0.0;
double reinf = nan_value();
double rCls = 0.0;
double matyp = ((truth((ctx.get(145) >= 1.0)) && truth((ctx.get(145) <= 3.0))) ? ctx.get(145) : 1.0);
double retyp = ((truth((ctx.get(164) >= 1.0)) && truth((ctx.get(164) <= 3.0))) ? ctx.get(164) : 1.0);
ctx.put(145, matyp);
ctx.put(164, retyp);
if (truth((matyp == 1.0))) {
(minf = ctx.get(147));
(mCls = 0.0);
} else if (truth((matyp == 2.0))) {
(minf = (ctx.get(147) / __builtin_sqrt(cla)));
(mCls = (((-0.5) * minf) / cla));
} else {
(minf = ctx.get(147));
(mCls = 0.0);
}
if (truth((retyp == 1.0))) {
(reinf = ctx.get(162));
(rCls = 0.0);
} else if (truth((retyp == 2.0))) {
(reinf = (ctx.get(162) / __builtin_sqrt(cla)));
(rCls = (((-0.5) * reinf) / cla));
} else {
(reinf = (ctx.get(162) / cla));
(rCls = ((-reinf) / cla));
}
if (truth((minf >= 0.99))) {
(minf = 0.99);
(mCls = 0.0);
}
if ((truth((ctx.get(162) > 0.0)) && truth(((reinf / ctx.get(162)) > 100.0)))) {
(reinf = (ctx.get(162) * 100.0));
(rCls = 0.0);
}
ctx.put(146, minf);
ctx.put(161, reinf);
ctx.put(148, mCls);
ctx.put(163, rCls);
return ([&]() { Record value{}; value.mCls = mCls; value.rCls = rCls; return value; }());
}

// js/xbl.js:102 (comset)
void fn_comset(Context ctx) {
double beta = __builtin_sqrt((1.0 - host_pow(ctx.get(146), 2.0)));
double betaMsq = ((-0.5) / beta);
ctx.put(192, (host_pow(ctx.get(146), 2.0) / host_pow((1.0 + beta), 2.0)));
ctx.put(194, ((1.0 / host_pow((1.0 + beta), 2.0)) - (((2.0 * ctx.get(192)) / (1.0 + beta)) * betaMsq)));
ctx.put(193, ctx.get(194));
if (truth((ctx.get(146) == 0.0))) {
ctx.put(35, (-999.0));
ctx.put(154, 999.0);
} else {
ctx.put(35, ((2.0 / (ctx.get(87) * host_pow(ctx.get(146), 2.0))) * (host_pow(((1.0 + ((0.5 * ctx.get(86)) * host_pow(ctx.get(146), 2.0))) / (1.0 + (0.5 * ctx.get(86)))), (ctx.get(87) / ctx.get(86))) - 1.0)));
ctx.put(154, ((ctx.get(152) / ctx.get(146)) * __builtin_sqrt(((1.0 + ((0.5 * ctx.get(86)) * host_pow(ctx.get(146), 2.0))) / (1.0 + (0.5 * ctx.get(86)))))));
}
}

// js/xbl.js:120 (iblsys)
void fn_iblsys(Context ctx) {
double iv = 0.0;
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
for (u32 ibl = 2u; truth((double(ibl) <= ctx.arr<1>(22).at(is))); advance_index(ibl)) {
(iv += 1.0);
(ctx.arr<2>(15).at(ibl).at(is) = iv);
}
}
ctx.put(151, iv);
if (truth((ctx.get(151) > (2.0 * ctx.get(136))))) {
unsupported();
}
}

// js/xbl.js:143 (mrchue)
void fn_mrchue(Context ctx) {
((void)0);
double hlmax = 3.8;
double htmax = 2.5;
if ((((truth((!bool(ctx.arr<2>(12)))) || truth((ctx.arr<2>(12).length < 3.0))) || truth((ctx.arr<2>(12)[1.0].length <= ctx.arr<1>(22)[1.0]))) || truth((ctx.arr<2>(12)[2.0].length <= ctx.arr<1>(22)[2.0])))) {
(ctx.arr<2>(12) = make_rows(V1{}, make_vec<1>((ctx.arr<1>(22)[1.0] + 1.0)), make_vec<1>((ctx.arr<1>(22)[2.0] + 1.0))));
} else {
ctx.arr<2>(12)[1.0].fill(0.0);
ctx.arr<2>(12)[2.0].fill(0.0);
}
auto resetThwaites = [&](double xsi, double uei) -> Record {
double ucon = (uei / host_pow(xsi, ctx.get(4)));
double tsq = ((0.45 / ((ucon * ((5.0 * ctx.get(4)) + 1.0)) * ctx.get(165))) * host_pow(xsi, (1.0 - ctx.get(4))));
double thi = __builtin_sqrt(tsq);
double dsi = (2.2 * thi);
return ([&]() { Record value{}; value.thi = thi; value.dsi = dsi; return value; }());
};
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
ctx.put(0, ctx.arr<1>(0).at(is));
fn_xifset(ctx, double(is));
double ibl = 2.0;
double xsi = ctx.arr<2>(60)[ibl].at(is);
double uei = ctx.arr<2>(30)[ibl].at(is);
ctx.put(4, 1.0);
Record init = resetThwaites(xsi, uei);
double thi = init.thi;
double dsi = init.dsi;
double ami = 0.0;
double cti = 0.03;
ctx.put(195, 0.0);
ctx.put(199, 0.0);
(ctx.arr<1>(16).at(is) = ctx.arr<1>(13).at(is));
for ((ibl = 2.0); truth((ibl <= ctx.arr<1>(22).at(is))); (ibl += 1.0)) {
double ibm = (ibl - 1.0);
double simi = (ibl == 2.0);
double wake = (ibl > ctx.arr<1>(13).at(is));
ctx.put(183, simi);
ctx.put(227, wake);
(xsi = ctx.arr<2>(60)[ibl].at(is));
(uei = ctx.arr<2>(30)[ibl].at(is));
double dswaki = 0.0;
if (truth(wake)) {
double iw = (ibl - ctx.arr<1>(13).at(is));
(dswaki = ctx.arr<1>(56)[iw]);
}
double direct = 1.0;
double htarg = 0.0;
double dmax = 0.0;
double msq = 0.0;
double hklim = 0.0;
for (u32 itbl = 1u; truth((double(itbl) <= 25.0)); advance_index(itbl)) {
fn_prepareBlStation(ctx, xsi, ami, cti, thi, dsi, dswaki, uei);
if ((truth((!truth(simi))) && truth((!truth(ctx.get(199)))))) {
fn_trchek(ctx);
(ami = ctx.get(2));
if (truth(ctx.get(195))) {
(ctx.arr<1>(16).at(is) = ibl);
if (truth((cti <= 0.0))) {
(cti = 0.03);
ctx.put(181, cti);
}
} else {
(ctx.arr<1>(16).at(is) = (ibl + 2.0));
}
}
if (truth((ibl == (ctx.arr<1>(13).at(is) + 1.0)))) {
double tte = (ctx.arr<2>(27)[ctx.arr<1>(13)[1.0]][1.0] + ctx.arr<2>(27)[ctx.arr<1>(13)[2.0]][2.0]);
double dte = ((ctx.arr<2>(11)[ctx.arr<1>(13)[1.0]][1.0] + ctx.arr<2>(11)[ctx.arr<1>(13)[2.0]][2.0]) + ctx.get(3));
double cte = (((ctx.arr<2>(5)[ctx.arr<1>(13)[1.0]][1.0] * ctx.arr<2>(27)[ctx.arr<1>(13)[1.0]][1.0]) + (ctx.arr<2>(5)[ctx.arr<1>(13)[2.0]][2.0] * ctx.arr<2>(27)[ctx.arr<1>(13)[2.0]][2.0])) / tte);
ctx.put(49, cte);
ctx.put(198, tte);
ctx.put(79, dte);
fn_tesys(cte, tte, dte, ctx);
} else {
fn_blsys(ctx);
}
if (truth(direct)) {
(ctx.arr<1>(45)[((4.0 * 6.0) + 1.0)] = 0.0);
(ctx.arr<1>(45)[((4.0 * 6.0) + 2.0)] = 0.0);
(ctx.arr<1>(45)[((4.0 * 6.0) + 3.0)] = 0.0);
(ctx.arr<1>(45)[((4.0 * 6.0) + 4.0)] = 1.0);
(ctx.arr<1>(48)[4.0] = 0.0);
fn_gauss1(4.0, ctx.arr<2>(44), ctx.arr<1>(48));
(dmax = js_max(__builtin_fabs((ctx.arr<1>(48)[2.0] / thi)), __builtin_fabs((ctx.arr<1>(48)[3.0] / dsi))));
if (truth((ibl < ctx.arr<1>(16).at(is)))) (dmax = js_max(dmax, __builtin_fabs((ctx.arr<1>(48)[1.0] / 10.0))));
if (truth((ibl >= ctx.arr<1>(16).at(is)))) (dmax = js_max(dmax, __builtin_fabs((ctx.arr<1>(48)[1.0] / cti))));
double rlx = 1.0;
if (truth((dmax > 0.3))) (rlx = (0.3 / dmax));
double hktest = nan_value();
double hmax = nan_value();
if (truth((ibl != (ctx.arr<1>(13).at(is) + 1.0)))) {
(msq = (((uei * uei) * ctx.get(130)) / (ctx.get(90) * (1.0 - (((0.5 * uei) * uei) * ctx.get(130))))));
double htest = ((dsi + (rlx * ctx.arr<1>(48)[3.0])) / (thi + (rlx * ctx.arr<1>(48)[2.0])));
Record hkin = fn_hkin(htest, msq);
(hktest = hkin.hk);
(hmax = htmax);
if (truth((ibl < ctx.arr<1>(16).at(is)))) (hmax = hlmax);
(direct = (hktest < hmax));
}
if (truth(direct)) {
if (truth((ibl >= ctx.arr<1>(16).at(is)))) (cti = (cti + (rlx * ctx.arr<1>(48)[1.0])));
(thi = (thi + (rlx * ctx.arr<1>(48)[2.0])));
(dsi = (dsi + (rlx * ctx.arr<1>(48)[3.0])));
} else {
if (truth((ibl < ctx.arr<1>(16).at(is)))) {
(htarg = (ctx.get(108) + ((0.03 * (ctx.get(229) - ctx.get(228))) / ctx.get(188))));
} else if (truth((ibl == ctx.arr<1>(16).at(is)))) {
(htarg = (ctx.get(108) + (((0.03 * (ctx.get(232) - ctx.get(228))) - (0.15 * (ctx.get(229) - ctx.get(232)))) / ctx.get(188))));
} else if (truth(wake)) {
double constant = ((0.03 * (ctx.get(229) - ctx.get(228))) / ctx.get(188));
double hk2 = ctx.get(108);
(hk2 = (hk2 - (((hk2 + (constant * host_pow((hk2 - 1.0), 3.0))) - ctx.get(108)) / (1.0 + ((3.0 * constant) * host_pow((hk2 - 1.0), 2.0))))));
(hk2 = (hk2 - (((hk2 + (constant * host_pow((hk2 - 1.0), 3.0))) - ctx.get(108)) / (1.0 + ((3.0 * constant) * host_pow((hk2 - 1.0), 2.0))))));
(hk2 = (hk2 - (((hk2 + (constant * host_pow((hk2 - 1.0), 3.0))) - ctx.get(108)) / (1.0 + ((3.0 * constant) * host_pow((hk2 - 1.0), 2.0))))));
(htarg = hk2);
} else {
(htarg = (ctx.get(108) - ((0.15 * (ctx.get(229) - ctx.get(228))) / ctx.get(188))));
}
if (truth(wake)) {
(htarg = js_max(htarg, 1.01));
} else {
double hmax = (truth((ibl < ctx.arr<1>(16).at(is))) ? hlmax : htmax);
(htarg = js_max(htarg, hmax));
}
(ctx.arr<2>(12).at(is)[ibl] = htarg);
record_log(0, ibl, htarg, 0.0);
(direct = 0.0);
continue;
}
} else {
(ctx.arr<1>(45)[((4.0 * 6.0) + 1.0)] = 0.0);
(ctx.arr<1>(45)[((4.0 * 6.0) + 2.0)] = ctx.get(116));
(ctx.arr<1>(45)[((4.0 * 6.0) + 3.0)] = ctx.get(114));
(ctx.arr<1>(45)[((4.0 * 6.0) + 4.0)] = ctx.get(117));
(ctx.arr<1>(48)[4.0] = (htarg - ctx.get(113)));
fn_gauss1(4.0, ctx.arr<2>(44), ctx.arr<1>(48));
(dmax = js_max(js_max(__builtin_fabs((ctx.arr<1>(48)[2.0] / thi)), __builtin_fabs((ctx.arr<1>(48)[3.0] / dsi))), __builtin_fabs((ctx.arr<1>(48)[4.0] / uei))));
if (truth((ibl >= ctx.arr<1>(16).at(is)))) (dmax = js_max(dmax, __builtin_fabs((ctx.arr<1>(48)[1.0] / cti))));
double rlx = 1.0;
if (truth((dmax > 0.3))) (rlx = (0.3 / dmax));
if (truth((ibl >= ctx.arr<1>(16).at(is)))) (cti = (cti + (rlx * ctx.arr<1>(48)[1.0])));
(thi = (thi + (rlx * ctx.arr<1>(48)[2.0])));
(dsi = (dsi + (rlx * ctx.arr<1>(48)[3.0])));
(uei = (uei + (rlx * ctx.arr<1>(48)[4.0])));
}
if (truth((ibl >= ctx.arr<1>(16).at(is)))) {
(cti = js_min(cti, 0.30));
(cti = js_max(cti, 0.0000001));
}
if (truth((ibl <= ctx.arr<1>(13).at(is)))) {
(hklim = 1.02);
} else {
(hklim = 1.00005);
}
(msq = (((uei * uei) * ctx.get(130)) / (ctx.get(90) * (1.0 - (((0.5 * uei) * uei) * ctx.get(130))))));
double dsw = (dsi - dswaki);
(dsw = fn_dslim(ctx, dsw, thi, uei, msq, hklim));
(dsi = (dsw + dswaki));
if (truth((dmax <= 1.0e-5))) {
break;
}
}
if (truth((dmax > 1.0e-5))) {
record_log(1, ibl, double(is), dmax);
if ((truth((dmax > 0.1)) && truth((ibl > 3.0)))) {
if (truth((ibl <= ctx.arr<1>(13).at(is)))) {
(thi = (ctx.arr<2>(27)[ibm].at(is) * host_pow((ctx.arr<2>(60)[ibl].at(is) / ctx.arr<2>(60)[ibm].at(is)), 0.5)));
(dsi = (ctx.arr<2>(11)[ibm].at(is) * host_pow((ctx.arr<2>(60)[ibl].at(is) / ctx.arr<2>(60)[ibm].at(is)), 0.5)));
} else if (truth((ibl == (ctx.arr<1>(13).at(is) + 1.0)))) {
(cti = ctx.get(49));
(thi = ctx.get(198));
(dsi = ctx.get(79));
} else {
(thi = ctx.arr<2>(27)[ibm].at(is));
double ratlen = ((ctx.arr<2>(60)[ibl].at(is) - ctx.arr<2>(60)[ibm].at(is)) / (10.0 * ctx.arr<2>(11)[ibm].at(is)));
(dsi = ((ctx.arr<2>(11)[ibm].at(is) + (thi * ratlen)) / (1.0 + ratlen)));
}
if (truth((ibl == ctx.arr<1>(16).at(is)))) (cti = 0.05);
if (truth((ibl > ctx.arr<1>(16).at(is)))) (cti = ctx.arr<2>(5)[ibm].at(is));
(uei = ctx.arr<2>(30)[ibl].at(is));
if ((truth((ibl > 2.0)) && truth((ibl < ctx.arr<1>(22).at(is))))) {
(uei = (0.5 * (ctx.arr<2>(30)[(ibl - 1.0)].at(is) + ctx.arr<2>(30)[(ibl + 1.0)].at(is))));
}
}
fn_prepareBlStation(ctx, xsi, ami, cti, thi, dsi, dswaki, uei);
if ((truth((!truth(simi))) && truth((!truth(ctx.get(199)))))) {
fn_trchek(ctx);
(ami = ctx.get(2));
if (truth(ctx.get(195))) (ctx.arr<1>(16).at(is) = ibl);
if (truth((!truth(ctx.get(195))))) (ctx.arr<1>(16).at(is) = (ibl + 2.0));
}
if (truth((ibl < ctx.arr<1>(16).at(is)))) fn_blvar(1.0, ctx);
if (truth((ibl >= ctx.arr<1>(16).at(is)))) fn_blvar(2.0, ctx);
if (truth(wake)) fn_blvar(3.0, ctx);
if (truth((ibl < ctx.arr<1>(16).at(is)))) fn_blmid(1.0, ctx);
if (truth((ibl >= ctx.arr<1>(16).at(is)))) fn_blmid(2.0, ctx);
if (truth(wake)) fn_blmid(3.0, ctx);
}
if (truth((ibl < ctx.arr<1>(16).at(is)))) (ctx.arr<2>(5)[ibl].at(is) = ami);
if (truth((ibl >= ctx.arr<1>(16).at(is)))) (ctx.arr<2>(5)[ibl].at(is) = cti);
(ctx.arr<2>(27)[ibl].at(is) = thi);
(ctx.arr<2>(11)[ibl].at(is) = dsi);
(ctx.arr<2>(30)[ibl].at(is) = uei);
(ctx.arr<2>(20)[ibl].at(is) = (dsi * uei));
(ctx.arr<2>(25)[ibl].at(is) = ((((0.5 * ctx.get(158)) * ctx.get(203)) * ctx.get(203)) * ctx.get(11)));
(ctx.arr<2>(9)[ibl].at(is) = ((((((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(203)) * ctx.get(71)) * ctx.get(124)) * 0.5));
(ctx.arr<2>(6)[ibl].at(is) = ctx.get(42));
(ctx.arr<2>(7)[ibl].at(is) = ctx.get(59));
(ctx.arr<2>(29)[ibl].at(is) = (ctx.get(124) * ctx.get(189)));
(ctx.arr<2>(10)[ibl].at(is) = (ctx.get(103) * ctx.get(189)));
(ctx.arr<2>(21)[ibl].at(is) = ((ctx.get(158) * ctx.get(203)) * ctx.get(53)));
(ctx.arr<2>(23)[ibl].at(is) = (((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(189)));
(ctx.arr<2>(19)[ibl].at(is) = ((((((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(203)) * ctx.get(189)) * ctx.get(124)) * 0.5));
fn_prepareBlStation(ctx, xsi, ami, cti, thi, dsi, dswaki, uei);
for (u32 icom = 1u; truth((double(icom) <= ctx.get(150))); advance_index(icom)) {
(ctx.arr<1>(3).at(icom) = ctx.arr<1>(4).at(icom));
}
if (truth(true)) {
fn_syncComToVars(ctx, 1.0);
}
if ((truth(ctx.get(195)) || truth((ibl == ctx.arr<1>(13).at(is))))) {
ctx.put(199, 1.0);
(ctx.arr<1>(26).at(is) = ctx.get(196));
(ctx.arr<1>(61).at(is) = ctx.get(232));
}
ctx.put(195, 0.0);
if (truth((ibl == ctx.arr<1>(13).at(is)))) {
(thi = (ctx.arr<2>(27)[ctx.arr<1>(13)[1.0]][1.0] + ctx.arr<2>(27)[ctx.arr<1>(13)[2.0]][2.0]));
(dsi = ((ctx.arr<2>(11)[ctx.arr<1>(13)[1.0]][1.0] + ctx.arr<2>(11)[ctx.arr<1>(13)[2.0]][2.0]) + ctx.get(3)));
}
}
}
}

// js/xbl.js:436 (mrchdu)
void fn_mrchdu(Context ctx) {
((void)0);
double deps = 5.0e-6;
double senswt = 1000.0;
V1 NBL = ctx.arr<1>(22);
V1 IBLTE = ctx.arr<1>(13);
V1 ITRAN = ctx.arr<1>(16);
V2 XSSI = ctx.arr<2>(60);
V2 UEDG = ctx.arr<2>(30);
V2 THET = ctx.arr<2>(27);
V2 DSTR = ctx.arr<2>(11);
V2 CTAU = ctx.arr<2>(5);
V2 MASS = ctx.arr<2>(20);
V2 TAU = ctx.arr<2>(25);
V2 DIS = ctx.arr<2>(9);
V2 CTQ = ctx.arr<2>(6);
V2 DELT = ctx.arr<2>(7);
V2 TSTR = ctx.arr<2>(29);
V1 WGAP = ctx.arr<1>(56);
V1 XSSITR = ctx.arr<1>(61);
V1 TFORCE = ctx.arr<1>(26);
V1 ACRIT = ctx.arr<1>(0);
V2 VS2 = ctx.arr<2>(44);
V1 VS2F = ctx.arr<1>(45);
V1 VSREZ = ctx.arr<1>(48);
auto blprv = fn_blprv;
auto blkin = fn_blkin;
double combinePredictor = ([&]() -> double { double value = (blprv == fn_blprv); return truth(value) ? (blkin == fn_blkin) : value; }());
auto trchek = fn_trchek;
auto blsys = fn_blsys;
auto tesys = fn_tesys;
auto hkin = fn_hkin;
auto blvar = fn_blvar;
auto blmid = fn_blmid;
auto syncComToVars = fn_syncComToVars;
double vsStride = 6.0;
V2 vtmp = ctx.arr<2>(68);
V1 vztmp = ctx.arr<1>(69);
if ((truth((!bool(vtmp))) || truth((vtmp.length < 5.0)))) {
(vtmp = make_matrix(5.0, 5.0));
(vztmp = make_vec<1>(5.0));
(ctx.arr<2>(68) = vtmp);
(ctx.arr<1>(69) = vztmp);
} else if ((truth((!bool(vztmp))) || truth((vztmp.length < 5.0)))) {
(vztmp = make_vec<1>(5.0));
(ctx.arr<1>(69) = vztmp);
}
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
ctx.put(0, ACRIT.at(is));
fn_xifset(ctx, double(is));
double ibl = 2.0;
double xsi = XSSI[ibl].at(is);
double uei = UEDG[ibl].at(is);
ctx.put(4, 1.0);
double itrold = ITRAN.at(is);
ctx.put(195, 0.0);
ctx.put(199, 0.0);
(ITRAN.at(is) = IBLTE.at(is));
double sens = 0.0;
double ami = 0.0;
for ((ibl = 2.0); truth((ibl <= NBL.at(is))); (ibl += 1.0)) {
double ibm = (ibl - 1.0);
double simi = (ibl == 2.0);
double wake = (ibl > IBLTE.at(is));
ctx.put(183, simi);
ctx.put(227, wake);
(xsi = XSSI[ibl].at(is));
(uei = UEDG[ibl].at(is));
double thi = THET[ibl].at(is);
double dsi = DSTR[ibl].at(is);
double cti = 0.03;
if (truth((ibl < itrold))) {
(ami = CTAU[ibl].at(is));
} else {
(cti = CTAU[ibl].at(is));
if (truth((cti <= 0.0))) (cti = 0.03);
}
double dswaki = 0.0;
if (truth(wake)) {
double iw = (ibl - IBLTE.at(is));
(dswaki = WGAP[iw]);
}
if (truth((ibl <= IBLTE.at(is)))) {
(dsi = (js_max((dsi - dswaki), (1.02000 * thi)) + dswaki));
}
if (truth((ibl > IBLTE.at(is)))) {
(dsi = (js_max((dsi - dswaki), (1.00005 * thi)) + dswaki));
}
double ueref = 0.0;
double hkref = 0.0;
double sennew = sens;
double dmax = 0.0;
for (u32 itbl = 1u; truth((double(itbl) <= 25.0)); advance_index(itbl)) {
if (truth(combinePredictor)) {
fn_blprvKin(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
} else {
fn_blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_blkin(ctx);
}
if ((truth((!truth(simi))) && truth((!truth(ctx.get(199)))))) {
fn_trchek(ctx);
(ami = ctx.get(2));
if (truth(ctx.get(195))) (ITRAN.at(is) = ibl);
if (truth((!truth(ctx.get(195))))) (ITRAN.at(is) = (ibl + 2.0));
}
if (truth((ibl == (IBLTE.at(is) + 1.0)))) {
double tte = (THET[IBLTE[1.0]][1.0] + THET[IBLTE[2.0]][2.0]);
double dte = ((DSTR[IBLTE[1.0]][1.0] + DSTR[IBLTE[2.0]][2.0]) + ctx.get(3));
double cte = (((CTAU[IBLTE[1.0]][1.0] * THET[IBLTE[1.0]][1.0]) + (CTAU[IBLTE[2.0]][2.0] * THET[IBLTE[2.0]][2.0])) / tte);
ctx.put(49, cte);
ctx.put(198, tte);
ctx.put(79, dte);
fn_tesys(cte, tte, dte, ctx);
} else {
fn_blsys(ctx);
}
if (truth((double(itbl) == 1.0))) {
(ueref = ctx.get(203));
(hkref = ctx.get(113));
if ((truth((ibl < ITRAN.at(is))) && truth((ibl >= itrold)))) {
double uem = UEDG[(ibl - 1.0)].at(is);
double dsm = DSTR[(ibl - 1.0)].at(is);
double thm = THET[(ibl - 1.0)].at(is);
double msq = (((uem * uem) * ctx.get(130)) / (ctx.get(90) * (1.0 - (((0.5 * uem) * uem) * ctx.get(130)))));
(hkref = fn_hkin((dsm / thm), msq).hk);
}
if (truth((ibl < itrold))) {
if (truth(ctx.get(195))) (CTAU[ibl].at(is) = 0.03);
if (truth(ctx.get(199))) (CTAU[ibl].at(is) = CTAU[(ibl - 1.0)].at(is));
if ((truth(ctx.get(195)) || truth(ctx.get(199)))) {
(cti = CTAU[ibl].at(is));
ctx.put(181, cti);
}
}
}
if ((truth(simi) || truth((ibl == (IBLTE.at(is) + 1.0))))) {
(VS2F[((4.0 * vsStride) + 1.0)] = 0.0);
(VS2F[((4.0 * vsStride) + 2.0)] = 0.0);
(VS2F[((4.0 * vsStride) + 3.0)] = 0.0);
(VS2F[((4.0 * vsStride) + 4.0)] = ctx.get(205));
(VSREZ[4.0] = (ueref - ctx.get(203)));
} else {
for (u32 k = 1u; truth((double(k) <= 4.0)); advance_index(k)) {
(vztmp.at(k) = VSREZ.at(k));
for (u32 l = 1u; truth((double(l) <= 4.0)); advance_index(l)) {
(vtmp.at(k).at(l) = VS2F[((double(k) * vsStride) + double(l))]);
}
}
(vtmp[4.0][1.0] = 0.0);
(vtmp[4.0][2.0] = ctx.get(116));
(vtmp[4.0][3.0] = ctx.get(114));
(vtmp[4.0][4.0] = (ctx.get(117) * ctx.get(205)));
(vztmp[4.0] = 1.0);
fn_gauss1(4.0, vtmp, vztmp);
(sennew = (((senswt * vztmp[4.0]) * hkref) / ueref));
if (truth((double(itbl) <= 5.0))) {
(sens = sennew);
} else if (truth((double(itbl) <= 15.0))) {
(sens = (0.5 * (sens + sennew)));
}
(VS2F[((4.0 * vsStride) + 1.0)] = 0.0);
(VS2F[((4.0 * vsStride) + 2.0)] = (ctx.get(116) * hkref));
(VS2F[((4.0 * vsStride) + 3.0)] = (ctx.get(114) * hkref));
(VS2F[((4.0 * vsStride) + 4.0)] = (((ctx.get(117) * hkref) + (sens / ueref)) * ctx.get(205)));
(VSREZ[4.0] = (((-host_pow(hkref, 2.0)) * ((ctx.get(113) / hkref) - 1.0)) - (sens * ((ctx.get(203) / ueref) - 1.0))));
}
fn_gauss1(4.0, VS2, VSREZ);
(dmax = js_max(js_max(__builtin_fabs((VSREZ[2.0] / thi)), __builtin_fabs((VSREZ[3.0] / dsi))), __builtin_fabs((VSREZ[4.0] / uei))));
if (truth((ibl >= ITRAN.at(is)))) (dmax = js_max(dmax, __builtin_fabs((VSREZ[1.0] / (10.0 * cti)))));
double rlx = 1.0;
if (truth((dmax > 0.3))) (rlx = (0.3 / dmax));
if (truth((ibl < ITRAN.at(is)))) (ami = (ami + (rlx * VSREZ[1.0])));
if (truth((ibl >= ITRAN.at(is)))) (cti = (cti + (rlx * VSREZ[1.0])));
(thi = (thi + (rlx * VSREZ[2.0])));
(dsi = (dsi + (rlx * VSREZ[3.0])));
(uei = (uei + (rlx * VSREZ[4.0])));
if (truth((ibl >= ITRAN.at(is)))) {
(cti = js_min(cti, 0.30));
(cti = js_max(cti, 0.0000001));
}
double hklim = 0.0;
if (truth((ibl <= IBLTE.at(is)))) {
(hklim = 1.02);
} else {
(hklim = 1.00005);
}
double msq = (((uei * uei) * ctx.get(130)) / (ctx.get(90) * (1.0 - (((0.5 * uei) * uei) * ctx.get(130)))));
double dsw = (dsi - dswaki);
(dsw = fn_dslim(ctx, dsw, thi, uei, msq, hklim));
(dsi = (dsw + dswaki));
if (truth((dmax <= deps))) {
break;
}
}
if (truth((dmax > deps))) {
record_log(2, ibl, double(is), dmax);
}
if (((truth((dmax > deps)) && truth((dmax > 0.1))) && truth((ibl > 3.0)))) {
if (truth((ibl <= IBLTE.at(is)))) {
(thi = (THET[ibm].at(is) * host_pow((XSSI[ibl].at(is) / XSSI[ibm].at(is)), 0.5)));
(dsi = (DSTR[ibm].at(is) * host_pow((XSSI[ibl].at(is) / XSSI[ibm].at(is)), 0.5)));
(uei = UEDG[ibm].at(is));
} else if (truth((ibl == (IBLTE.at(is) + 1.0)))) {
(cti = ctx.get(49));
(thi = ctx.get(198));
(dsi = ctx.get(79));
(uei = UEDG[ibm].at(is));
} else {
(thi = THET[ibm].at(is));
double ratlen = ((XSSI[ibl].at(is) - XSSI[ibm].at(is)) / (10.0 * DSTR[ibm].at(is)));
(dsi = ((DSTR[ibm].at(is) + (thi * ratlen)) / (1.0 + ratlen)));
(uei = UEDG[ibm].at(is));
}
if (truth((ibl == ITRAN.at(is)))) (cti = 0.05);
if (truth((ibl > ITRAN.at(is)))) (cti = CTAU[ibm].at(is));
}
if (truth((dmax > deps))) {
if (truth(combinePredictor)) {
fn_blprvKin(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
} else {
fn_blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_blkin(ctx);
}
if ((truth((!truth(simi))) && truth((!truth(ctx.get(199)))))) {
fn_trchek(ctx);
(ami = ctx.get(2));
if (truth(ctx.get(195))) (ITRAN.at(is) = ibl);
if (truth((!truth(ctx.get(195))))) (ITRAN.at(is) = (ibl + 2.0));
}
if (truth((ibl < ITRAN.at(is)))) fn_blvar(1.0, ctx);
if (truth((ibl >= ITRAN.at(is)))) fn_blvar(2.0, ctx);
if (truth(wake)) fn_blvar(3.0, ctx);
if (truth((ibl < ITRAN.at(is)))) fn_blmid(1.0, ctx);
if (truth((ibl >= ITRAN.at(is)))) fn_blmid(2.0, ctx);
if (truth(wake)) fn_blmid(3.0, ctx);
}
if (truth((ibl < ITRAN.at(is)))) (CTAU[ibl].at(is) = ami);
if (truth((ibl >= ITRAN.at(is)))) (CTAU[ibl].at(is) = cti);
(THET[ibl].at(is) = thi);
(DSTR[ibl].at(is) = dsi);
(UEDG[ibl].at(is) = uei);
(MASS[ibl].at(is) = (dsi * uei));
(TAU[ibl].at(is) = ((((0.5 * ctx.get(158)) * ctx.get(203)) * ctx.get(203)) * ctx.get(11)));
(DIS[ibl].at(is) = ((((((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(203)) * ctx.get(71)) * ctx.get(124)) * 0.5));
(CTQ[ibl].at(is) = ctx.get(42));
(DELT[ibl].at(is) = ctx.get(59));
(TSTR[ibl].at(is) = (ctx.get(124) * ctx.get(189)));
(ctx.arr<2>(10)[ibl].at(is) = (ctx.get(103) * ctx.get(189)));
(ctx.arr<2>(21)[ibl].at(is) = ((ctx.get(158) * ctx.get(203)) * ctx.get(53)));
(ctx.arr<2>(23)[ibl].at(is) = (((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(189)));
(ctx.arr<2>(19)[ibl].at(is) = ((((((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(203)) * ctx.get(189)) * ctx.get(124)) * 0.5));
if (truth(combinePredictor)) {
fn_blprvKin(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
} else {
fn_blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
fn_blkin(ctx);
}
for (u32 icom = 1u; truth((double(icom) <= ctx.get(150))); advance_index(icom)) {
(ctx.arr<1>(3).at(icom) = ctx.arr<1>(4).at(icom));
}
if (truth(true)) {
fn_syncComToVars(ctx, 1.0);
}
if ((truth(ctx.get(195)) || truth((ibl == IBLTE.at(is))))) {
ctx.put(199, 1.0);
(TFORCE.at(is) = ctx.get(196));
(XSSITR.at(is) = ctx.get(232));
}
ctx.put(195, 0.0);
}
}
}

// js/xbl.js:764 (assembleMassRows)
void fn_assembleMassRows(Context ctx, double iv, V1 d1M, V1 u1M, V1 d2M, V1 u2M, V1 ule1M, V1 ule2M, double xiUle1, double xiUle2) {
double nsys = ctx.get(151);
V1 vs1 = ctx.arr<1>(43);
V1 vs2 = ctx.arr<1>(45);
V2 vm1 = ctx.arr<3>(40)[1.0];
V2 vm2 = ctx.arr<3>(40)[2.0];
V2 vm3 = ctx.arr<3>(40)[3.0];
double r1d1 = vs1[9.0];
double r1u1 = vs1[10.0];
double r1d2 = vs2[9.0];
double r1u2 = vs2[10.0];
double r1x = ((vs1[11.0] + vs2[11.0]) + ctx.arr<1>(49)[1.0]);
double r2d1 = vs1[15.0];
double r2u1 = vs1[16.0];
double r2d2 = vs2[15.0];
double r2u2 = vs2[16.0];
double r2x = ((vs1[17.0] + vs2[17.0]) + ctx.arr<1>(49)[2.0]);
double r3d1 = vs1[21.0];
double r3u1 = vs1[22.0];
double r3d2 = vs2[21.0];
double r3u2 = vs2[22.0];
double r3x = ((vs1[23.0] + vs2[23.0]) + ctx.arr<1>(49)[3.0]);
V1 flat = ctx.arr<1>(41);
double stride = (nsys + 1.0);
double plane = (stride * stride);
if ((bool(flat) && truth((flat.length == (4.0 * plane))))) {
double index = (stride + iv);
for (double jv = 1.0; truth((jv <= nsys)); ((jv += 1.0), (index += stride))) {
double d1 = d1M[jv];
double u1 = u1M[jv];
double d2 = d2M[jv];
double u2 = u2M[jv];
double xi = ((xiUle1 * ule1M[jv]) + (xiUle2 * ule2M[jv]));
(flat[(plane + index)] = (((((r1d1 * d1) + (r1u1 * u1)) + (r1d2 * d2)) + (r1u2 * u2)) + (r1x * xi)));
(flat[((2.0 * plane) + index)] = (((((r2d1 * d1) + (r2u1 * u1)) + (r2d2 * d2)) + (r2u2 * u2)) + (r2x * xi)));
(flat[((3.0 * plane) + index)] = (((((r3d1 * d1) + (r3u1 * u1)) + (r3d2 * d2)) + (r3u2 * u2)) + (r3x * xi)));
}
return;
}
for (u32 jv = 1u; truth((double(jv) <= nsys)); advance_index(jv)) {
double d1 = d1M.at(jv);
double u1 = u1M.at(jv);
double d2 = d2M.at(jv);
double u2 = u2M.at(jv);
double xi = ((xiUle1 * ule1M.at(jv)) + (xiUle2 * ule2M.at(jv)));
(vm1.at(jv)[iv] = (((((r1d1 * d1) + (r1u1 * u1)) + (r1d2 * d2)) + (r1u2 * u2)) + (r1x * xi)));
(vm2.at(jv)[iv] = (((((r2d1 * d1) + (r2u1 * u1)) + (r2d2 * d2)) + (r2u2 * u2)) + (r2x * xi)));
(vm3.at(jv)[iv] = (((((r3d1 * d1) + (r3u1 * u1)) + (r3d2 * d2)) + (r3u2 * u2)) + (r3x * xi)));
}
}

// js/xbl.js:817 (setbl)
void fn_setbl(Context ctx) {
if ((((truth((!bool(ctx.arr<3>(34)))) || truth((!bool(ctx.arr<3>(36))))) || truth((!bool(ctx.arr<3>(38))))) || truth((!bool(ctx.arr<3>(40)))))) {
if ((truth((!truth(finite(ctx.get(151))))) || truth((ctx.get(151) <= 0.0)))) {
fn_iblsys(ctx);
}
if (truth(true)) {
unsupported();
}
}
double nsys = ctx.get(151);
V1 COM1 = ctx.arr<1>(3);
V1 COM2 = ctx.arr<1>(4);
V3 VA = ctx.arr<3>(34);
V3 VB = ctx.arr<3>(36);
V3 VDEL = ctx.arr<3>(38);
V3 VM = ctx.arr<3>(40);
V2 VS1 = ctx.arr<2>(42);
V2 VS2 = ctx.arr<2>(44);
V1 VS1F = ctx.arr<1>(43);
V1 VS2F = ctx.arr<1>(45);
V1 VSX = ctx.arr<1>(49);
V1 VSREZ = ctx.arr<1>(48);
V1 VSR = ctx.arr<1>(47);
V1 VSM = ctx.arr<1>(46);
V2 VZ = ctx.arr<2>(51);
V1 NBL = ctx.arr<1>(22);
V2 ISYS = ctx.arr<2>(15);
V2 IPAN = ctx.arr<2>(14);
V2 VTI = ctx.arr<2>(50);
V2 DIJ = ctx.arr<2>(8);
V2 UEDG = ctx.arr<2>(30);
V2 UINV = ctx.arr<2>(31);
V2 UINV_A = ctx.arr<2>(32);
V2 CTAU = ctx.arr<2>(5);
V2 THET = ctx.arr<2>(27);
V2 DSTR = ctx.arr<2>(11);
V2 MASS = ctx.arr<2>(20);
V1 WGAP = ctx.arr<1>(56);
V2 XSSI = ctx.arr<2>(60);
V1 ITRAN = ctx.arr<1>(16);
V1 IBLTE = ctx.arr<1>(13);
double vsStride = 6.0;
V1 vs1 = VS1F;
V1 vs2 = VS2F;
double vs1r1 = (1.0 * vsStride);
double vs1r2 = (2.0 * vsStride);
double vs1r3 = (3.0 * vsStride);
double vs2r1 = (1.0 * vsStride);
double vs2r2 = (2.0 * vsStride);
double vs2r3 = (3.0 * vsStride);
double qinf = ctx.get(152);
double gamma = ctx.get(87);
double gamm1 = ctx.get(86);
double hvrat = ctx.get(132);
COM1.fill(0.0, 1.0);
COM2.fill(0.0, 1.0);
if ((truth((!bool(ctx.arr<1>(18)))) || truth((ctx.arr<1>(18).length < (nsys + 1.0))))) {
(ctx.arr<1>(18) = make_vec<1>((nsys + 1.0)));
(ctx.arr<1>(17) = make_vec<1>((nsys + 1.0)));
} else {
ctx.arr<1>(18).fill(0.0);
ctx.arr<1>(17).fill(0.0);
}
V1 ivToIs = ctx.arr<1>(18);
V1 ivToIbl = ctx.arr<1>(17);
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
for (u32 ibl = 2u; truth((double(ibl) <= NBL.at(is))); advance_index(ibl)) {
double iv = ISYS.at(ibl).at(is);
(ivToIs[iv] = double(is));
(ivToIbl[iv] = double(ibl));
}
}
for (u32 k = 1u; truth((double(k) <= 3.0)); advance_index(k)) {
V2 VAk = VA.at(k);
V2 VBk = VB.at(k);
V2 VDELk = VDEL.at(k);
V2 VMk = VM.at(k);
for (u32 j = 1u; truth((double(j) <= 2.0)); advance_index(j)) {
VAk.at(j).fill(0.0, 1.0, (nsys + 1.0));
VBk.at(j).fill(0.0, 1.0, (nsys + 1.0));
VDELk.at(j).fill(0.0, 1.0, (nsys + 1.0));
}
for (u32 j = 1u; truth((double(j) <= nsys)); advance_index(j)) {
VMk.at(j).fill(0.0, 1.0, (nsys + 1.0));
}
}
double clmr = (truth(ctx.get(137)) ? ctx.get(33) : ctx.get(34));
auto record_28257 = fn_mrcl(ctx, clmr);
double maClmr = record_28257.mCls;
double reClmr = record_28257.rCls;
double minf = ctx.get(146);
double reinf = ctx.get(161);
double msqClmr = ((2.0 * ctx.get(146)) * maClmr);
fn_comset(ctx);
ctx.put(85, gamma);
ctx.put(90, gamm1);
ctx.put(153, qinf);
ctx.put(190, ctx.get(192));
ctx.put(191, ctx.get(194));
ctx.put(168, host_pow((1.0 + ((0.5 * ctx.get(90)) * host_pow(minf, 2.0))), (1.0 / ctx.get(90))));
ctx.put(169, ((0.5 * ctx.get(168)) / (1.0 + ((0.5 * ctx.get(90)) * host_pow(minf, 2.0)))));
ctx.put(130, ((ctx.get(90) * host_pow((minf / ctx.get(153)), 2.0)) / (1.0 + ((0.5 * ctx.get(90)) * host_pow(minf, 2.0)))));
ctx.put(131, (((ctx.get(90) * host_pow((1.0 / ctx.get(153)), 2.0)) / (1.0 + ((0.5 * ctx.get(90)) * host_pow(minf, 2.0)))) - (((0.5 * ctx.get(90)) * ctx.get(130)) / (1.0 + ((0.5 * ctx.get(90)) * host_pow(minf, 2.0))))));
double herat = (1.0 - ((0.5 * host_pow(ctx.get(153), 2.0)) * ctx.get(130)));
double heratMs = (((-0.5) * host_pow(ctx.get(153), 2.0)) * ctx.get(131));
ctx.put(165, (((reinf * __builtin_sqrt(host_pow(herat, 3.0))) * (1.0 + hvrat)) / (herat + hvrat)));
ctx.put(167, ((__builtin_sqrt(host_pow(herat, 3.0)) * (1.0 + hvrat)) / (herat + hvrat)));
ctx.put(166, ((ctx.get(165) * ((1.5 / herat) - (1.0 / (herat + hvrat)))) * heratMs));
ctx.put(134, (ctx.has(133) ? ctx.get(133) : (ctx.has(134) ? ctx.get(134) : 0.0)));
ctx.put(83, ctx.arr<1>(56).get_or(1.0, 0.0));
if (truth((!truth(ctx.get(138))))) {
fn_mrchue(ctx);
ctx.put(138, 1.0);
}
fn_mrchdu(ctx);
double maxNbl = js_max(NBL[1.0], NBL[2.0]);
V2 usav = ctx.arr<2>(86);
if ((truth((!bool(usav))) || truth((usav.length < (maxNbl + 1.0))))) {
(usav = make_matrix(maxNbl + 1, 2.0 + 1));
(ctx.arr<2>(86) = usav);
} else {
for (u32 i = 0u; truth((double(i) <= maxNbl)); advance_index(i)) {
if ((truth((!bool(usav.at(i)))) || truth((usav.at(i).length < 3.0)))) {
(usav.at(i) = make_vec<1>(3.0));
}
}
}
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
double nblIs = NBL.at(is);
for (u32 ibl = 2u; truth((double(ibl) <= nblIs)); advance_index(ibl)) {
(usav.at(ibl).at(is) = UEDG.at(ibl).at(is));
}
}
fn_ueset(ctx);
double badUeSet = 0.0;
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
double nblIs = NBL.at(is);
for (u32 ibl = 2u; truth((double(ibl) <= nblIs)); advance_index(ibl)) {
if (truth((!truth(finite(UEDG.at(ibl).at(is)))))) {
(badUeSet += 1.0);
}
}
}
if (truth((badUeSet == 0.0))) {
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
double nblIs = NBL.at(is);
for (u32 ibl = 2u; truth((double(ibl) <= nblIs)); advance_index(ibl)) {
double temp = usav.at(ibl).at(is);
(usav.at(ibl).at(is) = UEDG.at(ibl).at(is));
(UEDG.at(ibl).at(is) = temp);
}
}
} else {
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
double nblIs = NBL.at(is);
for (u32 ibl = 2u; truth((double(ibl) <= nblIs)); advance_index(ibl)) {
(UEDG.at(ibl).at(is) = usav.at(ibl).at(is));
}
}
}
double ile1 = IPAN[2.0][1.0];
double ile2 = IPAN[2.0][2.0];
double ite1 = IPAN[IBLTE[1.0]][1.0];
double ite2 = IPAN[IBLTE[2.0]][2.0];
double jvte1 = ISYS[IBLTE[1.0]][1.0];
double jvte2 = ISYS[IBLTE[2.0]][2.0];
double dule1 = (UEDG[2.0][1.0] - usav[2.0][1.0]);
double dule2 = (UEDG[2.0][2.0] - usav[2.0][2.0]);
V1 u1M = ctx.arr<1>(82);
V1 u2M = ctx.arr<1>(83);
V1 d1M = ctx.arr<1>(66);
V1 d2M = ctx.arr<1>(67);
V1 ule1M = ctx.arr<1>(84);
V1 ule2M = ctx.arr<1>(85);
V1 ute1M = ctx.arr<1>(87);
V1 ute2M = ctx.arr<1>(88);
if ((truth((!bool(u1M))) || truth((u1M.length < (nsys + 1.0))))) {
(u1M = make_vec<1>((nsys + 1.0)));
(u2M = make_vec<1>((nsys + 1.0)));
(d1M = make_vec<1>((nsys + 1.0)));
(d2M = make_vec<1>((nsys + 1.0)));
(ule1M = make_vec<1>((nsys + 1.0)));
(ule2M = make_vec<1>((nsys + 1.0)));
(ute1M = make_vec<1>((nsys + 1.0)));
(ute2M = make_vec<1>((nsys + 1.0)));
(ctx.arr<1>(82) = u1M);
(ctx.arr<1>(83) = u2M);
(ctx.arr<1>(66) = d1M);
(ctx.arr<1>(67) = d2M);
(ctx.arr<1>(84) = ule1M);
(ctx.arr<1>(85) = ule2M);
(ctx.arr<1>(87) = ute1M);
(ctx.arr<1>(88) = ute2M);
} else {
u1M.fill(0.0);
u2M.fill(0.0);
d1M.fill(0.0);
d2M.fill(0.0);
ule1M.fill(0.0);
ule2M.fill(0.0);
ute1M.fill(0.0);
ute2M.fill(0.0);
}
double sourceCount = ((NBL[1.0] + NBL[2.0]) - 2.0);
Sources sources = ctx.sources(0);
if ((truth((!bool(sources))) || truth((sources.arr(0).length < sourceCount)))) {
(sources = ([&]() { auto value = make_sources(); value.arr(0) = make_vec<1>(sourceCount); value.arr(1) = make_vec<1>(sourceCount); value.arr(2) = make_vec<1>(sourceCount); return value; }()));
(ctx.sources(0) = sources);
}
double source = 0.0;
for (u32 js = 1u; truth((double(js) <= 2.0)); advance_index(js)) {
for (u32 jbl = 2u; truth((double(jbl) <= NBL.at(js))); advance_index(jbl)) {
double j = IPAN.at(jbl).at(js);
double jv = ISYS.at(jbl).at(js);
(sources.arr(0)[source] = j);
(sources.arr(1)[source] = jv);
(sources.arr(2)[source] = VTI.at(jbl).at(js));
(source += 1.0);
(ule1M[jv] = (((-VTI[2.0][1.0]) * VTI.at(jbl).at(js)) * DIJ[ile1][j]));
(ule2M[jv] = (((-VTI[2.0][2.0]) * VTI.at(jbl).at(js)) * DIJ[ile2][j]));
(ute1M[jv] = (((-VTI[IBLTE[1.0]][1.0]) * VTI.at(jbl).at(js)) * DIJ[ite1][j]));
(ute2M[jv] = (((-VTI[IBLTE[2.0]][2.0]) * VTI.at(jbl).at(js)) * DIJ[ite2][j]));
}
}
double ule1A = UINV_A[2.0][1.0];
double ule2A = UINV_A[2.0][2.0];
(ctx.arr<1>(28)[1.0] = 0.0);
(ctx.arr<1>(28)[2.0] = 0.0);
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
for (u32 j = 0u; truth((double(j) < sourceCount)); advance_index(j)) {
double jv = sources.arr(1).at(j);
(u1M[jv] = 0.0);
(d1M[jv] = 0.0);
}
double u1A = 0.0;
double d1A = 0.0;
double due1 = 0.0;
double dds1 = 0.0;
ctx.put(4, 1.0);
ctx.put(0, ctx.arr<1>(0).at(is));
fn_xifset(ctx, double(is));
ctx.put(195, 0.0);
ctx.put(199, 0.0);
double ami = 0.0;
double cti = 0.0;
for (u32 ibl = 2u; truth((double(ibl) <= NBL.at(is))); advance_index(ibl)) {
double iv = ISYS.at(ibl).at(is);
ctx.put(183, (double(ibl) == 2.0));
ctx.put(227, (double(ibl) > IBLTE.at(is)));
ctx.put(195, (double(ibl) == ITRAN.at(is)));
ctx.put(199, (double(ibl) > ITRAN.at(is)));
double i = IPAN.at(ibl).at(is);
double xsi = XSSI.at(ibl).at(is);
if (truth((double(ibl) < ITRAN.at(is)))) (ami = CTAU.at(ibl).at(is));
if (truth((double(ibl) >= ITRAN.at(is)))) (cti = CTAU.at(ibl).at(is));
double uei = UEDG.at(ibl).at(is);
double thi = THET.at(ibl).at(is);
double mdi = MASS.at(ibl).at(is);
double dsi = (mdi / uei);
double dswaki = 0.0;
if (truth(ctx.get(227))) {
double iw = (double(ibl) - IBLTE.at(is));
(dswaki = WGAP[iw]);
}
double d2M2 = (1.0 / uei);
double d2U2 = ((-dsi) / uei);
V1 influence = DIJ[i];
double direction = (-VTI.at(ibl).at(is));
for (u32 j = 0u; truth((double(j) < sourceCount)); advance_index(j)) {
double jv = sources.arr(1).at(j);
(u2M[jv] = ((direction * sources.arr(2).at(j)) * influence[sources.arr(0).at(j)]));
(d2M[jv] = (d2U2 * u2M[jv]));
}
(d2M[iv] += d2M2);
double u2A = UINV_A.at(ibl).at(is);
double d2A = (d2U2 * u2A);
double due2 = (UEDG.at(ibl).at(is) - usav.at(ibl).at(is));
double dds2 = (d2U2 * due2);
fn_prepareBlStation(ctx, xsi, ami, cti, thi, dsi, dswaki, uei);
if (truth(ctx.get(195))) {
fn_trchek(ctx);
(ami = ctx.get(2));
}
if ((truth((double(ibl) == ITRAN.at(is))) && truth((!truth(ctx.get(195)))))) {

}
double xiUle1 = (truth((double(is) == 1.0)) ? ctx.get(186) : (-ctx.get(186)));
double xiUle2 = (truth((double(is) == 1.0)) ? (-ctx.get(187)) : ctx.get(187));
if (truth((double(ibl) == (IBLTE.at(is) + 1.0)))) {
double tte = (THET[IBLTE[1.0]][1.0] + THET[IBLTE[2.0]][2.0]);
double dte = ((DSTR[IBLTE[1.0]][1.0] + DSTR[IBLTE[2.0]][2.0]) + ctx.get(3));
double cte = (((CTAU[IBLTE[1.0]][1.0] * THET[IBLTE[1.0]][1.0]) + (CTAU[IBLTE[2.0]][2.0] * THET[IBLTE[2.0]][2.0])) / tte);
fn_tesys(cte, tte, dte, ctx);
double tteTte1 = 1.0;
double tteTte2 = 1.0;
double dteMte1 = (1.0 / UEDG[IBLTE[1.0]][1.0]);
double dteUte1 = ((-DSTR[IBLTE[1.0]][1.0]) / UEDG[IBLTE[1.0]][1.0]);
double dteMte2 = (1.0 / UEDG[IBLTE[2.0]][2.0]);
double dteUte2 = ((-DSTR[IBLTE[2.0]][2.0]) / UEDG[IBLTE[2.0]][2.0]);
double cteCte1 = (THET[IBLTE[1.0]][1.0] / tte);
double cteCte2 = (THET[IBLTE[2.0]][2.0] / tte);
double cteTte1 = ((CTAU[IBLTE[1.0]][1.0] - cte) / tte);
double cteTte2 = ((CTAU[IBLTE[2.0]][2.0] - cte) / tte);
for (u32 j = 0u; truth((double(j) < sourceCount)); advance_index(j)) {
double jv = sources.arr(1).at(j);
(d1M[jv] = ((dteUte1 * ute1M[jv]) + (dteUte2 * ute2M[jv])));
}
(d1M[jvte1] += dteMte1);
(d1M[jvte2] += dteMte2);
(due1 = 0.0);
(dds1 = ((dteUte1 * (UEDG[IBLTE[1.0]][1.0] - usav[IBLTE[1.0]][1.0])) + (dteUte2 * (UEDG[IBLTE[2.0]][2.0] - usav[IBLTE[2.0]][2.0]))));
fn_assembleMassRows(ctx, iv, d1M, u1M, d2M, u2M, ule1M, ule2M, xiUle1, xiUle2);
(VB[1.0][1.0][iv] = vs1[(vs1r1 + 1.0)]);
(VB[1.0][2.0][iv] = vs1[(vs1r1 + 2.0)]);
(VA[1.0][1.0][iv] = vs2[(vs2r1 + 1.0)]);
(VA[1.0][2.0][iv] = vs2[(vs2r1 + 2.0)]);
(VDEL[1.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[1.0] * reClmr) + (VSM[1.0] * msqClmr)) : (((vs1[(vs1r1 + 4.0)] * u1A) + (vs1[(vs1r1 + 3.0)] * d1A)) + ((vs2[(vs2r1 + 4.0)] * u2A) + (vs2[(vs2r1 + 3.0)] * d2A)))));
if (truth((!truth(ctx.get(137))))) {
(VDEL[1.0][2.0][iv] += (((vs1[(vs1r1 + 5.0)] + vs2[(vs2r1 + 5.0)]) + VSX[1.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))));
}
(VDEL[1.0][1.0][iv] = (((VSREZ[1.0] + ((vs1[(vs1r1 + 4.0)] * due1) + (vs1[(vs1r1 + 3.0)] * dds1))) + ((vs2[(vs2r1 + 4.0)] * due2) + (vs2[(vs2r1 + 3.0)] * dds2))) + (((vs1[(vs1r1 + 5.0)] + vs2[(vs2r1 + 5.0)]) + VSX[1.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
(VB[2.0][1.0][iv] = vs1[(vs1r2 + 1.0)]);
(VB[2.0][2.0][iv] = vs1[(vs1r2 + 2.0)]);
(VA[2.0][1.0][iv] = vs2[(vs2r2 + 1.0)]);
(VA[2.0][2.0][iv] = vs2[(vs2r2 + 2.0)]);
(VDEL[2.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[2.0] * reClmr) + (VSM[2.0] * msqClmr)) : (((vs1[(vs1r2 + 4.0)] * u1A) + (vs1[(vs1r2 + 3.0)] * d1A)) + ((vs2[(vs2r2 + 4.0)] * u2A) + (vs2[(vs2r2 + 3.0)] * d2A)))));
if (truth((!truth(ctx.get(137))))) {
(VDEL[2.0][2.0][iv] += (((vs1[(vs1r2 + 5.0)] + vs2[(vs2r2 + 5.0)]) + VSX[2.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))));
}
(VDEL[2.0][1.0][iv] = (((VSREZ[2.0] + ((vs1[(vs1r2 + 4.0)] * due1) + (vs1[(vs1r2 + 3.0)] * dds1))) + ((vs2[(vs2r2 + 4.0)] * due2) + (vs2[(vs2r2 + 3.0)] * dds2))) + (((vs1[(vs1r2 + 5.0)] + vs2[(vs2r2 + 5.0)]) + VSX[2.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
(VB[3.0][1.0][iv] = vs1[(vs1r3 + 1.0)]);
(VB[3.0][2.0][iv] = vs1[(vs1r3 + 2.0)]);
(VA[3.0][1.0][iv] = vs2[(vs2r3 + 1.0)]);
(VA[3.0][2.0][iv] = vs2[(vs2r3 + 2.0)]);
(VDEL[3.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[3.0] * reClmr) + (VSM[3.0] * msqClmr)) : (((vs1[(vs1r3 + 4.0)] * u1A) + (vs1[(vs1r3 + 3.0)] * d1A)) + ((vs2[(vs2r3 + 4.0)] * u2A) + (vs2[(vs2r3 + 3.0)] * d2A)))));
if (truth((!truth(ctx.get(137))))) {
(VDEL[3.0][2.0][iv] += (((vs1[(vs1r3 + 5.0)] + vs2[(vs2r3 + 5.0)]) + VSX[3.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))));
}
(VDEL[3.0][1.0][iv] = (((VSREZ[3.0] + ((vs1[(vs1r3 + 4.0)] * due1) + (vs1[(vs1r3 + 3.0)] * dds1))) + ((vs2[(vs2r3 + 4.0)] * due2) + (vs2[(vs2r3 + 3.0)] * dds2))) + (((vs1[(vs1r3 + 5.0)] + vs2[(vs2r3 + 5.0)]) + VSX[3.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
(VZ[1.0][1.0] = (vs1[(vs1r1 + 1.0)] * cteCte1));
(VZ[1.0][2.0] = ((vs1[(vs1r1 + 1.0)] * cteTte1) + (vs1[(vs1r1 + 2.0)] * tteTte1)));
(VB[1.0][1.0][iv] = (vs1[(vs1r1 + 1.0)] * cteCte2));
(VB[1.0][2.0][iv] = ((vs1[(vs1r1 + 1.0)] * cteTte2) + (vs1[(vs1r1 + 2.0)] * tteTte2)));
(VZ[2.0][1.0] = (vs1[(vs1r2 + 1.0)] * cteCte1));
(VZ[2.0][2.0] = ((vs1[(vs1r2 + 1.0)] * cteTte1) + (vs1[(vs1r2 + 2.0)] * tteTte1)));
(VB[2.0][1.0][iv] = (vs1[(vs1r2 + 1.0)] * cteCte2));
(VB[2.0][2.0][iv] = ((vs1[(vs1r2 + 1.0)] * cteTte2) + (vs1[(vs1r2 + 2.0)] * tteTte2)));
(VZ[3.0][1.0] = (vs1[(vs1r3 + 1.0)] * cteCte1));
(VZ[3.0][2.0] = ((vs1[(vs1r3 + 1.0)] * cteTte1) + (vs1[(vs1r3 + 2.0)] * tteTte1)));
(VB[3.0][1.0][iv] = (vs1[(vs1r3 + 1.0)] * cteCte2));
(VB[3.0][2.0][iv] = ((vs1[(vs1r3 + 1.0)] * cteTte2) + (vs1[(vs1r3 + 2.0)] * tteTte2)));
} else {
fn_blsys(ctx);
fn_assembleMassRows(ctx, iv, d1M, u1M, d2M, u2M, ule1M, ule2M, xiUle1, xiUle2);
(VB[1.0][1.0][iv] = vs1[(vs1r1 + 1.0)]);
(VB[1.0][2.0][iv] = vs1[(vs1r1 + 2.0)]);
(VA[1.0][1.0][iv] = vs2[(vs2r1 + 1.0)]);
(VA[1.0][2.0][iv] = vs2[(vs2r1 + 2.0)]);
(VDEL[1.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[1.0] * reClmr) + (VSM[1.0] * msqClmr)) : ((((vs1[(vs1r1 + 4.0)] * u1A) + (vs1[(vs1r1 + 3.0)] * d1A)) + ((vs2[(vs2r1 + 4.0)] * u2A) + (vs2[(vs2r1 + 3.0)] * d2A))) + (((vs1[(vs1r1 + 5.0)] + vs2[(vs2r1 + 5.0)]) + VSX[1.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))))));
(VDEL[1.0][1.0][iv] = (((VSREZ[1.0] + ((vs1[(vs1r1 + 4.0)] * due1) + (vs1[(vs1r1 + 3.0)] * dds1))) + ((vs2[(vs2r1 + 4.0)] * due2) + (vs2[(vs2r1 + 3.0)] * dds2))) + (((vs1[(vs1r1 + 5.0)] + vs2[(vs2r1 + 5.0)]) + VSX[1.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
(VB[2.0][1.0][iv] = vs1[(vs1r2 + 1.0)]);
(VB[2.0][2.0][iv] = vs1[(vs1r2 + 2.0)]);
(VA[2.0][1.0][iv] = vs2[(vs2r2 + 1.0)]);
(VA[2.0][2.0][iv] = vs2[(vs2r2 + 2.0)]);
(VDEL[2.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[2.0] * reClmr) + (VSM[2.0] * msqClmr)) : ((((vs1[(vs1r2 + 4.0)] * u1A) + (vs1[(vs1r2 + 3.0)] * d1A)) + ((vs2[(vs2r2 + 4.0)] * u2A) + (vs2[(vs2r2 + 3.0)] * d2A))) + (((vs1[(vs1r2 + 5.0)] + vs2[(vs2r2 + 5.0)]) + VSX[2.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))))));
(VDEL[2.0][1.0][iv] = (((VSREZ[2.0] + ((vs1[(vs1r2 + 4.0)] * due1) + (vs1[(vs1r2 + 3.0)] * dds1))) + ((vs2[(vs2r2 + 4.0)] * due2) + (vs2[(vs2r2 + 3.0)] * dds2))) + (((vs1[(vs1r2 + 5.0)] + vs2[(vs2r2 + 5.0)]) + VSX[2.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
(VB[3.0][1.0][iv] = vs1[(vs1r3 + 1.0)]);
(VB[3.0][2.0][iv] = vs1[(vs1r3 + 2.0)]);
(VA[3.0][1.0][iv] = vs2[(vs2r3 + 1.0)]);
(VA[3.0][2.0][iv] = vs2[(vs2r3 + 2.0)]);
(VDEL[3.0][2.0][iv] = (truth(ctx.get(137)) ? ((VSR[3.0] * reClmr) + (VSM[3.0] * msqClmr)) : ((((vs1[(vs1r3 + 4.0)] * u1A) + (vs1[(vs1r3 + 3.0)] * d1A)) + ((vs2[(vs2r3 + 4.0)] * u2A) + (vs2[(vs2r3 + 3.0)] * d2A))) + (((vs1[(vs1r3 + 5.0)] + vs2[(vs2r3 + 5.0)]) + VSX[3.0]) * ((xiUle1 * ule1A) + (xiUle2 * ule2A))))));
(VDEL[3.0][1.0][iv] = (((VSREZ[3.0] + ((vs1[(vs1r3 + 4.0)] * due1) + (vs1[(vs1r3 + 3.0)] * dds1))) + ((vs2[(vs2r3 + 4.0)] * due2) + (vs2[(vs2r3 + 3.0)] * dds2))) + (((vs1[(vs1r3 + 5.0)] + vs2[(vs2r3 + 5.0)]) + VSX[3.0]) * ((xiUle1 * dule1) + (xiUle2 * dule2)))));
}
(ctx.arr<2>(25).at(ibl).at(is) = ((((0.5 * ctx.get(158)) * ctx.get(203)) * ctx.get(203)) * ctx.get(11)));
(ctx.arr<2>(9).at(ibl).at(is) = ((((((ctx.get(158) * ctx.get(203)) * ctx.get(203)) * ctx.get(203)) * ctx.get(71)) * ctx.get(124)) * 0.5));
(ctx.arr<2>(6).at(ibl).at(is) = ctx.get(42));
(ctx.arr<2>(7).at(ibl).at(is) = ctx.get(59));
(ctx.arr<2>(33).at(ibl).at(is) = (1.60 / (1.0 + ctx.get(212))));
if (truth(ctx.get(195))) {
ctx.put(199, 1.0);
(ctx.arr<1>(16).at(is) = double(ibl));
(ctx.arr<1>(26).at(is) = ctx.get(196));
(ctx.arr<1>(61).at(is) = ctx.get(232));
double str = (truth((double(is) == 1.0)) ? (ctx.get(185) - ctx.get(232)) : (ctx.get(185) + ctx.get(232)));
double chx = (ctx.get(233) - ctx.get(231));
double chy = (ctx.get(248) - ctx.get(247));
double chsq = (host_pow(chx, 2.0) + host_pow(chy, 2.0));
V1 s = ctx.arr<1>(24).subarray(1.0, (ctx.get(149) + 1.0));
V1 x = ctx.arr<1>(57).subarray(1.0, (ctx.get(149) + 1.0));
V1 y = ctx.arr<1>(63).subarray(1.0, (ctx.get(149) + 1.0));
V1 xp = ctx.arr<1>(59).subarray(1.0, (ctx.get(149) + 1.0));
V1 yp = ctx.arr<1>(65).subarray(1.0, (ctx.get(149) + 1.0));
double xtr = fn_seval(str, x, xp, s, ctx.get(149));
double ytr = fn_seval(str, y, yp, s, ctx.get(149));
(ctx.arr<1>(58).at(is) = ((((xtr - ctx.get(231)) * chx) + ((ytr - ctx.get(247)) * chy)) / chsq));
(ctx.arr<1>(64).at(is) = ((((ytr - ctx.get(247)) * chx) - ((xtr - ctx.get(231)) * chy)) / chsq));
}
ctx.put(195, 0.0);
if (truth((double(ibl) == ctx.arr<1>(13).at(is)))) {
ctx.put(199, 1.0);
ctx.put(227, 1.0);
fn_blvar(3.0, ctx);
fn_blmid(3.0, ctx);
}
for (u32 j = 0u; truth((double(j) < sourceCount)); advance_index(j)) {
double jv = sources.arr(1).at(j);
(u1M[jv] = u2M[jv]);
(d1M[jv] = d2M[jv]);
}
(u1A = u2A);
(d1A = d2A);
(due1 = due2);
(dds1 = dds2);
if ((truth((double(ibl) == ctx.arr<1>(16).at(is))) && truth((ctx.get(229) > ctx.get(228))))) {
if (truth((double(is) == 1.0))) {
(ctx.arr<1>(28).at(is) = (((ctx.get(135) - ctx.arr<1>(16).at(is)) + 3.0) - ((ctx.get(232) - ctx.get(228)) / (ctx.get(229) - ctx.get(228)))));
} else {
(ctx.arr<1>(28).at(is) = (((ctx.get(135) + ctx.arr<1>(16).at(is)) - 2.0) + ((ctx.get(232) - ctx.get(228)) / (ctx.get(229) - ctx.get(228)))));
}
}
for (u32 icom = 1u; truth((double(icom) <= ctx.get(150))); advance_index(icom)) {
(COM1.at(icom) = COM2.at(icom));
}
if (truth(true)) {
fn_syncComToVars(ctx, 1.0);
}
}
}
}

// js/xbl.js:1378 (xifset)
void fn_xifset(Context ctx, double is) {
if (truth((ctx.arr<1>(62)[is] >= 1.0))) {
ctx.put(230, ctx.arr<2>(60)[ctx.arr<1>(13)[is]][is]);
return;
}
double chx = (ctx.get(233) - ctx.get(231));
double chy = (ctx.get(248) - ctx.get(247));
double chsq = (host_pow(chx, 2.0) + host_pow(chy, 2.0));
double dxNan = 0.0;
double agNan = 0.0;
for (u32 i = 1u; truth((double(i) <= ctx.get(149))); advance_index(i)) {
(ctx.arr<1>(52).at(i) = ((((ctx.arr<1>(57).at(i) - ctx.get(231)) * chx) + ((ctx.arr<1>(63).at(i) - ctx.get(247)) * chy)) / chsq));
(ctx.arr<1>(53).at(i) = ((((ctx.arr<1>(63).at(i) - ctx.get(247)) * chx) - ((ctx.arr<1>(57).at(i) - ctx.get(231)) * chy)) / chsq));
}
V1 s1 = ctx.arr<1>(24).subarray(1.0, (ctx.get(149) + 1.0));
V1 w1 = ctx.arr<1>(52).subarray(1.0, (ctx.get(149) + 1.0));
V1 w2 = ctx.arr<1>(53).subarray(1.0, (ctx.get(149) + 1.0));
V1 w3 = ctx.arr<1>(54).subarray(1.0, (ctx.get(149) + 1.0));
V1 w4 = ctx.arr<1>(55).subarray(1.0, (ctx.get(149) + 1.0));
fn_splind(w1, w3, s1, ctx.get(149), (-999.0), (-999.0));
fn_splind(w2, w4, s1, ctx.get(149), (-999.0), (-999.0));
if (truth((is == 1.0))) {
double str = (ctx.get(184) + ((ctx.arr<1>(24)[1.0] - ctx.get(184)) * ctx.arr<1>(62)[is]));
(str = fn_sinvrt(str, ctx.arr<1>(62)[is], w1, w3, s1, ctx.get(149)));
ctx.put(230, js_min((ctx.get(185) - str), ctx.arr<2>(60)[ctx.arr<1>(13)[is]][is]));
} else {
double str = (ctx.get(184) + ((ctx.arr<1>(24)[ctx.get(149)] - ctx.get(184)) * ctx.arr<1>(62)[is]));
(str = fn_sinvrt(str, ctx.arr<1>(62)[is], w1, w3, s1, ctx.get(149)));
ctx.put(230, js_min((str - ctx.get(185)), ctx.arr<2>(60)[ctx.arr<1>(13)[is]][is]));
}
if (truth((ctx.get(230) < 0.0))) {
ctx.put(230, ctx.arr<2>(60)[ctx.arr<1>(13)[is]][is]);
}
}

// js/xbl.js:1421 (velocitySources)
Sources fn_velocitySources(Context ctx, double corrections) {
double count = ((ctx.arr<1>(22)[1.0] + ctx.arr<1>(22)[2.0]) - 2.0);
Sources sources = ctx.sources(1);
if ((truth((!bool(sources))) || truth((sources.arr(0).length < count)))) {
(sources = ([&]() { auto value = make_sources(); value.arr(0) = make_vec<1>(count); value.arr(2) = make_vec<1>(count); value.arr(3) = make_vec<1>(count); value.arr(4) = make_vec<1>(count); return value; }()));
(ctx.sources(1) = sources);
}
(sources.count() = count);
double index = 0.0;
for (u32 side = 1u; truth((double(side) <= 2.0)); advance_index(side)) {
for (u32 ibl = 2u; truth((double(ibl) <= ctx.arr<1>(22).at(side))); advance_index(ibl)) {
(sources.arr(0)[index] = ctx.arr<2>(14).at(ibl).at(side));
(sources.arr(2)[index] = ctx.arr<2>(50).at(ibl).at(side));
if (truth(corrections)) {
double iv = ctx.arr<2>(15).at(ibl).at(side);
(sources.arr(3)[index] = (ctx.arr<2>(20).at(ibl).at(side) + ctx.arr<3>(38)[3.0][1.0][iv]));
(sources.arr(4)[index] = (-ctx.arr<3>(38)[3.0][2.0][iv]));
} else {
(sources.arr(3)[index] = ctx.arr<2>(20).at(ibl).at(side));
}
(index += 1.0);
}
}
return sources;
}

// js/xbl.js:1726 (ueset)
void fn_ueset(Context ctx) {
Sources sources = fn_velocitySources(ctx, 0.0);
for (u32 is = 1u; truth((double(is) <= 2.0)); advance_index(is)) {
for (u32 ibl = 2u; truth((double(ibl) <= ctx.arr<1>(22).at(is))); advance_index(ibl)) {
double i = ctx.arr<2>(14).at(ibl).at(is);
V1 influence = ctx.arr<2>(8)[i];
double direction = (-ctx.arr<2>(50).at(ibl).at(is));
double dui = 0.0;
for (u32 j = 0u; truth((double(j) < sources.count())); advance_index(j)) {
double ueM = ((direction * sources.arr(2).at(j)) * influence[sources.arr(0).at(j)]);
(dui += (ueM * sources.arr(3).at(j)));
}
(ctx.arr<2>(30).at(ibl).at(is) = (ctx.arr<2>(31).at(ibl).at(is) + dui));
}
}
}

// js/xsolve.js:92 (gauss1)
void fn_gauss1(double n, V2 z, V1 r) {
for (u32 np = 1u; truth((double(np) < n)); advance_index(np)) {
double np1 = (double(np) + 1.0);
double nx = double(np);
for (double nn = np1; truth((nn <= n)); (nn += 1.0)) {
if (truth((__builtin_fabs(z[nn].at(np)) > __builtin_fabs(z[nx].at(np))))) (nx = nn);
}
double pivot = (1.0 / z[nx].at(np));
(z[nx].at(np) = z.at(np).at(np));
for (double l = np1; truth((l <= n)); (l += 1.0)) {
double temp = (z[nx][l] * pivot);
(z[nx][l] = z.at(np)[l]);
(z.at(np)[l] = temp);
}
double temp = (r[nx] * pivot);
(r[nx] = r.at(np));
(r.at(np) = temp);
for (double k = np1; truth((k <= n)); (k += 1.0)) {
double ztmp = z[k].at(np);
for (double l = np1; truth((l <= n)); (l += 1.0)) {
(z[k][l] = (z[k][l] - (ztmp * z.at(np)[l])));
}
(r[k] = (r[k] - (ztmp * r.at(np))));
}
}
(r[n] = (r[n] / z[n][n]));
for (double np = (n - 1.0); truth((np >= 1.0)); (np -= 1.0)) {
double value = r[np];
for (double k = (np + 1.0); truth((k <= n)); (k += 1.0)) {
(value -= (z[np][k] * r[k]));
}
(r[np] = value);
}
}

// js/spline.js:21 (trisol)
void fn_trisol(V1 a, V1 b, V1 c, V1 d, double kk, double dOffset) {
for (u32 k = 1u; truth((double(k) < kk)); advance_index(k)) {
double km = (double(k) - 1.0);
(c[km] = (c[km] / a[km]));
double dk = (d[(dOffset + km)] / a[km]);
(d[(dOffset + km)] = dk);
(a.at(k) = (a.at(k) - (b.at(k) * c[km])));
(d[(dOffset + double(k))] = (d[(dOffset + double(k))] - (b.at(k) * dk)));
}
(d[((dOffset + kk) - 1.0)] = (d[((dOffset + kk) - 1.0)] / a[(kk - 1.0)]));
for (double k = (kk - 2.0); truth((k >= 0.0)); (k -= 1.0)) {
(d[(dOffset + k)] = (d[(dOffset + k)] - (c[k] * d[((dOffset + k) + 1.0)])));
}
}

// js/spline.js:68 (splindOffset)
void fn_splindOffset(V1 x, V1 xs, V1 s, double n, double xs1, double xs2, double offset) {
check_n(n);
V1 a = make_vec<1>(n);
V1 b = make_vec<1>(n);
V1 c = make_vec<1>(n);
for (u32 i = 1u; truth((double(i) <= (n - 2.0))); advance_index(i)) {
double im1 = ((offset + double(i)) - 1.0);
double i0 = (offset + double(i));
double ip1 = ((offset + double(i)) + 1.0);
double dsm = (s[i0] - s[im1]);
double dsp = (s[ip1] - s[i0]);
(b.at(i) = dsp);
(a.at(i) = (2.0 * (dsm + dsp)));
(c.at(i) = dsm);
(xs[(offset + double(i))] = (3.0 * ((((x[ip1] - x[i0]) * dsm) / dsp) + (((x[i0] - x[im1]) * dsp) / dsm))));
}
if (truth((xs1 == 999.0))) {
(a[0.0] = 2.0);
(c[0.0] = 1.0);
(xs[(offset + 0.0)] = ((3.0 * (x[(offset + 1.0)] - x[(offset + 0.0)])) / (s[(offset + 1.0)] - s[(offset + 0.0)])));
} else if (truth((xs1 == (-999.0)))) {
(a[0.0] = 1.0);
(c[0.0] = 1.0);
(xs[(offset + 0.0)] = ((2.0 * (x[(offset + 1.0)] - x[(offset + 0.0)])) / (s[(offset + 1.0)] - s[(offset + 0.0)])));
} else {
(a[0.0] = 1.0);
(c[0.0] = 0.0);
(xs[(offset + 0.0)] = xs1);
}
if (truth((xs2 == 999.0))) {
(b[(n - 1.0)] = 1.0);
(a[(n - 1.0)] = 2.0);
(xs[((offset + n) - 1.0)] = ((3.0 * (x[((offset + n) - 1.0)] - x[((offset + n) - 2.0)])) / (s[((offset + n) - 1.0)] - s[((offset + n) - 2.0)])));
} else if (truth((xs2 == (-999.0)))) {
(b[(n - 1.0)] = 1.0);
(a[(n - 1.0)] = 1.0);
(xs[((offset + n) - 1.0)] = ((2.0 * (x[((offset + n) - 1.0)] - x[((offset + n) - 2.0)])) / (s[((offset + n) - 1.0)] - s[((offset + n) - 2.0)])));
} else {
(a[(n - 1.0)] = 1.0);
(b[(n - 1.0)] = 0.0);
(xs[((offset + n) - 1.0)] = xs2);
}
if (((truth((n == 2.0)) && truth((xs1 == (-999.0)))) && truth((xs2 == (-999.0))))) {
(b[(n - 1.0)] = 1.0);
(a[(n - 1.0)] = 2.0);
(xs[((offset + n) - 1.0)] = ((3.0 * (x[((offset + n) - 1.0)] - x[((offset + n) - 2.0)])) / (s[((offset + n) - 1.0)] - s[((offset + n) - 2.0)])));
}
fn_trisol(a, b, c, xs, n, offset);
}

// js/spline.js:124 (splind)
void fn_splind(V1 x, V1 xs, V1 s, double n, double xs1, double xs2) {
fn_splindOffset(x, xs, s, n, xs1, xs2, 0.0);
}

// js/spline.js:155 (seval)
double fn_seval(double ss, V1 x, V1 xs, V1 s, double n) {
double ilow = 0.0;
double i = (n - 1.0);
while (truth(((i - ilow) > 1.0))) {
double imid = __builtin_floor(((i + ilow) / 2.0));
if (truth((ss < s[imid]))) {
(i = imid);
} else {
(ilow = imid);
}
}
double ds = (s[i] - s[(i - 1.0)]);
double t = ((ss - s[(i - 1.0)]) / ds);
double cx1 = (((ds * xs[(i - 1.0)]) - x[i]) + x[(i - 1.0)]);
double cx2 = (((ds * xs[i]) - x[i]) + x[(i - 1.0)]);
return (((t * x[i]) + ((1.0 - t) * x[(i - 1.0)])) + ((t - (t * t)) * (((1.0 - t) * cx1) - (t * cx2))));
}

// js/spline.js:176 (deval)
double fn_deval(double ss, V1 x, V1 xs, V1 s, double n) {
double ilow = 0.0;
double i = (n - 1.0);
while (truth(((i - ilow) > 1.0))) {
double imid = __builtin_floor(((i + ilow) / 2.0));
if (truth((ss < s[imid]))) {
(i = imid);
} else {
(ilow = imid);
}
}
double ds = (s[i] - s[(i - 1.0)]);
double t = ((ss - s[(i - 1.0)]) / ds);
double cx1 = (((ds * xs[(i - 1.0)]) - x[i]) + x[(i - 1.0)]);
double cx2 = (((ds * xs[i]) - x[i]) + x[(i - 1.0)]);
double devalValue = (((x[i] - x[(i - 1.0)]) + (((1.0 - (4.0 * t)) + ((3.0 * t) * t)) * cx1)) + ((t * ((3.0 * t) - 2.0)) * cx2));
(devalValue = (devalValue / ds));
return devalValue;
}

// js/spline.js:297 (sinvrt)
double fn_sinvrt(double si, double xi, V1 x, V1 xs, V1 s, double n) {
double sisav = si;
double siValue = si;
for (u32 iter = 0u; truth((double(iter) < 10.0)); advance_index(iter)) {
double res = (fn_seval(siValue, x, xs, s, n) - xi);
double resp = fn_deval(siValue, x, xs, s, n);
double ds = ((-res) / resp);
(siValue += ds);
if (truth((__builtin_fabs((ds / (s[(n - 1.0)] - s[0.0]))) < 1.0e-5))) {
return siValue;
}
}
return sisav;
}
enum : u32 {scalar_AMCRIT=0,scalar_AMPL1=1,scalar_AMPL2=2,scalar_ANTE=3,scalar_BULE=4,scalar_CF1=5,scalar_CF1_D1=6,scalar_CF1_MS=7,scalar_CF1_RE=8,scalar_CF1_T1=9,scalar_CF1_U1=10,scalar_CF2=11,scalar_CF2_D2=12,scalar_CF2_HK2=13,scalar_CF2_M2=14,scalar_CF2_MS=15,scalar_CF2_RE=16,scalar_CF2_RT2=17,scalar_CF2_T2=18,scalar_CF2_U2=19,scalar_CFFAC=20,scalar_CFM=21,scalar_CFM_D1=22,scalar_CFM_D2=23,scalar_CFM_HKA=24,scalar_CFM_MA=25,scalar_CFM_MS=26,scalar_CFM_RE=27,scalar_CFM_RTA=28,scalar_CFM_T1=29,scalar_CFM_T2=30,scalar_CFM_U1=31,scalar_CFM_U2=32,scalar_CL=33,scalar_CLSPEC=34,scalar_CPSTAR=35,scalar_CQ1=36,scalar_CQ1_D1=37,scalar_CQ1_MS=38,scalar_CQ1_RE=39,scalar_CQ1_T1=40,scalar_CQ1_U1=41,scalar_CQ2=42,scalar_CQ2_D2=43,scalar_CQ2_MS=44,scalar_CQ2_RE=45,scalar_CQ2_T2=46,scalar_CQ2_U2=47,scalar_CTCON=48,scalar_CTE=49,scalar_CTRCEX=50,scalar_CTRCON=51,scalar_D1=52,scalar_D2=53,scalar_DE1=54,scalar_DE1_D1=55,scalar_DE1_MS=56,scalar_DE1_T1=57,scalar_DE1_U1=58,scalar_DE2=59,scalar_DE2_D2=60,scalar_DE2_MS=61,scalar_DE2_T2=62,scalar_DE2_U2=63,scalar_DI1=64,scalar_DI1_D1=65,scalar_DI1_MS=66,scalar_DI1_RE=67,scalar_DI1_S1=68,scalar_DI1_T1=69,scalar_DI1_U1=70,scalar_DI2=71,scalar_DI2_D2=72,scalar_DI2_MS=73,scalar_DI2_RE=74,scalar_DI2_S2=75,scalar_DI2_T2=76,scalar_DI2_U2=77,scalar_DLCON=78,scalar_DTE=79,scalar_DUXCON=80,scalar_DW1=81,scalar_DW2=82,scalar_DWTE=83,scalar_GACON=84,scalar_GAMBL=85,scalar_GAMM1=86,scalar_GAMMA=87,scalar_GBCON=88,scalar_GCCON=89,scalar_GM1BL=90,scalar_H1=91,scalar_H1_D1=92,scalar_H1_T1=93,scalar_H2=94,scalar_H2_D2=95,scalar_H2_T2=96,scalar_H2_U2=97,scalar_HC1=98,scalar_HC1_D1=99,scalar_HC1_MS=100,scalar_HC1_T1=101,scalar_HC1_U1=102,scalar_HC2=103,scalar_HC2_D2=104,scalar_HC2_MS=105,scalar_HC2_T2=106,scalar_HC2_U2=107,scalar_HK1=108,scalar_HK1_D1=109,scalar_HK1_MS=110,scalar_HK1_T1=111,scalar_HK1_U1=112,scalar_HK2=113,scalar_HK2_D2=114,scalar_HK2_MS=115,scalar_HK2_T2=116,scalar_HK2_U2=117,scalar_HS1=118,scalar_HS1_D1=119,scalar_HS1_MS=120,scalar_HS1_RE=121,scalar_HS1_T1=122,scalar_HS1_U1=123,scalar_HS2=124,scalar_HS2_D2=125,scalar_HS2_MS=126,scalar_HS2_RE=127,scalar_HS2_T2=128,scalar_HS2_U2=129,scalar_HSTINV=130,scalar_HSTINV_MS=131,scalar_HVRAT=132,scalar_IDAMP=133,scalar_IDAMPV=134,scalar_IST=135,scalar_IVX=136,scalar_LALFA=137,scalar_LBLINI=138,scalar_M1=139,scalar_M1_MS=140,scalar_M1_U1=141,scalar_M2=142,scalar_M2_MS=143,scalar_M2_U2=144,scalar_MATYP=145,scalar_MINF=146,scalar_MINF1=147,scalar_MINF_CL=148,scalar_N=149,scalar_NCOM=150,scalar_NSYS=151,scalar_QINF=152,scalar_QINFBL=153,scalar_QSTAR=154,scalar_R1=155,scalar_R1_MS=156,scalar_R1_U1=157,scalar_R2=158,scalar_R2_MS=159,scalar_R2_U2=160,scalar_REINF=161,scalar_REINF1=162,scalar_REINF_CL=163,scalar_RETYP=164,scalar_REYBL=165,scalar_REYBL_MS=166,scalar_REYBL_RE=167,scalar_RSTBL=168,scalar_RSTBL_MS=169,scalar_RT1=170,scalar_RT1_MS=171,scalar_RT1_RE=172,scalar_RT1_T1=173,scalar_RT1_U1=174,scalar_RT2=175,scalar_RT2_MS=176,scalar_RT2_RE=177,scalar_RT2_T2=178,scalar_RT2_U2=179,scalar_S1=180,scalar_S2=181,scalar_SCCON=182,scalar_SIMI=183,scalar_SLE=184,scalar_SST=185,scalar_SST_GO=186,scalar_SST_GP=187,scalar_T1=188,scalar_T2=189,scalar_TKBL=190,scalar_TKBL_MS=191,scalar_TKLAM=192,scalar_TKL_MSQ=193,scalar_TK_MSQ=194,scalar_TRAN=195,scalar_TRFORC=196,scalar_TRFREE=197,scalar_TTE=198,scalar_TURB=199,scalar_U1=200,scalar_U1_MS=201,scalar_U1_UEI=202,scalar_U2=203,scalar_U2_MS=204,scalar_U2_UEI=205,scalar_US1=206,scalar_US1_D1=207,scalar_US1_MS=208,scalar_US1_RE=209,scalar_US1_T1=210,scalar_US1_U1=211,scalar_US2=212,scalar_US2_D2=213,scalar_US2_MS=214,scalar_US2_RE=215,scalar_US2_T2=216,scalar_US2_U2=217,scalar_V1=218,scalar_V1_MS=219,scalar_V1_RE=220,scalar_V1_U1=221,scalar_V2=222,scalar_V2_MS=223,scalar_V2_RE=224,scalar_V2_U2=225,scalar_VACCEL=226,scalar_WAKE=227,scalar_X1=228,scalar_X2=229,scalar_XIFORC=230,scalar_XLE=231,scalar_XT=232,scalar_XTE=233,scalar_XT_A1=234,scalar_XT_A2=235,scalar_XT_D1=236,scalar_XT_D2=237,scalar_XT_MS=238,scalar_XT_RE=239,scalar_XT_T1=240,scalar_XT_T2=241,scalar_XT_U1=242,scalar_XT_U2=243,scalar_XT_X1=244,scalar_XT_X2=245,scalar_XT_XF=246,scalar_YLE=247,scalar_YTE=248};
enum : u32 {array_ACRIT=0,array_C1SAV=1,array_C2SAV=2,array_COM1=3,array_COM2=4,array_CTAU=5,array_CTQ=6,array_DELT=7,array_DIJ=8,array_DIS=9,array_DSSR=10,array_DSTR=11,array_HTARG=12,array_IBLTE=13,array_IPAN=14,array_ISYS=15,array_ITRAN=16,array_IV_TO_IBL=17,array_IV_TO_IS=18,array_KDEF=19,array_MASS=20,array_MDEF=21,array_NBL=22,array_PDEF=23,array_S=24,array_TAU=25,array_TFORCE=26,array_THET=27,array_TINDEX=28,array_TSTR=29,array_UEDG=30,array_UINV=31,array_UINV_A=32,array_USLP=33,array_VA=34,array_VAF=35,array_VB=36,array_VBF=37,array_VDEL=38,array_VDELF=39,array_VM=40,array_VMF=41,array_VS1=42,array_VS1F=43,array_VS2=44,array_VS2F=45,array_VSM=46,array_VSR=47,array_VSREZ=48,array_VSX=49,array_VTI=50,array_VZ=51,array_W1=52,array_W2=53,array_W3=54,array_W4=55,array_WGAP=56,array_X=57,array_XOCTR=58,array_XP=59,array_XSSI=60,array_XSSITR=61,array_XSTRIP=62,array_Y=63,array_YOCTR=64,array_YP=65,array__D1M=66,array__D2M=67,array__MRCHDU_VTMP=68,array__MRCHDU_VZTMP=69,array__TRDIF_BL1=70,array__TRDIF_BL2=71,array__TRDIF_BLM=72,array__TRDIF_BLR=73,array__TRDIF_BLREZ=74,array__TRDIF_BLX=75,array__TRDIF_BT1=76,array__TRDIF_BT2=77,array__TRDIF_BTM=78,array__TRDIF_BTR=79,array__TRDIF_BTREZ=80,array__TRDIF_BTX=81,array__U1M=82,array__U2M=83,array__ULE1M=84,array__ULE2M=85,array__USAV=86,array__UTE1M=87,array__UTE2M=88};

// Share SETBL scratch storage with the existing strict-double BLSOLV kernel.
// Derived from XFOIL by Mark Drela. GNU General Public License v2 or later.
extern "C" int blsolv_flat(int, int, int, double, double, double,
                           double*, const double*, double*, double*, const double*);

V1 flat_tensor(V3 tensor, V1 flat, u32 middle, u32 stride) {
  if (tensor.length != 4 || flat.length != 4 * middle * stride) unsupported();
  for (u32 k = 0; k < 4; k++) {
    if (tensor[k].length != middle) unsupported();
    for (u32 j = 0; j < middle; j++) {
      auto row = tensor[k][j];
      if (row.length != stride || row.pointer != flat.pointer + 8 * (k * middle + j) * stride) unsupported();
    }
  }
  return flat;
}
bool overlaps(V1 a, V1 b) {
  return a.pointer < b.pointer + 8 * b.length && b.pointer < a.pointer + 8 * a.length;
}
void solve_boundary_layer_matrix(Context ctx) {
  u32 n = size_value(ctx.get(scalar_NSYS));
  if (n < 32 || n > 1024) unsupported();
  u32 stride = n + 1;
  V1 matrices[] = {
    flat_tensor(ctx.arr<3>(array_VA), ctx.arr<1>(array_VAF), 3, stride),
    flat_tensor(ctx.arr<3>(array_VB), ctx.arr<1>(array_VBF), 3, stride),
    flat_tensor(ctx.arr<3>(array_VM), ctx.arr<1>(array_VMF), stride, stride),
    flat_tensor(ctx.arr<3>(array_VDEL), ctx.arr<1>(array_VDELF), 3, stride),
  };
  for (u32 i = 0; i < 4; i++) for (u32 j = i + 1; j < 4; j++) {
    if (overlaps(matrices[i], matrices[j])) unsupported();
  }
  auto independent = [&](V1 a) {
    for (auto matrix : matrices) if (overlaps(a, matrix)) unsupported();
  };
  V1 s = ctx.arr<1>(array_S), iblte = ctx.arr<1>(array_IBLTE);
  independent(s); independent(iblte);
  u32 count = size_value(ctx.get(scalar_N));
  if (count < 1 || count >= s.length) unsupported();
  V1 edge = ctx.arr<2>(array_ISYS)[iblte[1]];
  independent(edge);
  u32 ivte1 = size_value(edge[1]);
  if (ivte1 < 1 || ivte1 > n) unsupported();
  u32 ivz = 0;
  auto vz = make_vec<1>(12);
  if (ivte1 < n) {
    V1 wake = ctx.arr<2>(array_ISYS)[iblte[2] + 1];
    independent(wake);
    ivz = size_value(wake[2]);
    if (ivz < 1 || ivz > n) unsupported();
    for (u32 k = 1; k <= 3; k++) {
      V1 row = ctx.arr<2>(array_VZ)[k];
      independent(row);
      vz[3 * k + 1] = row[1]; vz[3 * k + 2] = row[2];
    }
  }
  double vaccel = ctx.get(scalar_VACCEL);
  if (!finite(vaccel) || !finite(s[1]) || !finite(s[count])) unsupported();
  if (!blsolv_flat(n, ivte1, ivz, vaccel, s[1], s[count],
      (double*)matrices[0].pointer, (double*)matrices[1].pointer,
      (double*)matrices[2].pointer, (double*)matrices[3].pointer, (double*)vz.pointer)) unsupported();
}

extern "C" void boundary_layer(unsigned scalar, unsigned array, unsigned source, unsigned begin, unsigned end) {
 arena_begin(begin,end);
 fn_setbl(Context{(Slot*)scalar, (Descriptor*)array, (Sources*)source});
}
extern "C" void boundary_layer_step(unsigned scalar, unsigned array, unsigned source, unsigned begin, unsigned end) {
 arena_begin(begin,end);
 Context ctx{(Slot*)scalar, (Descriptor*)array, (Sources*)source};
 fn_setbl(ctx);
 solve_boundary_layer_matrix(ctx);
}
