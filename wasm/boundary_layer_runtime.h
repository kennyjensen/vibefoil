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
