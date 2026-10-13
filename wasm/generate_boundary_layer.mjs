// GPL-2.0-or-later. Generate a strict-double kernel from the maintained JS
// equations. This deliberately supports only the numeric subset used below;
// unsupported syntax is a build error, never an approximate translation.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { parse } from 'acorn';

const root = new URL('../', import.meta.url);
const selected = {
  'js/xblsys.js': ['syncComToVars', 'syncVarsToCom', 'copyCom', 'trchek', 'axset', 'trchek2', 'blsys', 'tesys', 'blprv', 'blprvKin', 'predictBlStation', 'blkin', 'blvar', 'blmid', 'trdif', 'bldif', 'dampl', 'dampl2', 'hkin', 'dil', 'dilw', 'hsl', 'cfl', 'dit', 'hst', 'cft', 'hct'],
  'js/xbl.js': ['prepareBlStation', 'dslim', 'mrcl', 'comset', 'iblsys', 'mrchue', 'mrchdu', 'assembleMassRows', 'setbl', 'xifset', 'velocitySources', 'ueset'],
  'js/xsolve.js': ['gauss1'],
  'js/spline.js': ['trisol', 'splindOffset', 'splind', 'seval', 'deval', 'sinvrt'],
};
const arrays = {};
for (const name of ('COM1 COM2 C1SAV C2SAV VS1F VS2F VSREZ VSR VSM VSX NBL IBLTE ITRAN TFORCE XSSITR WGAP ACRIT XSTRIP IV_TO_IS IV_TO_IBL X Y S XP YP W1 W2 W3 W4 XOCTR YOCTR TINDEX VMF _MRCHDU_VZTMP _U1M _U2M _D1M _D2M _ULE1M _ULE2M _UTE1M _UTE2M _TRDIF_BLREZ _TRDIF_BLM _TRDIF_BLR _TRDIF_BLX _TRDIF_BTREZ _TRDIF_BTM _TRDIF_BTR _TRDIF_BTX').split(' ')) arrays[name] = 1;
for (const name of ('VS1 VS2 VZ ISYS IPAN VTI XSSI UEDG UINV UINV_A THET DSTR CTAU MASS TAU DIS CTQ DELT TSTR DSSR USLP MDEF PDEF KDEF DIJ HTARG _USAV _MRCHDU_VTMP _TRDIF_BL1 _TRDIF_BL2 _TRDIF_BT1 _TRDIF_BT2').split(' ')) arrays[name] = 2;
for (const name of ['VA', 'VB', 'VDEL', 'VM']) arrays[name] = 3;
const sources = ['_setblSources', '_velocitySources'];
const callbacks = ['blprv', 'blkin', 'blvar', 'blmid', 'trchek', 'tesys', 'blsys', 'hkin', 'syncComToVars', 'syncVarsToCom', 'ensureViscousArrays'];
const booleanFields = ['LALFA', 'LBLINI', 'SIMI', 'TRAN', 'TURB', 'WAKE', 'TRFORC', 'TRFREE'];
const scalarFields = new Set(), recordFields = new Set(), functions = new Map(), files = {}, hashes = {};
function walk(n, fn) {
  if (!n || typeof n !== 'object') return;
  if (n.type) fn(n);
  for (const v of Object.values(n)) if (Array.isArray(v)) v.forEach(x => walk(x, fn)); else if (v && typeof v === 'object') walk(v, fn);
}
for (const [file, names] of Object.entries(selected)) {
  const source = files[file] = readFileSync(new URL(file, root), 'utf8');
  hashes[file] = createHash('sha256').update(source).digest('hex');
  const ast = parse(source, {ecmaVersion: 'latest', sourceType: 'module', locations: true});
  for (const name of names) {
    const node = ast.body.find(n => n.type === 'FunctionDeclaration' && n.id.name === name);
    assert.ok(node, `Missing ${name}`);
    functions.set(name, {node, file});
    walk(node, n => {
      if (n.type === 'MemberExpression' && !n.computed && n.object.type === 'Identifier' && n.object.name === 'ctx') {
        const key = n.property.name;
        if (!(key in arrays) && !sources.includes(key) && !callbacks.includes(key)) scalarFields.add(key);
      }
      if (n.type === 'ObjectExpression') for (const p of n.properties) recordFields.add(p.key.name);
    });
  }
}
scalarFields.add('VACCEL');
for (const name of ['VAF','VBF','VDELF']) arrays[name] = 1;
const scalars = [...scalarFields].sort();
const arrayNames = Object.keys(arrays).sort();
const scalarIndex = new Map(scalars.map((s, i) => [s, i]));
const arrayIndex = new Map(arrayNames.map((s, i) => [s, i]));
const outputRecords = new Set(['axset', 'dampl', 'dampl2', 'hkin', 'dil', 'dilw', 'hsl', 'cfl', 'dit', 'hst', 'cft', 'hct', 'mrcl', 'resetThwaites']);
const outputNumbers = new Set(['dslim', 'seval', 'deval', 'sinvrt']);
const signatures = {
  gauss1: {n:'double', z:'V2', r:'V1'},
  assembleMassRows: Object.fromEntries(['d1M','u1M','d2M','u2M','ule1M','ule2M'].map(n=>[n,'V1'])),
  trisol: Object.fromEntries(['a','b','c','d'].map(n=>[n,'V1'])),
  splindOffset: {x:'V1',xs:'V1',s:'V1'}, splind: {x:'V1',xs:'V1',s:'V1'},
  seval: {x:'V1',xs:'V1',s:'V1'}, deval: {x:'V1',xs:'V1',s:'V1'}, sinvrt: {x:'V1',xs:'V1',s:'V1'},
};
const logs = [];
const integerLoops = [];
const integerBindings = new Set();
let active, env;
function error(n, message) { throw new Error(`${active.file}:${n.loc?.start.line}: ${message} (${n.type})`); }
function field(n) { return n?.type === 'MemberExpression' && !n.computed && n.object.name === 'ctx' ? n.property.name : null; }
function ret(name) { return outputRecords.has(name) ? 'Record' : outputNumbers.has(name) ? 'double' : name === 'velocitySources' ? 'Sources' : 'void'; }
function callName(n) { return n.callee.type === 'Identifier' ? n.callee.name : field(n.callee); }
function type(n) {
  if (!n) return null;
  const f = field(n);
  if (f) return arrays[f] ? `V${arrays[f]}` : sources.includes(f) ? 'Sources' : callbacks.includes(f) ? `fn:${f}` : 'double';
  switch(n.type) {
    case 'Identifier': return env.get(n.name) ?? (functions.has(n.name) ? `fn:${n.name}` : null);
    case 'Literal': return typeof n.value === 'string' ? 'string' : n.value === null ? null : 'double';
    case 'BinaryExpression': case 'UnaryExpression': case 'UpdateExpression': return 'double';
    case 'LogicalExpression': return type(n.left) ?? type(n.right);
    case 'ConditionalExpression': return type(n.consequent) ?? type(n.alternate);
    case 'AssignmentExpression': return type(n.left) ?? type(n.right);
    case 'NewExpression': return n.callee.name === 'Float64Array' || n.callee.name === 'Int32Array' ? 'V1' : null;
    case 'ArrayExpression': return 'V2';
    case 'ObjectExpression': return n.properties.some(p=>p.key.name==='panel') ? 'Sources' : 'Record';
    case 'MemberExpression': {
      const t=type(n.object);
      if(n.computed) return t?.startsWith('V') ? (Number(t.slice(1))===1?'double':`V${Number(t.slice(1))-1}`) : null;
      if(n.property.name==='length') return 'double';
      if(t==='Sources') return n.property.name==='count'?'double':'V1';
      return t==='Record'?'double':null;
    }
    case 'CallExpression': {
      if(n.callee.type==='MemberExpression' && n.callee.property.name==='subarray') return 'V1';
      if(n.callee.type==='MemberExpression' && n.callee.object.name==='Math') return 'double';
      if(n.callee.type==='MemberExpression' && n.callee.object.name==='Array' && n.callee.property.name==='from') return 'V2';
      const name=callName(n);
      if(name==='create2d'||name==='createMatrix1Flat') return 'V2';
      if(functions.has(name)) return ret(name);
      if(env.get(name)?.startsWith('fn:')) return ret(env.get(name).slice(3));
      return null;
    }
    default:return null;
  }
}
function infer(node,name) {
  env=new Map([['ctx','Context'],['ctxIn','Context'],['defaultBlprv','fn:blprv'],['defaultBlkin','fn:blkin']]);
  for(const p of node.params) {
    const key=p.type==='AssignmentPattern'?p.left.name:p.name;
    env.set(key, signatures[name]?.[key] ?? (key==='ctx'||key==='ctxIn'?'Context':'double'));
  }
  const nodes=[];walk(node,n=>nodes.push(n));
  for(let pass=0;pass<12;pass++) for(const n of nodes) {
    if(n.type==='VariableDeclarator' && n.id.type==='Identifier') {
      const t=type(n.init); if(t && t!=='void') env.set(n.id.name,t);
    }
    if(n.type==='VariableDeclarator' && n.id.type==='ObjectPattern') for(const p of n.id.properties)env.set(p.value.name,'double');
    if(n.type==='AssignmentExpression') {
      if(n.left.type==='Identifier' && type(n.right)) env.set(n.left.name,type(n.right));
      if(n.left.type==='MemberExpression' && n.left.computed && n.left.object.type==='Identifier' && type(n.right)?.startsWith('V'))env.set(n.left.object.name,`V${Number(type(n.right).slice(1))+1}`);
    }
    if(n.type==='FunctionDeclaration' && n!==node) env.set(n.id.name,`fn:${n.id.name}`);
  }
}
const cppType=t=>t?.startsWith('fn:')?'auto':t??'double';
function truth(n) {
  if (n.type === 'LogicalExpression' && ['||', '&&'].includes(n.operator)) {
    return `(${truth(n.left)} ${n.operator} ${truth(n.right)})`;
  }
  const t=type(n);
  return t?.startsWith('V')||t==='Sources'?`bool(${expr(n)})`:`truth(${expr(n)})`;
}
function expr(n,expected=null) {
  if(!n)return '';
  const f=field(n);
  if(f) {
    if(arrays[f])return `ctx.arr<${arrays[f]}>(${arrayIndex.get(f)})`;
    if(sources.includes(f))return `ctx.sources(${sources.indexOf(f)})`;
    if(callbacks.includes(f))return f==='ensureViscousArrays'?'unsupported':`fn_${f}`;
    return `ctx.get(${scalarIndex.get(f)})`;
  }
  const e=x=>expr(x);
  switch(n.type) {
    case 'Identifier':
      // JS arithmetic remains double precision even when a loop's storage
      // has a proven integer representation. Only direct array subscripts
      // below use that representation without a double conversion.
      if(integerBindings.has(n.name))return `double(${n.name})`;
      if(n.name==='undefined'||n.name==='NaN')return 'nan_value()';
      if(n.name==='defaultBlprv')return 'fn_blprv';
      if(n.name==='defaultBlkin')return 'fn_blkin';
      if(n.name==='NCOM')return '73.0';
      if(functions.has(n.name) && !env.has(n.name))return `fn_${n.name}`;
      return n.name;
    case 'Literal':
      if(n.value===null)return expected?.startsWith('V')?`${expected}{}`:'nan_value()';
      if(typeof n.value==='string')return JSON.stringify(n.value);
      if(typeof n.value==='boolean')return n.value?'1.0':'0.0';
      return Number.isInteger(n.value)?`${n.value}.0`:n.raw;
    case 'ThisExpression':return 'ctxIn';
    case 'MemberExpression':
      if(n.computed && n.property.type==='Identifier' && integerBindings.has(n.property.name))return `${e(n.object)}.at(${n.property.name})`;
      if(n.computed)return `${e(n.object)}[${e(n.property)}]`;
      if(type(n.object)==='Sources')return n.property.name==='count'?`${e(n.object)}.count()`:`${e(n.object)}.arr(${['panel','system','sign','mass','sensitivity'].indexOf(n.property.name)})`;
      return `${e(n.object)}.${n.property.name}`;
    case 'BinaryExpression': {
      if(n.left.type==='UnaryExpression'&&n.left.operator==='typeof') {
        if(!type(n.left.argument)?.startsWith('fn:') || n.right.value!=='function'
          || !['===','!=='].includes(n.operator))error(n,'Unsupported typeof comparison');
        return n.operator==='==='?'true':'false';
      }
      if(n.operator==='**')return `host_pow(${e(n.left)}, ${e(n.right)})`;
      const op={'===':'==','!==':'!='}[n.operator]??n.operator;
      if(!['+','-','*','/','<','>','<=','>=','==','!='].includes(op))error(n,'Unsupported binary operator');
      return `(${e(n.left)} ${op} ${e(n.right)})`;
    }
    case 'LogicalExpression':
      if(n.operator==='??') {
        if(n.left.type==='LogicalExpression'&&n.left.operator==='??')return expr({...n,left:n.left.left,right:{...n,left:n.left.right,right:n.right}});
        const k=field(n.left);
        if(k)return `(ctx.has(${scalarIndex.get(k)}) ? ${e(n.left)} : ${e(n.right)})`;
        if(n.left.type==='MemberExpression'&&n.left.computed)return `${e(n.left.object)}.get_or(${e(n.left.property)}, ${e(n.right)})`;
        error(n,'Unsupported nullish expression');
      }
      if(n.left.name==='ctxIn')return 'ctxIn';
      if(type(n.left)!=='double'||type(n.right)!=='double')error(n,'Unsupported logical value types');
      // JS short circuit operators return an operand, not a boolean. Evaluate
      // the left operand once and evaluate the right only when required.
      return `([&]() -> double { double value = ${e(n.left)}; return truth(value) ? ${n.operator==='||'?'value':e(n.right)} : ${n.operator==='||'?e(n.right):'value'}; }())`;
    case 'UnaryExpression':return n.operator==='!'?`(!${truth(n.argument)})`:`(${n.operator}${e(n.argument)})`;
    case 'ConditionalExpression':return `(${truth(n.test)} ? ${e(n.consequent)} : ${e(n.alternate)})`;
    case 'SequenceExpression':return `(${n.expressions.map(e).join(', ')})`;
    case 'AssignmentExpression': {
      if(n.left.type==='ObjectPattern')return `([&]() { auto value = ${e(n.right)}; ${n.left.properties.map(p=>`${e(p.value)} = value.${p.key.name};`).join(' ')} return value; }())`;
      const k=field(n.left);const rhs=expr(n.right,type(n.left));
      if(k&&scalarIndex.has(k))return `ctx.put(${scalarIndex.get(k)}, ${n.operator==='='?rhs:`(${e(n.left)} ${n.operator.slice(0,-1)} ${rhs})`})`;
      return `(${e(n.left)} ${n.operator} ${rhs})`;
    }
    case 'UpdateExpression':return n.prefix?`${n.operator}${e(n.argument)}`:`${e(n.argument)}${n.operator}`;
    case 'NewExpression': {
      if(n.callee.type!=='Identifier'||!['Float64Array','Int32Array','Array'].includes(n.callee.name)
        || n.arguments.length!==1)error(n,'Unsupported allocation');
      const rank=expected?.startsWith('V')?expected.slice(1):'1';
      return `make_vec<${rank}>(${e(n.arguments[0])})`;
    }
    case 'ArrayExpression':return `make_rows(${n.elements.map(x=>x.value===null?'V1{}':expr(x,'V1')).join(', ')})`;
    case 'ObjectExpression': {
      const t=type(n);
      if(t==='Sources')return `([&]() { auto value = make_sources(); ${n.properties.map(p=>`value.${p.key.name==='count'?'count()':`arr(${['panel','system','sign','mass','sensitivity'].indexOf(p.key.name)})`} = ${e(p.value)};`).join(' ')} return value; }())`;
      return `([&]() { ${t} value{}; ${n.properties.map(p=>`value.${p.key.name} = ${e(p.value)};`).join(' ')} return value; }())`;
    }
    case 'CallExpression': {
      const name=callName(n);
      if(n.callee.type==='MemberExpression'&&n.callee.object.name==='Array'&&n.callee.property.name==='from') {
        const size=n.arguments[0].properties.find(p=>p.key.name==='length').value;
        const body=n.arguments[1].body;
        assert.equal(body.type,'NewExpression');
        return `make_matrix(${e(size)}, ${e(body.arguments[0])})`;
      }
      const args=n.arguments.map(e);
      if(name==='ensureCtx'||name==='ensureDefectArrays')return '((void)0)';
      if(name==='checkN')return `check_n(${args[0]})`;
      if(name==='create2d'||name==='createMatrix1Flat')return `make_matrix(${args[0]} + 1, ${args[1]} + 1)`;
      if(n.callee.type==='MemberExpression'&&n.callee.object.name==='Number') {
        assert.equal(n.callee.property.name,'isFinite');return `finite(${args[0]})`;
      }
      if(n.callee.type==='MemberExpression'&&n.callee.object.name==='Math') {
        const m=n.callee.property.name;
        if(m==='min'||m==='max')return args.slice(1).reduce((a,b)=>`js_${m}(${a}, ${b})`,args[0]);
        if(['abs','sqrt','floor'].includes(m))return `__builtin_${m==='abs'?'fabs':m}(${args.join(', ')})`;
        if(['exp','log','log10','sin','cos','pow','atan','tanh'].includes(m))return `host_${m}(${args.join(', ')})`;
        error(n,'Unsupported Math function');
      }
      if(n.callee.type==='MemberExpression'&&['fill','subarray'].includes(n.callee.property.name))return `${e(n.callee.object)}.${n.callee.property.name}(${args.join(', ')})`;
      if(name==='ensureViscousArrays')return 'unsupported()';
      if(functions.has(name))return `fn_${name}(${args.join(', ')})`;
      if(n.callee.type==='Identifier'&&env.get(name)?.startsWith('fn:'))return `${name}(${args.join(', ')})`;
      error(n,`Unsupported call ${name??n.callee.type}`);
      break;
    }
    default:error(n,'Unsupported expression');
  }
}
function logging(n) {
  assert.equal(n.arguments.length,1);
  const template=n.arguments[0];assert.equal(template.type,'TemplateLiteral');
  let js=files[active.file].slice(template.start,template.end); const values=[]; const replacements=[];
  for(const expression of template.expressions) {
    let base=expression;
    while(base.type==='CallExpression'&&base.callee.type==='MemberExpression')base=base.callee.object;
    if(base.type==='CallExpression'&&base.callee.name==='String')base=base.arguments[0];
    replacements.push({start:base.start-template.start,end:base.end-template.start,text:`values[${values.length}]`});
    values.push(expr(base));
  }
  for(const r of replacements.reverse())js=js.slice(0,r.start)+r.text+js.slice(r.end);
  const id=logs.length;logs.push(`(values) => ${js}`);
  assert.ok(values.length<=3);while(values.length<3)values.push('0.0');
  return `record_log(${id}, ${values.join(', ')});`;
}
function declaration(n,forLoop=false) {
  return n.declarations.map(d=>{
    if(d.id.type==='ObjectPattern') {
      const tmp=`record_${d.start}`;
      return `auto ${tmp} = ${expr(d.init)};\n`+d.id.properties.map(p=>`double ${p.value.name} = ${tmp}.${p.key.name};`).join('\n');
    }
    const t=env.get(d.id.name)??'double';
    const value=d.init?expr(d.init,t):t==='double'?'nan_value()':`${t}{}`;
    return `${cppType(t)} ${d.id.name} = ${value}${forLoop?'':';'}`;
  }).join(forLoop?', ':'\n');
}
function integerLoop(n) {
  // A nonnegative integer literal advanced by exactly one stays integral.
  // Reject any body that assigns/shadows the binding, or contains a nested
  // function whose captures would need separate lexical analysis. Overflow
  // traps transactionally, so it can never wrap to a different JS index.
  if(n.init?.type!=='VariableDeclaration' || n.init.kind!=='let' || n.init.declarations.length!==1)return null;
  const d=n.init.declarations[0];
  if(d.id.type!=='Identifier' || d.init?.type!=='Literal' || !Number.isInteger(d.init.value)
    || d.init.value<0 || d.init.value>0xffffffff)return null;
  const name=d.id.name,u=n.update;
  const increment=u?.type==='UpdateExpression' && u.operator==='++' && u.argument.type==='Identifier' && u.argument.name===name
    || u?.type==='AssignmentExpression' && u.operator==='+=' && u.left.type==='Identifier' && u.left.name===name
      && u.right.type==='Literal' && u.right.value===1;
  if(!increment)return null;
  let safe=true;
  walk(n.body,child=>{
    if(child.type==='AssignmentExpression' && child.left.type==='Identifier' && child.left.name===name)safe=false;
    if(child.type==='UpdateExpression' && child.argument.type==='Identifier' && child.argument.name===name)safe=false;
    if(child.type==='VariableDeclarator')walk(child.id,id=>{if(id.type==='Identifier' && id.name===name)safe=false;});
    if(['FunctionDeclaration','FunctionExpression','ArrowFunctionExpression'].includes(child.type))safe=false;
  });
  // The condition is also allowed to be an expression in the JS subset;
  // ensure it cannot modify or shadow the integer binding either.
  walk(n.test,child=>{
    if(child.type==='AssignmentExpression' || child.type==='UpdateExpression'
      || ['FunctionExpression','ArrowFunctionExpression'].includes(child.type))safe=false;
  });
  return safe?d:null;
}
function statement(n) {
  switch(n.type) {
    case 'BlockStatement':return `{\n${n.body.map(statement).join('\n')}\n}`;
    case 'VariableDeclaration':return declaration(n);
    case 'ExpressionStatement':
      if(n.expression.type==='CallExpression'&&n.expression.callee.object?.name==='console')return logging(n.expression);
      return expr(n.expression)+';';
    case 'IfStatement':return `if (${truth(n.test)}) ${statement(n.consequent)}${n.alternate?' else '+statement(n.alternate):''}`;
    case 'ForStatement': {
      const d=integerLoop(n);
      if(!d)return `for (${n.init?.type==='VariableDeclaration'?declaration(n.init,true):expr(n.init)}; ${n.test?truth(n.test):''}; ${expr(n.update)}) ${statement(n.body)}`;
      const name=d.id.name,wasInteger=integerBindings.has(name);
      integerBindings.add(name);
      const result=`for (u32 ${name} = ${d.init.value}u; ${n.test?truth(n.test):''}; advance_index(${name})) ${statement(n.body)}`;
      if(!wasInteger)integerBindings.delete(name);
      integerLoops.push({file:active.file,line:n.loc.start.line,name});
      return result;
    }
    case 'WhileStatement':return `while (${truth(n.test)}) ${statement(n.body)}`;
    case 'ReturnStatement':return `return${n.argument?' '+expr(n.argument):''};`;
    case 'BreakStatement':return 'break;';
    case 'ContinueStatement':return 'continue;';
    case 'ThrowStatement':return 'unsupported();';
    case 'EmptyStatement':return ';';
    case 'FunctionDeclaration':return `auto ${n.id.name} = [&](${n.params.map(p=>`double ${p.name}`).join(', ')}) -> ${ret(n.id.name)} ${statement(n.body)};`;
    default:error(n,'Unsupported statement');
  }
}
function signature(name,node,defaults=false) {
  return `${ret(name)} fn_${name}(${node.params.map(p=>{
    const key=p.type==='AssignmentPattern'?p.left.name:p.name;
    return `${signatures[name]?.[key]??(key==='ctx'||key==='ctxIn'?'Context':'double')} ${key}${defaults&&p.type==='AssignmentPattern'?' = '+expr(p.right):''}`;
  }).join(', ')})`;
}
const prototypes=[],bodies=[];
for(const [name,entry]of functions) {
  active=entry;infer(entry.node,name);
  prototypes.push(signature(name,entry.node,true)+';');
  bodies.push(`// ${entry.file}:${entry.node.loc.start.line} (${name})\n${signature(name,entry.node)} ${statement(entry.node.body)}`);
}
const prefix=readFileSync(new URL('./boundary_layer_runtime.h',import.meta.url),'utf8');
const records=`struct Record {\n${[...recordFields].filter(n=>!['panel','system','sign','mass','sensitivity','count'].includes(n)).sort().map(n=>`double ${n} = nan_value();`).join('\n')}\n};`;
const accessors = 'enum : u32 {' + scalars.map((name,i)=>`scalar_${name}=${i}`).join(',') + '};\n'
  + 'enum : u32 {' + arrayNames.map((name,i)=>`array_${name}=${i}`).join(',') + '};\n';
const step = readFileSync(new URL('./boundary_layer_step.h',import.meta.url),'utf8');
const code=`// Generated from current JavaScript by generate_boundary_layer.mjs.\n// Derived from XFOIL by Mark Drela. GNU General Public License v2 or later.\n// https://web.mit.edu/drela/Public/web/xfoil/\n${prefix}\n${records}\n${prototypes.join('\n')}\n${bodies.join('\n\n')}\n${accessors}\n${step}\nextern "C" void boundary_layer(unsigned scalar, unsigned array, unsigned source, unsigned begin, unsigned end) {\n arena_begin(begin,end);\n fn_setbl(Context{(Slot*)scalar, (Descriptor*)array, (Sources*)source});\n}\nextern "C" void boundary_layer_step(unsigned scalar, unsigned array, unsigned source, unsigned begin, unsigned end) {\n arena_begin(begin,end);\n Context ctx{(Slot*)scalar, (Descriptor*)array, (Sources*)source};\n fn_setbl(ctx);\n solve_boundary_layer_matrix(ctx);\n}\n`;
const metadata={hashes,scalars,arrays:Object.fromEntries(arrayNames.map(n=>[n,arrays[n]])),sources,booleanFields,callbacks,functions:[...functions.keys()],integerLoops};
writeFileSync(new URL('./boundary_layer_generated.cpp',import.meta.url),code);
writeFileSync(new URL('../js/boundary_layer_layout.js',import.meta.url),`// Generated ABI and diagnostic formatting. GPL-2.0-or-later.\nexport const layout = ${JSON.stringify(metadata,null,2)};\nexport const logFormats = [\n${logs.join(',\n')}\n];\n`);
process.stdout.write(`Generated ${functions.size} functions, ${scalars.length} scalar fields, ${arrayNames.length} array fields; ${integerLoops.length} proven integer loops.\n`);
