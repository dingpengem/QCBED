"""Read-only fixed-form source inventory; no compiler or binary evidence.

Handles tab-form, continuation records and long lines in this source tree.
Expression calls are recognized from known functions and EXTERNAL declarations;
this is not a full Fortran compiler or an ABI-equivalence checker.
"""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
UNIT = re.compile(r'^(?:(REAL(?:\*\d+)?|DOUBLE PRECISION|INTEGER(?:\*\d+)?|COMPLEX(?:\*\d+)?|LOGICAL)\s+)?(PROGRAM|SUBROUTINE|FUNCTION|BLOCK DATA)\s*(\w+)?\s*(?:\((.*)\))?$', re.I)
DECL = re.compile(r'^(DOUBLE PRECISION|REAL(?:\*\d+)?|INTEGER(?:\*\d+)?|COMPLEX(?:\*\d+)?|LOGICAL|CHARACTER(?:\*\d+)?)\s+(.+)$', re.I)
INTRINSICS = set('ABS DABS DSQRT SQRT DBLE REAL INT IDINT MOD DMOD EXP DEXP LOG DLOG SIN COS DSIN DCOS ASIN DASIN ATAN ATAN2 ACOS DACOS DCMPLX CMPLX DIMAG AIMAG CONJG MAX MIN MAX0 MIN0 SIGN ISHFT IAND IOR TIME ETIME DTIME'.split())


def split_args(s):
    parts, start, depth, quote = [], 0, 0, None
    for i, c in enumerate(s):
        if quote:
            if c == quote: quote = None
        elif c in "'\"": quote = c
        elif c == '(': depth += 1
        elif c == ')': depth -= 1
        elif c == ',' and depth == 0:
            parts.append(s[start:i].strip()); start = i+1
    if s[start:].strip(): parts.append(s[start:].strip())
    return parts


def statements(text):
    current = None
    for n, raw in enumerate(text.splitlines(), 1):
        if not raw.strip() or raw[0] in 'cC*!': continue
        line = raw.expandtabs(8)
        continuation = len(line) > 5 and line[5] not in ' 0'
        body = line[6:] if len(line) > 6 else ''
        # Strip inline comments, preserving quoted literals.
        body = re.split(r"!(?=(?:[^']*'[^']*')*[^']*$)", body, maxsplit=1)[0].strip()
        if not body: continue
        if continuation and current:
            current['text'] += ' '+body; current['end_line'] = n
        else:
            if current: yield current
            current = {'line': n, 'end_line': n, 'text': body}
    if current: yield current


def declarations(stmts):
    result = {}
    for s in stmts:
        m = DECL.match(s['text'])
        if m and not re.search(r'\bFUNCTION\b', m[2], re.I):
            for item in split_args(m[2]):
                name = re.match(r'\w+', item)
                if name: result[name[0].upper()] = {'type': m[1].upper(), 'declaration': item, 'line': s['line']}
    return result


def expression_type(expr, unit):
    """Conservative visible scalar/element type; no compiler flag assumptions."""
    e = expr.strip().upper()
    if re.fullmatch(r'[+-]?\d+', e): return 'INTEGER'
    if re.fullmatch(r'[+-]?(?:\d+\.\d*|\.\d+)(?:[ED][+-]?\d+)?', e):
        return 'REAL*8' if 'D' in e else 'REAL'
    m = re.match(r'(DBLE|REAL|INT|IDINT)\s*\(', e)
    if m: return {'DBLE':'REAL*8','REAL':'REAL','INT':'INTEGER','IDINT':'INTEGER'}[m[1]]
    m = re.fullmatch(r'(\w+)(?:\(.*\))?', e)
    if m:
        return unit['declarations'].get(m[1],{}).get('type', 'implicit INTEGER' if m[1][0] in 'IJKLMN' else 'implicit REAL')
    return 'expression: unresolved'


def build_manifest(root=ROOT):
    files, units = [], []
    source_root = root/'legacy_fortran'
    for path in sorted(source_root.rglob('*'), key=lambda p: str(p).lower()):
        if path.suffix.lower() not in ('.f', '.for', '.ftn'): continue
        data = path.read_bytes()
        info = {'file': path.relative_to(root).as_posix(), 'sha256': hashlib.sha256(data).hexdigest(), 'nul': b'\0' in data}
        files.append(info)
        if info['nul'] or data.startswith(b'\x7fELF'):
            info.update(readable=False, encoding=None, lines=None, units=[], truncation_flags=[], note='binary content; not parsed as source')
            continue
        try: text = data.decode('utf-8-sig'); info['encoding'] = 'utf-8'
        except UnicodeDecodeError:
            text = data.decode('cp1252'); info['encoding'] = 'cp1252'
        info.update(readable=not info['nul'], lines=len(text.splitlines()), units=[], truncation_flags=[])
        if info['nul']: continue
        active = None
        for stmt in statements(text):
            match = UNIT.match(stmt['text'])
            if match:
                if active: info['truncation_flags'].append('new unit before END')
                kind, name = match[2].upper(), (match[3] or '<unnamed>').upper()
                active = {'name': name, 'type': ((match[1]+' ') if match[1] else '')+kind,
                          'file': info['file'], 'start_line': stmt['line'], 'end_line': None,
                          'arguments': split_args((match[4] or '').upper()), 'signature': stmt['text'],
                          'filename_matches': path.stem.upper() == name, 'statements': []}
                units.append(active); info['units'].append(name)
            if active:
                active['statements'].append(stmt)
                if re.fullmatch(r'END(?:\s+(?:SUBROUTINE|FUNCTION|PROGRAM)(?:\s+\w+)?)?', stmt['text'], re.I):
                    active['end_line'] = stmt['end_line']; active = None
        if active: info['truncation_flags'].append('unit missing terminal END')
        if not info['units']:
            info['note'] = 'declaration include, not a program unit' if path.suffix.lower()=='.ftn' else 'no recognizable unit'
    include_stmts = list(statements((source_root/'MULTIS.FTN').read_text()))
    functions = {u['name'] for u in units if 'FUNCTION' in u['type']}
    by_name = {}
    calls = []
    for u in units:
        by_name.setdefault(u['name'], []).append(u)
        own = declarations(u['statements'])
        u['includes'] = [s for s in u['statements'] if re.match(r'INCLUDE\b',s['text'],re.I)]
        u['declarations'] = (declarations(include_stmts) if u['includes'] else {}) | own
        u['external'] = [v.strip().upper() for s in u['statements'] if re.match(r'EXTERNAL\s',s['text'],re.I) for v in s['text'].split(None,1)[1].split(',')]
        u['argument_types'] = {a:u['declarations'].get(a, {'type':'implicit/unresolved'}) for a in u['arguments']}
        if 'FUNCTION' in u['type']:
            u['return_type'] = u['declarations'].get(u['name'],{}).get('type', u['type'].replace('FUNCTION','').strip() or 'implicit/unresolved')
        for s in u['statements'][1:]:
            t = s['text']
            if DECL.match(t) or re.match(r'(EXTERNAL|DIMENSION|COMMON|PARAMETER|FORMAT|INCLUDE)\b',t,re.I): continue
            # Quoted output text is never a call expression.
            t = re.sub(r"'[^']*'|\"[^\"]*\"", "''", t)
            patterns = [(m, 'CALL') for m in re.finditer(r'\bCALL\s+(\w+)\s*',t,re.I)]
            for m in re.finditer(r'\b(\w+)\s*(?=\()',t):
                name = m[1].upper()
                if name in functions | set(u['external']) | INTRINSICS and not re.search(r'CALL\s+$',t[:m.start()],re.I):
                    patterns.append((m, 'function'))
            for m, syntax in patterns:
                pos=m.end(); args=[]
                if pos < len(t) and t[pos]=='(':
                    depth=1; end=pos+1
                    while end<len(t) and depth:
                        depth += (t[end]=='(')-(t[end]==')'); end+=1
                    args=split_args(t[pos+1:end-1])
                name=m[1].upper()
                calls.append({'caller':u['name'],'caller_file':u['file'],'line':s['line'], 'end_line':s['end_line'],
                              'callee':name,'syntax':syntax,'arguments':args,'argument_count':len(args),
                              'callback':name in u['arguments'] and name in u['external']})
        normalized='\n'.join(re.sub(r'\s+','',s['text']).upper() for s in u['statements'])
        u['statement_sha256']=hashlib.sha256(normalized.encode()).hexdigest()
    for c in calls:
        defs=by_name.get(c['callee'],[])
        caller = next(u for u in by_name[c['caller']] if u['file']==c['caller_file'])
        c['actual_types'] = [expression_type(a,caller) for a in c['arguments']]
        c['formal_types'] = [[u['argument_types'][a]['type'] for a in u['arguments']] for u in defs]
        c['visible_kind_disagreements'] = []
        for u in defs:
            for i,(actual,formal) in enumerate(zip(c['actual_types'], c['formal_types'][0])):
                if actual in ('INTEGER','INTEGER*8','REAL','REAL*8','COMPLEX','COMPLEX*16') and formal in ('INTEGER','INTEGER*8','REAL','REAL*8','COMPLEX','COMPLEX*16') and actual != formal:
                    c['visible_kind_disagreements'].append({'argument':i+1,'actual':actual,'formal':formal})
        c['caller_return_declaration'] = caller['declarations'].get(c['callee']) if c['syntax']=='function' else None
        c['resolution']= ('uncertain: callback formal' if c['callback'] else
                         'definition found exactly once' if len(defs)==1 else
                         'definition found multiple times' if defs else
                         'intrinsic/runtime/library routine' if c['callee'] in INTRINSICS else 'definition missing')
        c['signature_status']= ('callback interface: inspect binding' if c['callback'] else
                               'argument-count mismatch' if defs and any(len(u['arguments'])!=c['argument_count'] for u in defs) else
                               'potentially kind-mismatched' if c['visible_kind_disagreements'] else
                               'count compatible; kinds require review' if defs else c['resolution'])
    for u in units:
        u['callers']=[c for c in calls if c['callee']==u['name']]
        u['callback_bindings']=[c for c in calls if u['name'] in [a.upper() for a in c['arguments']]]
    return {'scope':'current readable source only; extended fixed form; not a compiler', 'files':files,
            'routines':{k:v for k,v in sorted(by_name.items())},'calls':calls}


if __name__=='__main__':
    manifest=build_manifest()
    (ROOT/'docs/legacy_source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f"{len(manifest['files'])} files; {sum(map(len,manifest['routines'].values()))} units; {len(manifest['calls'])} calls")
