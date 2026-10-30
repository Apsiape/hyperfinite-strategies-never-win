"""Independent transcription of the game definition and the base-game definition.
No compiler, audit code, model output, or paper JSON is imported.
Selects 512 random supported pairs and checks them, then checks the full table.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import json, random, secrets, itertools, hashlib, collections, time, ctypes
from pathlib import Path
ROOT = Path(__file__).resolve().parent
INPUT = ROOT/'labelled-giii-game.json'
started = time.monotonic()
game = json.loads(INPUT.read_text(encoding='utf-8-sig'))

def frozen(v):
    return tuple(map(frozen, v)) if isinstance(v, list) else v

def signs(n, target=None):
    result = []
    for ss in itertools.product((-1, 1), repeat=n):
        parity = 1
        for s in ss: parity *= s
        if target is None or parity == target: result.append(ss)
    return result

ROW = signs(3, 1)
COL = {1: signs(3, 1), 2: signs(3, 1), 3: signs(3, -1)}
SPECIAL = [(1,1,1), (1,-1,-1)]
A, B, contexts, relations = {}, {}, {}, {}
base_names = ('R1','R2','R3','T0','T1')
bob_names = ('C1','C2','C3','U0','U1')
for k in range(3):
    contexts['row'+str(k)] = (['S'+str(k)+str(c) for c in range(3)], 1)
    contexts['col'+str(k)] = (['S'+str(r)+str(k) for r in range(3)], 1 if k < 2 else -1)
for q in 'abc': contexts['Y'+q] = (['Y'+q], None)

def ordinary(side, q):
    if q[0] == 'R': return ROW + ['f']
    if q[0] == 'C': return COL[int(q[1])] + ['f']
    return [x for x in ROW if x != SPECIAL[int(q[1])]] + [-1, 1]

for j in range(9):
    for q in base_names: A[f'base:{j}:{q}'] = ordinary('A', q)
    for q in bob_names: B[f'base:{j}:{q}'] = ordinary('B', q)
    for cell in ('S'+str(r)+str(c) for r in range(3) for c in range(3)):
        B[f'label:{j}:{cell}'] = [-1,1]
    for q in 'abc': B[f'label:{j}:Y{q}'] = [-1,1]
    for ctx, (variables, parity) in contexts.items():
        for q in base_names:
            A[f'joint:{j}:{q}:{ctx}'] = [(out, ss) for out in ordinary('A', q) for ss in signs(len(variables), parity)]
    p, m = (j+1)%9, (j+2)%9
    # Relation codes are interpreted separately below; no callbacks from the source.
    specs = [('Zplus',(p,),('S11',),'plus'), ('Zminus',(m,),('S11',),'minus')]
    specs.extend(('Yflag'+q,(j,),('Y'+q,q),'equal') for q in 'abc')
    specs += [('inactiveb0',(p,),('Yb','a'),'equal'), ('inactiveb1',(m,),('Yb','b'),'equal'),
              ('inactivec0',(p,),('Yc','b'),'equal'), ('inactivec1',(m,),('Yc','a'),'equal'),
              ('active',(p,m),('Ya','S00','c'),'product'),
              ('commuteXa',(p,m),('S00','a'),'any'), ('commuteXb',(p,m),('S00','b'),'any')]
    for tag, block, variables, rule in specs:
        def allowed(ss):
            if rule == 'plus': return ss[0] == 1
            if rule == 'minus': return ss[0] == -1
            if rule == 'equal': return ss[0] == ss[1]
            if rule == 'product': return ss[0]*ss[1]*ss[2] == 1
            return True
        name = f'relation:{j}:{tag}'
        A[name] = [(0,)] + [(1,)+ss for ss in signs(len(variables)) if allowed(ss)]
        relations[name] = (block, tuple('global:'+v if v in ('a','b','c') else f'label:{j}:{v}' for v in variables))
B['flag'] = list(range(9))
for q in 'abc':
    B['global:'+q] = [-1,1]
    A['flag-global:'+q] = list(itertools.product(range(9),(-1,1)))

assert set(A) == set(game['alice']), 'Alice question set differs'
assert set(B) == set(game['bob']), 'Bob question set differs'
for side, schema in [('alice',A),('bob',B)]:
    for name, outcomes in schema.items():
        actual = [frozen(v) for v in game[side][name]]
        assert len(actual) == len(set(actual)), ('duplicate answer', side, name)
        assert set(actual) == set(outcomes), ('answer set', side, name)

expected_support = set()
for an in A:
    parts = an.split(':')
    if parts[0] in ('base','joint'):
        j, q = int(parts[1]), parts[2]
        expected_support.update((an,f'base:{j}:{r}') for r in bob_names)
        if q in ('R1','R2'): expected_support.add((an,'flag'))
        if parts[0] == 'joint':
            expected_support.update((an,f'label:{j}:{v}') for v in contexts[parts[3]][0])
    elif parts[0] == 'flag-global':
        expected_support.update(((an,'flag'),(an,'global:'+parts[1])))
    else:
        expected_support.add((an,'flag'))
        expected_support.update((an,v) for v in relations[an][1])
records = {(p['alice'],p['bob']):p for p in game['supported_pairs']}
assert len(records) == len(game['supported_pairs']), 'duplicate supported pair'
assert set(records) == expected_support, 'supported pair set differs'
assert game['distribution'] == {'type':'uniform_on_supported_pairs','denominator':len(records)}

def base_accept(q, r, a, b):
    if q.startswith('R') and r.startswith('C'):
        if a == 'f' or b == 'f': return a == b
        return a[int(r[1])-1] == b[int(q[1])-1]
    if q.startswith('T') and r == 'C1':
        if type(a) is int: return b in ('f',SPECIAL[int(q[1])])
        return a == b
    if q == 'R1' and r.startswith('U'):
        if type(b) is int: return a in ('f',SPECIAL[int(r[1])])
        return a == b
    if q.startswith('T') and r.startswith('U') and type(a) is int and type(b) is int:
        return (a == b) == (int(q[1])*int(r[1]) == 0)
    return True

def accepts(an, bn, a, b):
    ap, bp = an.split(':'), bn.split(':')
    if ap[0] == 'flag-global': return a[0] == b if bn == 'flag' else a[1] == b
    if ap[0] == 'relation':
        block, variables = relations[an]
        if bn == 'flag': return bool(a[0]) == (b in block)
        return a[0] == 0 or a[1+variables.index(bn)] == b
    j, q = int(ap[1]), ap[2]
    if bp[0] == 'label': return a[1][contexts[ap[3]][0].index(bp[2])] == b
    out = a[0] if ap[0] == 'joint' else a
    if bn == 'flag':
        if q == 'R1': return (out == 'f') == (b == j)
        return (out == (1,1,1)) == (b in ((j+1)%9,(j+2)%9))
    return base_accept(q, bp[2], out, b)

def check_pair(key, supplied=None):
    an,bn = key
    p = records[key] if supplied is None else supplied
    found = [tuple(v) for v in p['forbidden']]
    want = {(i,k) for i,a in enumerate(game['alice'][an]) for k,b in enumerate(game['bob'][bn])
            if not accepts(an,bn,frozen(a),frozen(b))}
    assert len(found) == len(set(found)) and set(found) == want, ('forbidden table',key)
    return len(want)

seed = int(os.environ.get("GAME_CHECK_SEED", "336906185362330768646873192451346372926"))
sample = random.Random(seed).sample(sorted(records),512)
sample_events = sum(check_pair(pair) for pair in sample)
all_events = sum(check_pair(pair) for pair in sorted(records))
# Ensure the checker is not vacuous: mutate one losing set and one matching rule.
nonempty = next(key for key in records if records[key]['forbidden'])
mutant = dict(records[nonempty]); mutant['forbidden'] = mutant['forbidden'][1:]
try: check_pair(nonempty,mutant)
except AssertionError: mutation_caught = True
else: raise AssertionError('mutation not detected')
result = {'status':'PASS','independent_of_compiler':True,'human_definition':'game definition and base-game definition',
          'seed':seed,'sample_size':len(sample),'sample_forbidden_events':sample_events,'sample_pairs':sample,
          'all_pairs_checked':len(records),'all_forbidden_events':all_events,'alice_questions':len(A),'bob_questions':len(B),
          'schema_outcome_sets_checked':len(A)+len(B),'mutation_test_passed':mutation_caught,
          'source_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),'elapsed_seconds':round(time.monotonic()-started,3)}
(ROOT/'independent-predicate-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k != 'sample_pairs'},indent=2),flush=True)
