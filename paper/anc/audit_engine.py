"""Independent action audit of the symbolic reduction and table structure.
Does not import/execute the compiler or use its act/mono functions for expected values.
"""
import ast, itertools, random, json, hashlib, ctypes
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parent
src=(ROOT/'compile_and_check.py').read_text(encoding='utf-8-sig')
tree=ast.parse(src)
# Load only the engine, stopping before construction of the model.
cut=next(i for i,n in enumerate(tree.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='R' for t in n.targets))
engine=ast.Module(body=tree.body[:cut],type_ignores=[])
ns={'__file__':str(ROOT/'compile_and_check.py')}; exec(compile(engine,'engine-prefix','exec'),ns)
OUT={'a':(1,0),'b':(0,1),'c':(0,1)}; RES={'a':('c','c'),'b':('a','b'),'c':('b','a')}
def automaton(q,w):
    if not w:return ()
    return (OUT[q][w[0]],)+automaton(RES[q][w[0]],w[1:])
def group_action(g,w):
    for q in g[::-1]: w=automaton(q,w)
    return w

def apply(m,w):
    mu,g,nu=m
    return mu+group_action(g,w[len(nu):]) if w[:len(nu)]==nu else None
rng=random.Random(30102026)
def word(n): return tuple(rng.randrange(2) for _ in range(n))
def monomial(): return (word(rng.randrange(4)),''.join(rng.choice('abc') for _ in range(rng.randrange(5))),word(rng.randrange(4)))
trials=0
for _ in range(600):
    a,b=monomial(),monomial(); c=ns['mono'](*a,*b)
    for w in itertools.product((0,1),repeat=7):
        v=apply(b,w); left=None if v is None else apply(a,v)
        right=None if c is None else apply(c,w)
        assert left==right,(a,b,c,w,left,right)
        trials+=1
# The normal form's Cuntz expansion is checked on words long enough for every cylinder.
for _ in range(120):
    m=monomial(); mu,g,nu=m; D=len(nu)+rng.randrange(4)
    for tail in itertools.product((0,1),repeat=D-len(nu)):
        ow,sec=ns['act'](g,tail)
        for w in itertools.product((0,1),repeat=5):
            assert group_action(g,tail+w)==ow+group_action(sec,w)
            trials+=1
# Independent validation of JSON references, outcome indices and uniform distribution.
game=json.loads((ROOT/'labelled-giii-game.json').read_text())
A,B=game['alice'],game['bob']; pairs=game['supported_pairs']
assert len(A)==561 and len(B)==157 and len(pairs)==3696
seen=set(); bad=0
for p in pairs:
    a,b=p['alice'],p['bob']; assert a in A and b in B
    assert (a,b) not in seen;seen.add((a,b))
    assert len({tuple(t) for t in p['forbidden']})==len(p['forbidden'])
    for ia,ib in p['forbidden']: assert 0<=ia<len(A[a]) and 0<=ib<len(B[b])
    bad+=len(p['forbidden'])
assert bad==90630 and game['distribution']['denominator']==3696
# Evidence that even a nontrivial group action is not erased by the zero tester.
assert ns['normal'](ns['sub'](ns['UG']['a'],ns['ID']))
# Portable peak-memory observation.
from resource_usage import peaks
from types import SimpleNamespace
_peak,_commit=peaks()
assert _peak<1_000_000_000 and _commit<1_000_000_000
pm=SimpleNamespace(PeakWorkingSetSize=_peak,PeakPagefileUsage=_commit)
r={'status':'PASS','independent_action_tests':trials,'table_pairs':len(pairs),'forbidden_events':bad,'peak_working_set_bytes':pm.PeakWorkingSetSize,'game_sha256':hashlib.sha256((ROOT/'labelled-giii-game.json').read_bytes()).hexdigest()}
(ROOT/'engine-audit-results.json').write_text(json.dumps(r,indent=2),encoding='utf8');print(json.dumps(r,indent=2),flush=True)

