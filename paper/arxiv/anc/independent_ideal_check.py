"""Independent exact ideal checker on a random sample selected by the predicate audit.
Implements local word rewriting on generators, not the compiler's monomial product.
Builds base effects as outer products / 4x4 spectral projectors, not signray().
Only inputs: the game JSON and the previously recorded random sample.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import json, itertools, time, ctypes, hashlib
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from dataclasses import dataclass
ROOT = Path(__file__).resolve().parent
started = time.monotonic()

@dataclass(frozen=True, slots=True)
class Q:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    def __add__(self,other):
        other = other if isinstance(other,Q) else Q(Fraction(other))
        return Q(self.a+other.a,self.b+other.b)
    __radd__ = __add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,other): return self + (-other if isinstance(other,Q) else -Q(Fraction(other)))
    def __mul__(self,other):
        other = other if isinstance(other,Q) else Q(Fraction(other))
        return Q(self.a*other.a+2*self.b*other.b,self.a*other.b+self.b*other.a)
    __rmul__ = __mul__
    def __bool__(self): return bool(self.a or self.b)

ONE = Q(Fraction(1)); ZERO = Q(); HALF = Q(Fraction(1,2)); SQHALF = Q(Fraction(0),Fraction(1,2))
SECTIONS = {'a':('c','c'),'b':('a','b'),'c':('b','a')}
# Digits 0,1 are S0,S1; digits 2,3 are their adjoints. a,b,c are involutions.
@lru_cache(maxsize=60000)
def rewrite(word):
    while True:
        changed = False
        for k in range(len(word)-1):
            left,right = word[k:k+2]
            replacement = None
            if left in '23' and right in '01':
                if int(left)-2 != int(right): return None
                replacement = ''
            elif left in 'abc' and right in '01':
                x = int(right)
                replacement = str(x ^ (left == 'a')) + SECTIONS[left][x]
            elif left in '23' and right in 'abc':
                y = (int(left)-2) ^ (right == 'a')
                replacement = SECTIONS[right][y] + str(y+2)
            elif left == right and left in 'abc': replacement = ''
            if replacement is not None:
                word = word[:k] + replacement + word[k+2:]
                changed = True
                break
        if not changed: return word

adjmap = str.maketrans('0123','2301')
def adjword(w): return w.translate(adjmap)[::-1]
def addterm(out,key,value):
    value = out.get(key,ZERO) + value
    if value: out[key] = value
    elif key in out: del out[key]
def plus(*operators):
    out = {}
    for op in operators:
        for key,value in op.items(): addterm(out,key,value)
    return out
def times(op,c):
    return {key:value*c for key,value in op.items() if value*c}
def minus(a,b): return plus(a,times(b,-1))
def adj(op): return {(j,i,adjword(w)):c for (i,j,w),c in op.items()}
def multiply(a,b):
    result,right = {},{}
    for (j,k,w),c in b.items(): right.setdefault(j,[]).append((k,w,c))
    for (i,j,v),d in a.items():
        for k,w,c in right.get(j,()):
            reduced = rewrite(v+w)
            if reduced is not None: addterm(result,(i,k,reduced),d*c)
    return result
def chain(*args):
    out = ID
    for arg in args: out = multiply(out,arg)
    return out
def atom(i,j,w=''): return {(i,j,w):ONE}
def flagunit(j,k): return plus(*(atom(2*(j%9)+a,2*(k%9)+a) for a in range(2)))
ID = plus(*(atom(i,i) for i in range(18)))
F = [flagunit(j,j) for j in range(9)]
def diagonal(w): return plus(*(atom(i,i,w) for i in range(18)))
S = [diagonal('0'),diagonal('1')]
U = {q:diagonal(q) for q in 'abc'}

def expanded(op):
    if not op: return {}
    depth = max(sum(c in '23' for c in w) for i,j,w in op)
    result = {}
    for (i,j,w),c in op.items():
        cut = next((k for k,s in enumerate(w) if s in '23'),len(w))
        for bits in itertools.product('01',repeat=depth-(len(w)-cut)):
            tail = ''.join(bits)
            nw = rewrite(w[:cut]+tail+adjword(tail)+w[cut:])
            if nw is not None: addterm(result,(i,j,nw),c)
    return result

def iszero(op,why):
    rem = expanded(op)
    if rem: raise AssertionError((why,len(rem),list(rem.items())[:2]))
def equal(a,b,why): iszero(minus(a,b),why)
def mirror_shift(op):
    result = {}
    for (i,j,w),c in op.items():
        degree = sum(x in '01' for x in w)-sum(x in '23' for x in w)
        # 2^(-degree/2); integer and radical cases also handle negative powers.
        n,r = divmod(-degree,2)
        multiplier = Q(Fraction(2)**n,0) if r == 0 else Q(0,Fraction(2)**n)
        result[i,j,w] = c*multiplier
    return result

equal(plus(multiply(S[0],adj(S[0])),multiply(S[1],adj(S[1]))),ID,'Cuntz sum')
for x in range(2):
    for y in range(2): equal(multiply(adj(S[x]),S[y]),ID if x==y else {},'Cuntz orthogonality')
assert expanded(minus(U['a'],ID)), 'nonidentity group word erased'

# Ordinary real matrices are used to obtain the magic eigenprojections independently.
I2 = [[1,0],[0,1]]; X2 = [[0,1],[1,0]]; Z2 = [[1,0],[0,-1]]
J2 = [[0,-1],[1,0]]
def tensor(a,b): return [[a[i//2][j//2]*b[i%2][j%2] for j in range(4)] for i in range(4)]
magic = [[tensor(X2,I2),tensor(I2,X2),tensor(X2,X2)],
         [tensor(I2,Z2),tensor(Z2,I2),tensor(Z2,Z2)],
         [tensor(X2,Z2),tensor(Z2,X2),[[-v for v in row] for row in tensor(J2,J2)]]]
def matmul(a,b): return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def magic_projector(ops,ss):
    out = [[Fraction(i==j) for j in range(4)] for i in range(4)]
    for op,s in zip(ops,ss):
        out = matmul(out,[[Fraction(int(i==j)+s*op[i][j],2) for j in range(4)] for i in range(4)])
    return out

Vs, units, labels = {},{},{}
auxX = plus(*(atom(2*j+a,2*j+1-a) for j in range(9) for a in range(2)))
auxZ = plus(*(times(atom(2*j+a,2*j+a),(-1)**a) for j in range(9) for a in range(2)))

def prepare(j):
    if j in Vs: return
    v = [F[j]]
    for k in range(1,5):
        p = j+2*k-1
        v.append(plus(multiply(flagunit(p,j),adj(S[0])),multiply(flagunit(p+1,j),adj(S[1]))))
    Vs[j] = v
    units[j] = [[multiply(va,adj(vb)) for vb in v] for va in v]
    def transport(l): return plus(*(chain(vk,l,adj(vk)) for vk in v))
    x = transport(plus(multiply(S[0],adj(S[1])),multiply(S[1],adj(S[0]))))
    z = transport(minus(multiply(S[0],adj(S[0])),multiply(S[1],adj(S[1]))))
    grid = [[x,auxX,multiply(x,auxX)], [auxZ,z,multiply(auxZ,z)],
            [multiply(x,auxZ),multiply(auxX,z),chain(x,auxZ,auxX,z)]]
    labels[j] = {'S'+str(r)+str(c):grid[r][c] for r in range(3) for c in range(3)}
    labels[j].update({'Y'+q:transport(U[q]) for q in 'abc'})

def embed(j,m):
    prepare(j)
    return plus(*(times(units[j][a][b],m[a][b]) for a in range(5) for b in range(5) if m[a][b]))

def base_effect(j,q,outcome):
    prepare(j)
    if outcome == 'f': return F[j]
    if type(outcome) is int:
        t = int(q[1]); s = outcome
        if q[0] == 'T':
            ray = [SQHALF if k in (t,t+2) else ZERO for k in range(4)]
        else: ray = [Q(Fraction((-1)**(t*k),2)) for k in range(4)]
        vector = [SQHALF] + [s*SQHALF*r for r in ray]
        m = [[vector[a]*vector[b] for b in range(5)] for a in range(5)]
        return embed(j,m)
    ops = magic[int(q[1])-1] if q[0]=='R' else [magic[r][int(q[1])-1] for r in range(3)]
    if q[0] == 'T': ops = [magic[r][0] for r in range(3)]
    if q[0] == 'U': ops = magic[0]
    p = magic_projector(ops,outcome)
    m = [[ZERO]*5 for _ in range(5)]
    for a in range(4):
        for b in range(4): m[a+1][b+1] = Q(p[a][b])
    return embed(j,m)

def spectral(op,s): return times(plus(ID,times(op,s)),HALF)
def frozen(v): return tuple(map(frozen,v)) if isinstance(v,list) else v

game_path = ROOT/'labelled-giii-game.json'
game = json.loads(game_path.read_text(encoding='utf-8-sig'))
sample_info = json.loads((ROOT/'independent-predicate-results.json').read_text())
pairs = {(p['alice'],p['bob']):p for p in game['supported_pairs']}
model = {}

def get_model(side,name):
    key = (side,name)
    if key in model: return model[key]
    parts = name.split(':')
    ps = []
    for answer in game['alice' if side=='A' else 'bob'][name]:
        answer = frozen(answer)
        if name == 'flag': op = F[answer]
        elif parts[0] == 'global': op = spectral(U[parts[1]],answer)
        elif parts[0] == 'label':
            j = int(parts[1]); prepare(j)
            op = spectral(labels[j][parts[2]],answer)
        elif parts[0] == 'base': op = base_effect(int(parts[1]),parts[2],answer)
        elif parts[0] == 'flag-global': op = multiply(F[answer[0]],spectral(U[parts[1]],answer[1]))
        elif parts[0] == 'joint':
            j,q,ctx = int(parts[1]),parts[2],parts[3]; prepare(j)
            if ctx.startswith('row'): variables = ['S'+ctx[3]+str(c) for c in range(3)]
            elif ctx.startswith('col'): variables = ['S'+str(r)+ctx[3] for r in range(3)]
            else: variables = [ctx]
            op = base_effect(j,q,answer[0])
            for var,s in zip(variables,answer[1]): op = multiply(op,spectral(labels[j][var],s))
        else:
            j,tag = int(parts[1]),parts[2]; prepare(j)
            p,m = (j+1)%9,(j+2)%9
            if tag == 'Zplus': block,variables = (p,),['S11']
            elif tag == 'Zminus': block,variables = (m,),['S11']
            elif tag.startswith('Yflag'): block,variables = (j,),['Y'+tag[-1],tag[-1]]
            elif tag.startswith('inactive'):
                q,x = tag[-2],int(tag[-1]); block = ((p,m)[x],)
                variables = ['Y'+q, {'b':('a','b'),'c':('b','a')}[q][x]]
            elif tag == 'active': block,variables = (p,m),['Ya','S00','c']
            else: block,variables = (p,m),['S00',tag[-1]]
            support = plus(*(F[f] for f in block))
            if answer[0] == 0: op = minus(ID,support)
            else:
                op = support
                for var,s in zip(variables,answer[1:]):
                    observable = U[var] if var in ('a','b','c') else labels[j][var]
                    op = multiply(op,spectral(observable,s))
        ps.append(op)
    model[key] = ps
    return ps

sample = [tuple(x) for x in sample_info['sample_pairs']]
losing = 0
for index,(an,bn) in enumerate(sample):
    aa,bb = get_model('A',an),get_model('B',bn)
    bb_shifted = [mirror_shift(b) for b in bb]
    for a,b in pairs[an,bn]['forbidden']:
        iszero(multiply(aa[a],bb_shifted[b]),('losing',an,bn,a,b)); losing += 1
    if index%100 == 0: print('SAMPLE',index,'losing identities',losing,flush=True)
for key,ps in model.items():
    equal(plus(*ps),ID,('sum',key))
    for n,p in enumerate(ps):
        equal(p,adj(p),('self-adjoint',key,n))
        equal(multiply(p,p),p,('projection',key,n))
from resource_usage import peaks
from types import SimpleNamespace
_peak,_commit=peaks()
assert _peak<1_000_000_000 and _commit<1_000_000_000
mem=SimpleNamespace(peak=_peak,peak_commit=_commit)
result={'status':'PASS','independent_of_supplied_code':True,'sample_size':len(sample),'seed':sample_info['seed'],
        'losing_identities_checked':losing,'question_pvms_checked':len(model),'arithmetic':'exact Q(sqrt(2))',
        'method':'local generator rewriting; Cuntz expansion to common annihilation depth; separately constructed rank-one base effects',
        'finite_truncation':False,'elapsed_seconds':round(time.monotonic()-started,3),
        'peak_working_set_bytes':mem.peak,'peak_commit_bytes':mem.peak_commit,
        'game_sha256':hashlib.sha256(game_path.read_bytes()).hexdigest()}
(ROOT/'independent-ideal-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2),flush=True)
