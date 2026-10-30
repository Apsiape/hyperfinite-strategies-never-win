"""Exact labelled-GIII compiler/checker. Python -u; standard library only.
Coefficients are Q(sqrt(2)); operators are 18x18 matrices of S_mu U_g S_nu*.
No finite-depth/truncation test is used to establish a losing identity.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import json, itertools, hashlib, time, ctypes
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parent
STARTED=time.monotonic()
ZERO=(F(0),F(0)); ONE=(F(1),F(0)); HALF=(F(1,2),F(0)); SQ=(F(0),F(1))
def co(x): return x if isinstance(x,tuple) else (F(x),F(0))
def ca(x,y): return (x[0]+y[0],x[1]+y[1])
def cm(x,y): return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def cn(x): return (-x[0],-x[1])
def pow2half(n):
    k,r=divmod(n,2)
    return (F(2)**k,F(0)) if not r else (F(0),F(2)**k)
def red(w):
    out=[]
    for c in w:
        if out and out[-1]==c: out.pop()
        else: out.append(c)
    return ''.join(out)
OUT={'a':(1,0),'b':(0,1),'c':(0,1)}
RES={'a':('c','c'),'b':('a','b'),'c':('b','a')}
@lru_cache(maxsize=10000)
def act(g,w):
    v=w; section=''
    for c in reversed(g):
        q=c; out=[]
        for x in v:
            out.append(OUT[q][x]); q=RES[q][x]
        v=tuple(out); section=red(q+section)
    return v,section
@lru_cache(maxsize=30000)
def mono(mu,g,nu,alpha,h,beta):
    if alpha[:len(nu)]==nu:
        v,k=act(g,alpha[len(nu):])
        return mu+v,red(k+h),beta
    if nu[:len(alpha)]==alpha:
        v,k=act(h[::-1],nu[len(alpha):])
        return mu,red(g+k[::-1]),beta+v
    return None

def clean(d): return {k:v for k,v in d.items() if v!=ZERO}
def add(*xs):
    d={}
    for x in xs:
        for k,v in x.items(): d[k]=ca(d.get(k,ZERO),v)
    return clean(d)
def scale(x,c):
    c=co(c)
    return clean({k:cm(v,c) for k,v in x.items()})
def sub(x,y): return add(x,scale(y,-1))
def star(x): return {(j,i,nu,g[::-1],mu):v for (i,j,mu,g,nu),v in x.items()}
def mul(x,y):
    d={}; right={}
    for (k,j,a,h,b),v in y.items(): right.setdefault(k,[]).append((j,a,h,b,v))
    for (i,k,mu,g,nu),v in x.items():
        for j,a,h,b,w in right.get(k,[]):
            m=mono(mu,g,nu,a,h,b)
            if m is not None:
                key=(i,j,*m); d[key]=ca(d.get(key,ZERO),cm(v,w))
    return clean(d)
def prod(xs):
    xs=list(xs)
    if not xs: return ID
    out=xs[0]
    for x in xs[1:]: out=mul(out,x)
    return out

def term(i,j,mu=(),g='',nu=(),c=1): return {(i,j,tuple(mu),g,tuple(nu)):co(c)}
def eflag(i,j): return add(*(term(2*(i%9)+a,2*(j%9)+a) for a in range(2)))
ID=add(*(eflag(j,j) for j in range(9)))
FLAGS=[eflag(j,j) for j in range(9)]
def scalar(mu=(),g='',nu=()): return add(*(term(i,i,mu,g,nu) for i in range(18)))
T=[scalar((x,)) for x in range(2)]
UG={q:scalar(g=q) for q in 'abc'}
AUX_X=add(*(term(2*j+a,2*j+1-a) for j in range(9) for a in range(2)))
AUX_Z=add(*(term(2*j+a,2*j+a,c=(-1)**a) for j in range(9) for a in range(2)))

def normal(x):
    if not x: return {}
    D=max(len(k[4]) for k in x); d={}
    for (i,j,mu,g,nu),v in x.items():
        for tail in itertools.product((0,1),repeat=D-len(nu)):
            out,sec=act(g,tail); key=(i,j,mu+out,sec,nu+tail)
            d[key]=ca(d.get(key,ZERO),v)
    return clean(d)
def assertzero(x,what):
    z=normal(x)
    if z: raise AssertionError((what,'nonzero',len(z),list(z.items())[:3]))
def eq(x,y,what): assertzero(sub(x,y),what)
def sig_half(x):
    return clean({k:cm(v,pow2half(len(k[4])-len(k[2]))) for k,v in x.items()})
def spectral(o,s): return scale(add(ID,scale(o,s)),HALF)

R=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
C=[[t for t in itertools.product((1,-1),repeat=3) if t[0]*t[1]*t[2]==eta] for eta in (1,1,-1)]
AQ=['R1','R2','R3','T0','T1']; BQ=['C1','C2','C3','U0','U1']
def alans(q): return R+['f'] if q[0]=='R' else [t for t in C[0] if t!=R[int(q[1])]]+[1,-1]
def boans(q): return C[int(q[1])-1]+['f'] if q[0]=='C' else [t for t in R if t!=R[int(q[1])]]+[1,-1]
def fms(q,r,a,b):
    if q[0]=='R' and r[0]=='C':
        if a=='f' or b=='f': return a==b
        return a[int(r[1])-1]==b[int(q[1])-1]
    if q[0]=='T' and r=='C1':
        return (b=='f' or b==R[int(q[1])]) if isinstance(a,int) else a==b
    if q=='R1' and r[0]=='U':
        return (a=='f' or a==R[int(r[1])]) if isinstance(b,int) else a==b
    if q[0]=='T' and r[0]=='U' and isinstance(a,int) and isinstance(b,int):
        return a*b==(-1)**(int(q[1])*int(r[1]))
    return True
# Exact 4x4 real magic square, entry O33 = Y tensor Y.
X=((0,1),(1,0)); Z=((1,0),(0,-1)); I=((1,0),(0,1))
def km(a,b): return tuple(tuple(a[i//2][j//2]*b[i%2][j%2] for j in range(4)) for i in range(4))
SQ4=[[km(X,I),km(I,X),km(X,X)],[km(I,Z),km(Z,I),km(Z,Z)],[km(X,Z),km(Z,X),tuple(tuple(-v for v in row) for row in km(((0,-1),(1,0)),((0,-1),(1,0))))]]

def emb4(j,m):
    return add(*(scale(mul(VS[j][a+1],star(VS[j][b+1])),m[a][b]) for a in range(4) for b in range(4) if m[a][b]))
VS=[]; LAB=[]; BASEA=[]; BASEB=[]; square=[]
for j in range(9):
    vs=[FLAGS[j]]
    for k in range(4):
        p=j+1+2*k
        vs.append(add(mul(eflag(p,j),star(T[0])),mul(eflag(p+1,j),star(T[1]))))
    VS.append(vs)
    def transport(l): return add(*(mul(mul(v,l),star(v)) for v in vs))
    lx=transport(add(mul(T[0],star(T[1])),mul(T[1],star(T[0]))))
    lz=transport(sub(mul(T[0],star(T[0])),mul(T[1],star(T[1]))))
    ly={q:transport(UG[q]) for q in 'abc'}
    a,b,c,d=lx,AUX_X,AUX_Z,lz
    osq=[[a,b,mul(a,b)],[c,d,mul(c,d)],[mul(a,c),mul(b,d),prod([a,c,b,d])]]
    obs={f'S{r}{c}':osq[r][c] for r in range(3) for c in range(3)}
    obs.update({f'Y{q}':ly[q] for q in 'abc'})
    LAB.append(obs); square.append(osq)
    ms=[[emb4(j,SQ4[r][c]) for c in range(3)] for r in range(3)]
    magic=sub(ID,FLAGS[j])
    def ray(ops,s):
        p=magic
        for o,v in zip(ops,s): p=scale(mul(p,add(magic,scale(o,v))),HALF)
        return p
    def signray(vec,s):
        # vec is list of coefficients in Q(sqrt2), in the four magic coordinates.
        out=scale(FLAGS[j],HALF)
        for a,va in enumerate(vec):
            out=add(out,scale(add(VS[j][a+1],star(VS[j][a+1])),cm(co(s),cm(HALF,va))))
            for b,vb in enumerate(vec): out=add(out,scale(mul(VS[j][a+1],star(VS[j][b+1])),cm(HALF,cm(va,vb))))
        return out
    am={}; bm={}
    for q in AQ:
        ps=[]
        for ans in alans(q):
            if q[0]=='R': p=FLAGS[j] if ans=='f' else ray(ms[int(q[1])-1],ans)
            elif isinstance(ans,int):
                x=int(q[1]); v=[ZERO]*4
                v[x]=(F(0),F(1,2)); v[x+2]=(F(0),F(1,2)); p=signray(v,ans)
            else: p=ray([ms[r][0] for r in range(3)],ans)
            ps.append(normal(p))
        am[q]=ps
    for q in BQ:
        ps=[]
        for ans in boans(q):
            if q[0]=='C': p=FLAGS[j] if ans=='f' else ray([ms[r][int(q[1])-1] for r in range(3)],ans)
            elif isinstance(ans,int):
                y=int(q[1]); v=[(F((-1)**(y*k),2),F(0)) for k in range(4)]; p=signray(v,ans)
            else: p=ray(ms[0],ans)
            ps.append(normal(p))
        bm[q]=ps
    BASEA.append(am); BASEB.append(bm)

A={}; B={}; AP={}; BP={}; pairs=[]
def question(side,name,answers,ps):
    target,pt=(A,AP) if side=='A' else (B,BP)
    assert name not in target
    target[name]=answers; pt[name]=ps

def pair(a,b,rule,ok):
    bad=[[ia,ib] for ia,aa in enumerate(A[a]) for ib,bb in enumerate(B[b]) if not ok(aa,bb)]
    pairs.append({'alice':a,'bob':b,'rule':rule,'forbidden':bad})

for j in range(9):
    for q in AQ: question('A',f'base:{j}:{q}',alans(q),BASEA[j][q])
    for q in BQ: question('B',f'base:{j}:{q}',boans(q),BASEB[j][q])
question('B','flag',list(range(9)),FLAGS)
for q in 'abc': question('B','global:'+q,[1,-1],[spectral(UG[q],s) for s in (1,-1)])
for j in range(9):
    for v,o in LAB[j].items(): question('B',f'label:{j}:{v}',[1,-1],[spectral(o,s) for s in (1,-1)])

def basepairs(j,q,name,getbase=lambda a:a):
    for r in BQ: pair(name,f'base:{j}:{r}','base',lambda a,b,q=q,r=r:fms(q,r,getbase(a),b))
    if q=='R1': pair(name,'flag','flag-link',lambda a,b:(getbase(a)=='f')==(b==j))
    if q=='R2': pair(name,'flag','corner-link',lambda a,b:(getbase(a)==R[0])==(b in ((j+1)%9,(j+2)%9)))
for j in range(9):
    for q in AQ: basepairs(j,q,f'base:{j}:{q}')
    contexts=[]
    for r in range(3): contexts.append((f'row{r}',[f'S{r}{c}' for c in range(3)],1))
    for c in range(3): contexts.append((f'col{c}',[f'S{r}{c}' for r in range(3)],1 if c<2 else -1))
    contexts += [(f'Y{q}',[f'Y{q}'],None) for q in 'abc']
    for ctx,vs,parity in contexts:
        signs=[]
        for ss in itertools.product((1,-1),repeat=len(vs)):
            p=1
            for s in ss:p*=s
            if parity is None or p==parity: signs.append(ss)
        cp=[prod([spectral(LAB[j][v],s) for v,s in zip(vs,ss)]) for ss in signs]
        for q in AQ:
            name=f'joint:{j}:{q}:{ctx}'
            answers=[(ans,ss) for ans in alans(q) for ss in signs]
            ps=[mul(p,l) for p in BASEA[j][q] for l in cp]
            question('A',name,answers,ps); basepairs(j,q,name,lambda a:a[0])
            for k,v in enumerate(vs): pair(name,f'label:{j}:{v}','label-match',lambda a,b,k=k:a[1][k]==b)
for q in 'abc':
    name='flag-global:'+q
    answers=list(itertools.product(range(9),(1,-1)))
    question('A',name,answers,[mul(FLAGS[f],spectral(UG[q],s)) for f,s in answers])
    pair(name,'flag','global-flag',lambda a,b:a[0]==b)
    pair(name,'global:'+q,'global-match',lambda a,b:a[1]==b)

# Conditional contexts: d=0 has just one answer; d=1 obeys the local equation.
def conditional(j,tag,block,vs,allowed):
    names=[]; ops=[]
    for v in vs:
        if v in ('a','b','c'): names.append('global:'+v); ops.append(UG[v])
        else: names.append(f'label:{j}:{v}'); ops.append(LAB[j][v])
    ss=[s for s in itertools.product((1,-1),repeat=len(vs)) if allowed(s)]
    answers=[(0,)]+[(1,)+s for s in ss]
    Cj=add(*(FLAGS[f] for f in block))
    ps=[sub(ID,Cj)]+[mul(Cj,prod([spectral(o,s) for o,s in zip(ops,signs)])) for signs in ss]
    name=f'relation:{j}:{tag}'
    question('A',name,answers,ps)
    pair(name,'flag','conditional-support',lambda a,b:(a[0]==1)==(b in block))
    for k,v in enumerate(names): pair(name,v,'conditional-match',lambda a,b,k=k:a[0]==0 or a[k+1]==b)
for j in range(9):
    p=(j+1)%9; r=(j+2)%9; cj=(p,r)
    conditional(j,'Zplus',(p,),['S11'],lambda s:s[0]==1)
    conditional(j,'Zminus',(r,),['S11'],lambda s:s[0]==-1)
    for q in 'abc': conditional(j,'Yflag'+q,(j,),['Y'+q,q],lambda s:s[0]==s[1])
    for q,sections in [('b',('a','b')),('c',('b','a'))]:
        for x,f in enumerate((p,r)): conditional(j,'inactive'+q+str(x),(f,),['Y'+q,sections[x]],lambda s:s[0]==s[1])
    conditional(j,'active',cj,['Ya','S00','c'],lambda s:s[0]*s[1]*s[2]==1)
    for q in 'ab': conditional(j,'commuteX'+q,cj,['S00',q],lambda s:True)

assert len(A)==561 and len(B)==157 and len(pairs)==3696,(len(A),len(B),len(pairs))
assert len({(p['alice'],p['bob']) for p in pairs})==len(pairs)
game={'name':'Bellaterra labelled game','version':1,'coefficient_field':'rational predicate; operator model Q(sqrt(2))','distribution':{'type':'uniform_on_supported_pairs','denominator':len(pairs)},'alice':A,'bob':B,'supported_pairs':pairs,'unsupported_pairs':'zero probability'}
with open(ROOT/'labelled-giii-game.json','w',encoding='utf8') as f: json.dump(game,f,separators=(',',':'))
print('TABLE',len(A),len(B),len(pairs),'max answers',max(map(len,list(A.values())+list(B.values()))),flush=True)

# Every nonzero model entry is exact. Use cyclic flag symmetry to check all instances
# against j=0 representatives; mapping shifts matrix indices and leaves aux/tree unchanged.
def shiftop(x,j):
    return {(2*((i//2-j)%9)+i%2,2*((k//2-j)%9)+k%2,mu,g,nu):v for (i,k,mu,g,nu),v in x.items()}
def anchor(name):
    parts=name.split(':')
    return int(parts[1]) if parts[0] in ('base','joint','label','relation') else 0

def fingerprint(x): return tuple(sorted(normal(x).items()))
checked_pvms=set(); checks={'pvm_templates':0,'povm_questions':0,'pair_templates':0,'pairs':0,'forbidden_events':0,'losing_products':0}
for side,questions,models in [('A',A,AP),('B',B,BP)]:
    for name,ans in questions.items():
        j=anchor(name); ps=[shiftop(x,j) for x in models[name]]
        sig=tuple(fingerprint(x) for x in ps)
        checks['povm_questions']+=1
        if sig not in checked_pvms:
            eq(add(*ps),ID,('sum',name))
            for k,p in enumerate(ps):
                eq(p,star(p),('self-adjoint',name,k)); eq(mul(p,p),p,('idempotent',name,k))
            checked_pvms.add(sig); checks['pvm_templates']+=1
    print('PVM DONE',side,checks,flush=True)
checked_pairs=set()
for num,p in enumerate(pairs):
    a,b=p['alice'],p['bob']; j=anchor(a)
    aps=[shiftop(x,j) for x in AP[a]]; bps=[sig_half(shiftop(x,j)) for x in BP[b]]
    sig=(tuple(fingerprint(x) for x in aps),tuple(fingerprint(x) for x in bps),tuple(map(tuple,p['forbidden'])))
    checks['pairs']+=1; checks['forbidden_events']+=len(p['forbidden'])
    if sig not in checked_pairs:
        for ia,ib in p['forbidden']:
            assertzero(mul(aps[ia],bps[ib]),('losing',a,b,ia,ib)); checks['losing_products']+=1
        checked_pairs.add(sig); checks['pair_templates']+=1
    if num%500==0: print('PAIR PROGRESS',num,checks,flush=True)
# Direct extraction identities, all nine copies, including active recursion.
S=[]
for x in range(2):
    chunks=[]
    for j in range(9):
        w=VS[j][1]
        z=mul(FLAGS[(j+1)%9],w) if x==0 else prod([LAB[j]['S00'],FLAGS[(j+2)%9],w])
        chunks.append(star(z))
    S.append(add(*chunks))
for x in range(2):
    eq(mul(star(S[x]),S[x]),ID,('isometry',x))
    for y in range(2):
        if x!=y: assertzero(mul(star(S[x]),S[y]),('orthogonal',x,y))
eq(add(*(mul(s,star(s)) for s in S)),ID,'Cuntz cover')
for q in 'abc':
    for x in range(2): eq(mul(UG[q],S[x]),mul(S[OUT[q][x]],UG[RES[q][x]]),('recursion',q,x))
from resource_usage import peaks
from types import SimpleNamespace
_peak,_commit=peaks()
assert _peak<1_000_000_000 and _commit<1_000_000_000
pm=SimpleNamespace(PeakWorkingSetSize=_peak,PeakPagefileUsage=_commit)
checks.update({'peak_working_set_bytes':pm.PeakWorkingSetSize,'peak_commit_bytes':pm.PeakPagefileUsage,'elapsed_seconds':round(time.monotonic()-STARTED,3)})
checks.update({'status':'PASS','exact':True,'finite_truncation_used':False,'arithmetic':'Q(sqrt(2)), Fraction','proof_rule':'A sigma_{-i/2}(B) = 0 in the self-similar Cuntz algebra','game_sha256':hashlib.sha256((ROOT/'labelled-giii-game.json').read_bytes()).hexdigest()})
(ROOT/'exact-check-results.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
print(json.dumps(checks,indent=2),flush=True)

