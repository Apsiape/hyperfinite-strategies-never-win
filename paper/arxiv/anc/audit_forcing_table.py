"""Check that the emitted table contains every hypothesis used by the forcing proof.
Read-only inputs: emitted JSON plus the paper's independently supplied GIII JSON.
"""
import json,itertools,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
D=json.loads((ROOT/'labelled-giii-game.json').read_text())
P=json.loads((ROOT/'base-giii-game.json').read_text())
A,B=D['alice'],D['bob']; cells={(p['alice'],p['bob']):p for p in D['supported_pairs']}
def freeze(a):return tuple(freeze(v) for v in a) if isinstance(a,list) else a
def semantic_bad(a,b,p):return {(freeze(a[i]),freeze(b[j])) for i,j in p['forbidden']}
def bname(p):return 'flag' if p=='Lab' else 'base:'+p
base=0; inherited=0; matches=0; cond=0
for p in P['supported_pairs']:
    an='base:'+p['alice'];bn=bname(p['bob'])
    expected=semantic_bad(P['alice'][p['alice']],P['bob'][p['bob']],p)
    actual=semantic_bad(A[an],B[bn],cells[(an,bn)])
    assert actual==expected,(an,bn)
    base+=1
contexts={}
for r in range(3):contexts['row'+str(r)]=([f'S{r}{c}' for c in range(3)],1)
for c in range(3):contexts['col'+str(c)]=([f'S{r}{c}' for r in range(3)],1 if c<2 else -1)
for q in 'abc':contexts['Y'+q]=(['Y'+q],None)

def exactly(a,b,accept):
    global matches
    p=cells[(a,b)]
    bad={(i,j) for i,x in enumerate(A[a]) for j,y in enumerate(B[b]) if not accept(x,y)}
    assert bad=={tuple(t) for t in p['forbidden']},(a,b)
    matches+=1
for j in range(9):
    for ctx,(variables,parity) in contexts.items():
        for q in ('R1','R2','R3','T0','T1'):
            an=f'joint:{j}:{q}:{ctx}'
            signs=[]
            for s in itertools.product((1,-1),repeat=len(variables)):
                p=1
                for v in s:p*=v
                if parity is None or p==parity:signs.append(s)
            assert {freeze(x) for x in A[an]}=={(freeze(a),s) for a in P['alice'][f'{j}:{q}'] for s in signs}
            for pp in P['supported_pairs']:
                if pp['alice']!=f'{j}:{q}':continue
                bn=bname(pp['bob']); target=semantic_bad(P['alice'][pp['alice']],P['bob'][pp['bob']],pp)
                actual={(freeze(A[an][i][0]),freeze(A[an][i][1]),freeze(B[bn][k])) for i,k in cells[(an,bn)]['forbidden']}
                expected={(a,s,b) for a,b in target for s in signs}
                assert actual==expected,(an,bn)
                inherited+=1
            for k,v in enumerate(variables): exactly(an,f'label:{j}:{v}',lambda a,b,k=k:a[1][k]==b)
    # Every conditional context invoked in the proof, checked from its logical meaning.
    plus,minus=(j+1)%9,(j+2)%9
    specifications=[('Zplus',(plus,),['S11'],lambda s:s[0]==1),('Zminus',(minus,),['S11'],lambda s:s[0]==-1)]
    specifications += [('Yflag'+q,(j,),['Y'+q,q],lambda s:s[0]==s[1]) for q in 'abc']
    specifications += [('inactiveb0',(plus,),['Yb','a'],lambda s:s[0]==s[1]),('inactiveb1',(minus,),['Yb','b'],lambda s:s[0]==s[1]),('inactivec0',(plus,),['Yc','b'],lambda s:s[0]==s[1]),('inactivec1',(minus,),['Yc','a'],lambda s:s[0]==s[1]),('active',(plus,minus),['Ya','S00','c'],lambda s:s[0]*s[1]*s[2]==1),('commuteXa',(plus,minus),['S00','a'],lambda s:True),('commuteXb',(plus,minus),['S00','b'],lambda s:True)]
    for tag,block,variables,allowed in specifications:
        an=f'relation:{j}:{tag}'
        assert {freeze(a) for a in A[an]}=={(0,)}|{(1,)+s for s in itertools.product((1,-1),repeat=len(variables)) if allowed(s)}
        exactly(an,'flag',lambda a,b:(a[0]==1)==(b in block))
        for k,v in enumerate(variables): exactly(an,'global:'+v if v in 'abc' else f'label:{j}:{v}',lambda a,b,k=k:a[0]==0 or a[k+1]==b)
        cond+=1
for q in 'abc':
    an='flag-global:'+q
    assert {freeze(a) for a in A[an]}==set(itertools.product(range(9),(1,-1)))
    exactly(an,'flag',lambda a,b:a[0]==b);exactly(an,'global:'+q,lambda a,b:a[1]==b)
result={'status':'PASS','paper_base_pairs_preserved':base,'inherited_base_pairs_checked':inherited,'label_and_conditional_matches_checked':matches,'conditional_contexts_checked':cond,'all_supported_pairs_accounted_for':base+inherited+matches,'game_sha256':hashlib.sha256((ROOT/'labelled-giii-game.json').read_bytes()).hexdigest()}
assert result['all_supported_pairs_accounted_for']==3696
(ROOT/'forcing-table-audit-results.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2),flush=True)
