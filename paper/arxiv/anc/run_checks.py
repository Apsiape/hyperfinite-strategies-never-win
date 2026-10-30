"""Reproduce certificates in a separate directory; stop on any failure."""
from pathlib import Path
import argparse,os,sys,subprocess,json,hashlib,shutil
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output-dir',required=True)
args=parser.parse_args()
OUT=Path(args.output_dir).resolve()
if OUT==ROOT or ROOT in OUT.parents: raise ValueError('Use an output directory outside the ancillary source tree')
if not __debug__: raise RuntimeError('Assertions must be enabled')
os.environ['PYTHONDONTWRITEBYTECODE']='1'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
if (ROOT/'SHA256SUMS').exists():
    subprocess.run([sys.executable,'-u',str(ROOT/'verify_hashes.py')],check=True,cwd=ROOT)
OUT.mkdir(parents=True,exist_ok=True)
for name in ('compile_and_check.py','independent_predicate_check.py','independent_ideal_check.py','audit_engine.py','audit_forcing_table.py','resource_usage.py','labelled-giii-game.json','base-giii-game.json'):
    shutil.copyfile(ROOT/name,OUT/name)
for name in ('compile_and_check.py','independent_predicate_check.py','independent_ideal_check.py','audit_engine.py','audit_forcing_table.py'):
    print('RUN',name,flush=True)
    p=subprocess.run([sys.executable,'-u',str(OUT/name)],cwd=OUT,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=600)
    (OUT/(name.removesuffix('.py')+'.log')).write_text(p.stdout+p.stderr,encoding='utf-8',newline='\n')
    if p.returncode:
        print(p.stdout+p.stderr,flush=True)
        raise RuntimeError('Certificate failed: '+name)
    print('PASS',name,flush=True)
for name in ('exact-check-results.json','independent-predicate-results.json','independent-ideal-results.json','engine-audit-results.json','forcing-table-audit-results.json'):
    d=json.loads((OUT/name).read_text());assert d['status']=='PASS',name
summary={'status':'PASS','scripts':5,'assertions':True,'processes':'sequential','memory_limit_bytes':1_000_000_000,'game_sha256':hashlib.sha256((OUT/'labelled-giii-game.json').read_bytes()).hexdigest()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary),flush=True)
