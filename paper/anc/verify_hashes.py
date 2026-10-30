"""Verify every listed ancillary byte sequence before reproduction."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
manifest=(ROOT/'SHA256SUMS').read_bytes()
if b'\r' in manifest: raise RuntimeError('Manifest must use LF line endings')
count=0
for line in manifest.decode('utf-8').splitlines():
    digest,name=line.split('  ',1)
    path=(ROOT/name).resolve()
    if ROOT not in path.parents or not path.is_file(): raise RuntimeError('Invalid manifest entry')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=digest: raise RuntimeError('Hash mismatch: '+name)
    count+=1
print(json.dumps({'status':'PASS','files':count,'manifest':'SHA-256, LF'}))
