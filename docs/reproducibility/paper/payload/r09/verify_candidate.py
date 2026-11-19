from pathlib import Path
import hashlib,json,zipfile,tarfile
P=Path(__file__).resolve().parent;H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((P/'NEW_MANIFEST.json').read_text());e={x['path']:x for x in m['files']};actual={p.relative_to(P).as_posix() for p in P.rglob('*') if p.is_file()};assert actual==set(e)|{'NEW_MANIFEST.json','SHA256SUMS'}
for n,x in e.items():assert H(P/n)==x['sha256'] and (P/n).stat().st_size==x['bytes'],n
for l in (P/'SHA256SUMS').read_text().splitlines():d,n=l.split('  ',1);assert H(P/n)==d,n
f=json.loads((P/'r09/R09_CANDIDATE/FREEZE_BEFORE_EXECUTION.json').read_text());assert all(H(P/'r09/R09_CANDIDATE'/n)==d for n,d in f['files'].items())
print(json.dumps({'status':'PASS_CANDIDATE_BYTES_ONLY','files':len(actual),'science_executed':False},indent=2))
