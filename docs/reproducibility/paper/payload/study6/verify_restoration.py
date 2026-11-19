#!/usr/bin/env python3
"""Private byte/binding verification only; public license pending. Never runs science."""
from pathlib import Path
import hashlib,json,tarfile
P=Path(__file__).resolve().parent;H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((P/'NEW_CANDIDATE_MANIFEST.json').read_text());e={x['path']:x for x in m['files']};actual={p.relative_to(P).as_posix() for p in P.rglob('*') if p.is_file()}
assert actual==set(e)|{'NEW_CANDIDATE_MANIFEST.json','SHA256SUMS'}
for n,x in e.items():assert H(P/n)==x['sha256'] and (P/n).stat().st_size==x['bytes'],n
for line in (P/'SHA256SUMS').read_text().splitlines():d,n=line.split('  ',1);assert H(P/n)==d,n
original=json.loads((P/'ORIGINAL_SOURCE503_EXPECTED_HASHES.json').read_text());missing={'docs/developer/bounded-repair.md'}
with tarfile.open(P/'study6/sources/patched.tar.gz') as t:
 source={}
 for m in t.getmembers():
  assert not Path(m.name).is_absolute() and '..' not in Path(m.name).parts and not m.issym() and not m.islnk()
  if m.isfile():source[m.name.split('/',1)[1]]=t.extractfile(m).read()
assert set(source)==set(original)-missing and all(hashlib.sha256(b).hexdigest()==original[n] for n,b in source.items())
assert sum(l.startswith(b'export(') for l in source['NAMESPACE'].splitlines())==20
freeze=json.loads((P/'study6/FREEZE_BEFORE_FRESH_V3.json').read_text())
assert all(H(P/'study6'/n)==d for n,d in freeze['files'].items())
cases=json.loads((P/'study6/inputs/cases.json').read_text());assert len(cases)==6
print(json.dumps({'status':'PASS_BYTES_AND_BINDING_ONLY','source_files':len(source),'API_count':20,'cases_declared':6,'frozen_V3_files_verified':len(freeze['files']),'science_executed':False},indent=2))
