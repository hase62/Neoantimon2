#!/usr/bin/env python3
"""NEOANTIMON2_TEST_MOCK: synthetic software fixture, not DTU code/models."""
import json,sys,time
from pathlib import Path
a=sys.argv[1:]
def arg(k):return a[a.index(k)+1]
if '-h' in a:
 print('# NetMHCIIpan version 4.3k\n# NEOANTIMON2_TEST_MOCK')
 print('[-f filename] input\n[-a allele] input\n[-inptype int] input\n[-rdir path] data\n[-list] 0 list\n[-BA] 0 include BA\n[-context] 0 context\n[-inv_all] 0 inversion\n[-tdir path] temp')
 sys.exit()
root=Path(arg('-rdir'));mode=(root/'mode.txt').read_text().strip()
if '-list' in a:
 print('# NetMHCIIpan version 4.3k\nDRB1_0101\nDRB1_0301');sys.exit()
if mode=='exit-zero-error':print('# NetMHCIIpan version 4.3k\nError. Cannot load model. Exit');sys.exit(0)
if mode=='failure':print('controlled failure',file=sys.stderr);sys.exit(7)
if mode=='timeout':time.sleep(4)
peps=Path(arg('-f')).read_text().splitlines();mhc=arg('-a')
rows=json.loads((root/'data/binding-map.json').read_text());by={r['peptide_sequence']:r for r in rows}
ba='-BA' in a and mode!='el-only'
fields=['Pos','MHC','Peptide','Of','Core','Core_Rel','Inverted','Identity','Score_EL','%Rank_EL','Exp_Bind']
if ba:fields+=['Score_BA','%Rank_BA','Affinity(nM)']
fields+=['BindLevel']
if mode=='badcolumns':fields.remove('Inverted')
print('# NetMHCIIpan version 4.3k\n# NEOANTIMON2_TEST_MOCK; no genuine prediction')
print(' '.join(fields));print('-'*140)
for i,p in enumerate(peps):
 if mode=='missing-self' and p=='ACDEFGHIKAAAAAA':continue
 if mode=='truncated' and i==len(peps)-1:continue
 r=by[p];of=r['cores'][0]['offset'];inv=1 if mode=='inverted' else 0
 core=(p[::-1] if inv else p)[of:of+9]
 c={'Pos':str(i+1),'MHC':mhc,'Peptide':p,'Of':str(of),'Core':core,'Core_Rel':'0.750','Inverted':str(inv),'Identity':'Synthetic','Score_EL':'0.123456','%Rank_EL':'12.34','Exp_Bind':'NA','Score_BA':'0.700000','%Rank_BA':'1.50','Affinity(nM)':str(r['affinity_nM'])}
 if mode=='bad-inverted':c['Inverted']='2'
 if mode=='badcore':c['Core']='XXXXXXXXX'
 if mode=='wrongallele':c['MHC']='DRB1_9999'
 if mode=='badnumber':c['Affinity(nM)']='NA'
 print(' '.join(c[k] for k in fields if k!='BindLevel'))
 if mode=='conflict' and i==0:
  c['Score_EL']='0.999999';print(' '.join(c[k] for k in fields if k!='BindLevel'))
print('-'*140);print('Number of strong binders: 0 Number of weak binders: 0');print('-'*140)
