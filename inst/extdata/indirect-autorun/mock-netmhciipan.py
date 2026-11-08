#!/usr/bin/env python3
"""NEOANTIMON2_TEST_MOCK: controlled stdout, not NetMHCIIpan or real scores."""
import json,sys,time
from pathlib import Path
args=sys.argv[1:]
def option(k): return args[args.index(k)+1]
if '-h' in args:
 print('# NetMHCIIpan version 3.2\n# NEOANTIMON2_TEST_MOCK\n-f FILE -a ALLELE -inptype TYPE -rdir ROOT -listMHC')
 sys.exit(0)
root=Path(option('-rdir')); mode=(root/'mode.txt').read_text().strip()
if '-listMHC' in args:
 print('DRB1_9999' if mode=='unsupported' else 'DRB1_0101\nDRB1_0301');sys.exit(0)
if mode=='failure':print('Controlled failure',file=sys.stderr);sys.exit(7)
if mode=='timeout':time.sleep(4)
peptides=Path(option('-f')).read_text().splitlines();allele=option('-a')
records=json.loads((root/'data/binding-map.json').read_text());by={r['peptide_sequence']:r for r in records}
print('# NetMHCIIpan version 3.2\n# NEOANTIMON2_TEST_MOCK; no real inference')
header='Seq Allele Peptide Identity Pos Core Core_Rel 1-log50k(aff) Affinity(nM) %Rank Exp_Bind BindingLevel'
if mode=='badcolumns':header=header.replace('Affinity(nM)','Affinity_EL')
if mode=='no-rank':header=header.replace(' %Rank','')
if mode=='pos-layout':header='Pos MHC Peptide Of Core Core_Rel Identity Score Affinity(nM) %Rank BindLevel'
print(header);print('-'*100)
for i,pep in enumerate(peptides):
 if mode=='truncated' and i==len(peptides)-1:continue
 if mode=='missing-self' and pep=='ACDEFGHIKAAAAAA':continue
 r=by[pep];offset=r['cores'][0]['offset'];core=r['cores'][0]['core_sequence'];aff=r['affinity_nM']
 if mode=='badcore':core='XXXXXXXXX'
 if mode=='bad-offset':offset=7
 if mode=='wrongallele':allele='DRB1_9999'
 if mode=='badnumber':aff='NA'
 if mode=='conflict' and i==0:
  print(i,allele,pep,'Sequence',1,pep[1:10],1,0.5,aff,0.5,9.999)
 if mode=='pos-layout':fields=[i,allele,pep,offset,core,1,'Sequence',0.5,aff,0.5]
 else:
  fields=[i,allele,pep,'Sequence',offset,core,1,0.5,aff]
  if mode!='no-rank':fields.append(0.5)
  fields.append(9.999)
 print(*fields)
print('-'*100)
