"""Independent read-only raw-score, complete-window and provenance adjudication.
Uses Python standard library only. Does not execute Neoantimon2 or a predictor; independently checks saved protein windows.
"""
from pathlib import Path
import csv,json,hashlib,sys,re
from decimal import Decimal
W=Path(sys.argv[1]);base=W/'real-final';cache=Path(sys.argv[2]) if len(sys.argv)>2 else W/'reference-cache/3.66.0'
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source=read(cache/'source.json')
assert source['commit']=='b2055ea33ef464781b9e8ecab20b275b6704f878'
for f,m in source['files'].items():
 assert sha(cache/f)==m['sha256']
 assert m['url']==f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{source["commit"]}/{f}'
fasta={};name=None
for line in (cache/'hla_prot.fasta').read_text().splitlines():
 if line.startswith('>'):name=line.split()[1];fasta[name]=''
 elif name:fasta[name]+=line.strip()
status={r['Allele']:r for r in csv.DictReader(l for l in (cache/'Allele_status.txt').read_text().splitlines() if not l.startswith('#'))}
reports={};all_loci=set();expected_files={};fingerprints={}
for root in sorted(p for p in base.iterdir() if p.is_dir()):
 manifest=read(root/'run.json');plan=read(root/'plan.json');result=read(root/'result.json');ex=read(root/'explanation.json');refs=read(root/'reference.json');rr={r['allele']:r for r in refs['records']}
 assert manifest['status']=='COMPLETED' and manifest['real_predictor_executed']
 assert result['status']=='COMPLETED' and ex['result']==result
 assert not result['official_pirche_reproduction'] and result['clinical_risk_score'] is None
 for r in refs['records']:
  assert r['sequence_status']=='complete'
  support=json.loads(r['resolution_provenance'][r['resolution_provenance'].index('{'):])
  assert support['commit']==source['commit'] and support['source_alleles']
  for a in support['source_alleles']:
   assert fasta[a]==r['protein_sequence'];st=status[a]
   assert st['Partial']=='Full' and int(st['Start'])==1 and int(st['End'])==3*(len(fasta[a])+1)
   raw=r['allele'].removeprefix('HLA-');assert a==raw or a.startswith(raw+':')
 origins=ex['origins'];oi={o['origin_id']:o for o in origins};assert len(oi)==len(origins)
 expected=set()
 for role,scope in [('donor',result['source_scope']),('recipient',result['self_scope'])]:
  typed=manifest['request'][role]
  for t in typed:
   if t['locus'] not in scope or t['status']!='typed':continue
   for allele in set(t['alleles']):
    protein=rr[allele]['protein_sequence']
    for start in range(len(protein)-14):expected.add((role,t['locus'],allele,start+1,start+15,protein[start:start+15]))
 actual={(o['role'],o['locus'],o['allele'],o['start'],o['end'],o['peptide_sequence']) for o in origins}
 assert expected==actual and len(actual)==len(origins)
 for o in origins:
  assert rr[o['allele']]['protein_sequence'][o['start']-1:o['end']]==o['peptide_sequence']
  if o['role']=='donor':all_loci.add(o['locus'])
 pairs={p['presenter']:p for p in plan['presenters']}
 assert all(p['status']=='evaluated' for p in pairs.values())
 for pair,p in pairs.items():
  alpha,beta=pair.split('~');strip=lambda s:s.removeprefix('HLA-')
  vendor=(strip(beta).replace('*','_').replace(':','') if p['group']=='DR' else 'HLA-'+strip(alpha).replace('*','').replace(':','')+'-'+strip(beta).replace('*','').replace(':',''))
  assert vendor==p['vendor']
  assert all(any(t['status']=='typed' and c in t['alleles'] for t in manifest['request']['recipient']) for c in [alpha,beta])
 qlinks=list(csv.DictReader((root/'query_origin_links.tsv').open(),delimiter='\t'))
 assert len(qlinks)==len(origins)*len(pairs)
 assert {(l['origin_id'],l['presenter']) for l in qlinks}=={(o,p) for o in oi for p in pairs}
 for l in qlinks:assert oi[l['origin_id']]['peptide_sequence']==l['peptide_sequence'] and oi[l['origin_id']]['role']==l['role']
 queries=list(csv.DictReader((root/'queries.tsv').open(),delimiter='\t'))
 assert {(q['peptide_sequence'],q['presenter']) for q in queries}=={(o['peptide_sequence'],p) for o in origins for p in pairs}
 assert len(queries)==len({q['query_id'] for q in queries})
 parsed={};raw_rows=0;inverted={p:0 for p in pairs};raw_hashes={}
 for command in sorted((root/'raw').glob('*.command.json')):
  cmd=read(command);argv=cmd['argv'];assert cmd['exit_status']==0 and '-BA' in argv and '-context' not in argv
  vendor=argv[argv.index('-a')+1];pair=next(p for p,v in pairs.items() if v['vendor']==vendor)
  pepfile=command.with_name(command.name.replace('.command.json','.pep'));peps=set(pepfile.read_text().splitlines());observed=set()
  stdout=command.with_name(command.name.replace('.command.json','.stdout'));raw_hashes[stdout.name]=sha(stdout)
  lines=stdout.read_text().splitlines();assert '# NetMHCIIpan version 4.3k' in lines
  header=None
  for line in lines:
   f=line.split()
   if not f or f[0].startswith('#') or set(f[0])=={'-'}:continue
   if f[0]=='Pos':header=[c for c in f if c!='BindLevel'];continue
   if not f[0].isdigit():assert line.startswith('Number of strong binders:');continue
   assert header and len(f)>=len(header)
   d=dict(zip(header,f[:len(header)]));assert d['MHC']==vendor and d['Peptide'] in peps
   p=d['Peptide'];off=int(d['Of']);inv=int(d['Inverted']);assert inv in [0,1] and 0<=off<=6
   assert (p[::-1] if inv else p)[off:off+9]==d['Core']
   assert 0<=Decimal(d['Core_Rel'])<=1 and 0<=Decimal(d['Score_EL'])<=1
   key=(p,pair);assert key not in parsed;parsed[key]=d;observed.add(p);raw_rows+=1;inverted[pair]+=inv
  assert observed==peps
 assert set(parsed)=={(q['peptide_sequence'],q['presenter']) for q in queries}
 imported=read(root/'predictions.json');assert len(imported)==len(parsed)
 metrics={(m['peptide_sequence'],m['presenter']):m for m in csv.DictReader((root/'model_metrics.tsv').open(),delimiter='\t')}
 for r in imported:
  key=(r['peptide_sequence'],r['presenter']);d=parsed[key];assert r['raw_output']['columns']==d
  assert Decimal(str(r['affinity_nM']))==Decimal(d['Affinity(nM)']) and Decimal(str(r['percentile_rank']))==Decimal(d['%Rank_BA'])
  assert r['cores'][0]['orientation']==('reverse' if d['Inverted']=='1' else 'forward')
  assert r['raw_output']['forward_ligand_core_start0']==(6-int(d['Of']) if d['Inverted']=='1' else int(d['Of']))
  for target,col in [('Core','Core'),('Core_Rel','Core_Rel'),('Inverted','Inverted'),('affinity_nM','Affinity(nM)'),('rank_BA','%Rank_BA'),('score_BA','Score_BA'),('score_EL','Score_EL'),('rank_EL','%Rank_EL')]:assert metrics[key][target]==d[col]
 donor={o['peptide_sequence'] for o in origins if o['role']=='donor'};self_peps={o['peptide_sequence'] for o in origins if o['role']=='recipient'}
 key=lambda d,p:(d['Core'],int(d['Inverted']),p)
 units={};selfkeys={}
 for p,meta in pairs.items():
  threshold=Decimal(str(result['thresholds'][meta['group']]))
  selfkeys[p]={key(d,p) for (pep,pair),d in parsed.items() if pair==p and pep in self_peps and Decimal(d['Affinity(nM)'])<threshold}
  units[p]={key(d,p) for (pep,pair),d in parsed.items() if pair==p and pep in donor and Decimal(d['Affinity(nM)'])<threshold}-selfkeys[p]
  summary=next(s for s in result['presenter_summary'] if s['presenter']==p);assert summary['complete_count']==len(units[p])
 assert {(u['core'],int(u['orientation']=='reverse'),u['presenter']) for u in result['units']}==set().union(*units.values())
 candidates=result['candidates'];ci={c['candidate_id']:c for c in candidates}
 assert len(ci)==len(donor)*len(pairs)
 assert {(c['peptide_sequence'],c['presenter']) for c in candidates}=={(p,a) for p in donor for a in pairs}
 for c in candidates:
  p=c['presenter'];d=parsed[(c['peptide_sequence'],p)];threshold=Decimal(str(result['thresholds'][pairs[p]['group']]))
  expected_decision='excluded_nonbinder' if Decimal(d['Affinity(nM)'])>=threshold else 'self_match' if key(d,p) in selfkeys[p] else 'supported_nonself_within_scope'
  assert c['decision']==expected_decision and Decimal(str(c['affinity_nM']))==Decimal(d['Affinity(nM)'])
  assert set(c['origin_ids'])=={o['origin_id'] for o in origins if o['role']=='donor' and o['peptide_sequence']==c['peptide_sequence']}
 links=result['candidate_origin_links'];assert len(links)==sum(len(c['origin_ids']) for c in candidates)
 assert {(l['candidate_id'],l['origin_id']) for l in links}=={(c['candidate_id'],o) for c in candidates for o in c['origin_ids']}
 for l in links:assert oi[l['origin_id']]['role']=='donor' and ci[l['candidate_id']]['presenter']==l['presenter']
 for g in result['group_summary']:
  union=set().union(*(v for p,v in units.items() if pairs[p]['group']==g['group']));assert g['complete_count']==len(union) and g['clinical_risk_score'] is None
 for s in result['locus_summary']:
  selected={o['peptide_sequence'] for o in origins if o['role']=='donor' and o['locus']==s['locus']}
  expected_units={key(parsed[(p,s['presenter'])],s['presenter']) for p in selected}&units[s['presenter']]
  assert s['complete_count']==len(expected_units)
 for line in (root/'SHA256SUMS.txt').read_text().splitlines():
  h,f=line.split('  ',1);assert sha(root/f)==h
 fingerprints[root.name]={str(f.relative_to(root)):{'sha256':sha(f),'mtime_ns':f.stat().st_mtime_ns} for f in root.rglob('*') if f.is_file()}
 reports[root.name]={'status':'PASS','source_loci_with_origins':sorted({o['locus'] for o in origins if o['role']=='donor'}),'donor_origins':result['donor_origin_count'],'self_origins':result['self_origin_count'],'unique_queries':len(queries),'pmhc_candidates':len(candidates),'candidate_origin_links':len(links),'raw_prediction_rows':raw_rows,'inverted_rows':inverted,'group_summary':result['group_summary'],'raw_stdout_sha256':raw_hashes}
assert all_loci=={'A','B','C','DRB1','DRB3','DRB4','DRB5','DQA1','DQB1','DPA1','DPB1','MICA','MICB'}
assert any(n>0 for p,n in reports['DQ_DP']['inverted_rows'].items() if 'DPB1' in p)
(W/'independent-verification.json').write_text(json.dumps({'status':'PASS','all_13_source_loci_verified':True,'cases':reports,'evidence':'Authentic pretrained 4.3k on official IPD proteins, fictional declared typing; no clinical validation'},indent=2))
(W/'before-reopen.json').write_text(json.dumps(fingerprints,indent=2));print(json.dumps(reports,indent=2))
