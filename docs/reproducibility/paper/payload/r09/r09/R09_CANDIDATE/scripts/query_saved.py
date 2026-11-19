"""Common read-only R09 queries. Does not call Core, predictors or pVAC inference."""
from pathlib import Path
from decimal import Decimal
import argparse,csv,json,hashlib,os,sys,datetime

def table(p):
 with p.open()as f:return list(csv.DictReader(f,delimiter='\t'))
def num(v):
 if v is None or v in ['NA','X','']:return None
 d=Decimal(str(v));assert d.is_finite();return format(d.normalize(),'f')
def ar(x):return []if x is None else x if isinstance(x,list)else[x]
def key(r):return tuple(r[k]for k in ['sample_id','contig','pos','ref','alt','transcript','peptide_start','peptide_end','peptide_sequence','HLA'])
def make_context(r):return 'R09:'+':'.join([r['contig'],str(r['pos']),r['ref'],r['alt'],r['transcript'],str(r['peptide_start'])+'-'+str(r['peptide_end']),r['peptide_sequence']])
def neo(E,reread=False):
 N=E/'neo/observed/R09-resolved-contexts';ex=json.loads((N/('reopened_target_native_explain.json'if reread else'target_all_native_explain.json')).read_text());er=(json.loads((N/'reopened_public_views.json').read_text())['external']if reread else json.loads((N/'external_results.json').read_text()));ext={r['candidate_id']:r for r in er};rows=[]
 for c in ex['candidates']:
  r=c['result'];assert r['peptide_sequence']=='MKACDFGH';oarr=c['origins'];assert set(r['origin_ids'])=={o['origin_id']for o in oarr};ei={e['event_id']:e for e in c['events']};ev=ext[r['candidate_id']]['external_evidence'];assert len(ev)==1;ev=ev[0]
  for o in oarr:
   assert len(o['source_event_ids'])==1;e=ei[o['source_event_ids'][0]]
   rows.append(dict(sample_id=r['sample_id'],contig=e['contig'],pos=e['pos'],ref=e['ref'],alt=e['alt'],transcript=o['transcript_id']+'.'+o['transcript_version'],peptide_start=o['peptide_start_aa'],peptide_end=o['peptide_end_aa'],peptide_sequence=r['peptide_sequence'],HLA=r['mhc_key'],score=num(ev['score']),raw_score=ev['score'],metric='controlled_affinity_nm',unit=ev['unit'],direction=ev['score_direction'],status=ev['status'],source_locator=f'{N.name}/{"reopened_target_native_explain.json"if reread else"target_all_native_explain.json"}#candidates[{r["candidate_id"]}]/origins[{o["origin_id"]}];external_results:{r["candidate_id"]}',native_Origin_ID=o['origin_id'],native_Event_ID=e['event_id'],model=ev['model_version'],provider=ev['provider'],run_id=ev['run_id'],retrieval='NATIVE neo_explain/neo_external_results; shared normalization EXTERNAL_QUERY'))
 return rows

def pvac_table(E):
 rows=[]
 for i,r in enumerate(table(E/'pvac/R09.binding.filtered.tsv'),2):
  x=dict(sample_id=r['R09 Sample'],contig=r['Chromosome'],pos=int(r['Start'])+1,ref=r['Reference'],alt=r['Variant'],transcript=r['Transcript'],peptide_start=int(r['Sub-peptide Position']),peptide_end=int(r['Sub-peptide Position'])+int(r['Peptide Length'])-1,peptide_sequence=r['MT Epitope Seq'],HLA=r['HLA Allele'],score=num(r['MHCflurry MT IC50 Score']),raw_score=r['MHCflurry MT IC50 Score'],metric='controlled_affinity_nm',unit='nM',direction='lower_is_better',status='evaluated',source_locator=f'pvac/R09.binding.filtered.tsv:row{i}',native_Origin_ID=None,native_Event_ID=None,model='NOT_EXPOSED_IN_TESTED_STANDARD_FIELD',provider='MHCflurry schema carrier ONLY; CONTROL_NOT_INFERENCE',run_id='FILE_LEVEL_CONTROL_MANIFEST',retrieval='NATIVE binding_filter output; EXTERNAL_QUERY normalization')
  assert make_context(x)==r['R09 Canonical Context'];assert int(r['Stop'])==x['pos'];rows.append(x)
 return rows

def pvac_metrics(E,mode,sample):
 p=E/f'pvac/{mode}/R09.aggregated.metrics.json';data=json.loads(p.read_text());rows=[]
 for eventkey,event in data.items():
  if not isinstance(event,dict)or'good_binders'not in event:continue
  chrom,start,stop,ref,alt=eventkey.rsplit('-',4);assert int(stop)==int(start)+1
  for setkey,ts in event['good_binders'].items():
   for ann in ts['transcripts']:
    tx,gene,mut,protpos=ann.rsplit('-',3)
    for sequence,peptide in ts['peptides'].items():
     calls=peptide['individual_ic50_calls'];assert calls['algorithms']==['MHCflurry'];mutation_pos=int(peptide['mutation_position']);window=int(protpos)-mutation_pos+1
     for allele in peptide['hla_types']:
      vv=calls['MT'][allele];assert len(vv)==1
      rows.append(dict(sample_id=sample,contig=chrom,pos=int(start)+1,ref=ref,alt=alt,transcript=tx,peptide_start=window,peptide_end=window+len(sequence)-1,peptide_sequence=sequence,HLA=allele,score=num(vv[0]),raw_score=repr(vv[0]),metric='controlled_affinity_nm',unit='nM',direction='lower_is_better',status='evaluated',source_locator=f'pvac/{mode}/R09.aggregated.metrics.json#{eventkey}/good_binders/{setkey}/transcripts[{ann}]/peptides/{sequence}/individual_ic50_calls/MT/{allele}',native_Origin_ID=None,native_Event_ID=None,model='NOT_EXPOSED_IN_TESTED_STANDARD_FIELD',provider='MHCflurry schema carrier ONLY; CONTROL_NOT_INFERENCE',run_id='FILE_LEVEL_CONTROL_MANIFEST',retrieval='NATIVE aggregate metrics; EXTERNAL_QUERY across Event keys/Transcript Sets; sample from single-sample manifest; register derived from protein/mutation position for this SNV fixture'))
 return rows

def signatures(rows):return sorted([list(key(r))+[r['score'],r['unit'],r['direction'],r['status']]for r in rows])
def summary(rows):
 return dict(contexts=sorted({make_context(r)for r in rows}),Events=sorted({f'{r["contig"]}:{r["pos"]}:{r["ref"]}:{r["alt"]}'for r in rows}),unique_sequences=sorted({r['peptide_sequence']for r in rows}),pmhc=sorted({r['peptide_sequence']+'|'+r['HLA']for r in rows}),relations=len({key(r)for r in rows}))
def apply_policy(rows,excluded):return[r for r in rows if f'{r["contig"]}:{r["pos"]}:{r["ref"]}:{r["alt"]}'not in excluded]
def queries(rows,expect):
 byhla={h:sorted({make_context(r)for r in rows if r['HLA']==h})for h in expect['HLA_scores']};scores={h:sorted({r['score']for r in rows if r['HLA']==h})for h in expect['HLA_scores']};out={'Q1':{'contexts_per_HLA':byhla,'Events':summary(rows)['Events']},'Q2':{'unique_synthesized_sequences':summary(rows)['unique_sequences'],'HLA_scores':scores,'unit':'nM','direction':'lower_is_better','metric_scope':'CONTROL_NOT_INFERENCE'}}
 for q,name in [('Q3a','exclude_chr1'),('Q3b','exclude_chrB'),('Q4','exclude_all')]:out[q]=summary(apply_policy(rows,expect['exclusion_policies'][name]))
 return out

def main():
 a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--execution',required=True);a.add_argument('--out',required=True);a.add_argument('--reread',action='store_true');a.add_argument('--baseline');v=a.parse_args();P=Path(v.root);E=Path(v.execution);O=Path(v.out);O.mkdir(parents=True,exist_ok=False);expected=json.loads((P/'expected/query_expected.json').read_text());canonical=json.loads((P/'inputs/canonical.json').read_text());expected_rows=[]
 for c in canonical['contexts']:
  for m in canonical['measurements']:expected_rows.append(dict(**{k:c[k]for k in ['sample_id','contig','pos','ref','alt','transcript','peptide_start','peptide_end','peptide_sequence']},HLA=m['HLA'],score=m['score_string'],unit=m['unit'],direction=m['direction'],status=m['status']))
 base_queries=queries(expected_rows,expected);datasets={'Neoantimon2':neo(E,v.reread),'pVACtools_full_filtered':pvac_table(E),'pVACtools_default_metrics':pvac_metrics(E,'default',canonical['sample_id']),'pVACtools_affinity_metrics':pvac_metrics(E,'affinity-explicit',canonical['sample_id'])};checks=[];qrows=[];allrows=[];results={};sigs={}
 def ck(k,ok,detail=None):checks.append(dict(check=k,PASS=bool(ok),detail=detail))
 for name,rows in datasets.items():
  ck('exact_predeclared_relation_and_measurement_set',signatures(rows)==signatures(expected_rows),name);ck('no_duplicate_biological_context_rows',len(rows)==len({key(r)for r in rows})==6,name);results[name]=queries(rows,expected);sigs[name]=signatures(rows)
  for r in rows:allrows.append(dict(dataset=name,canonical_context_id=make_context(r),**r))
  for q,answer in results[name].items():
   passed=answer==base_queries[q];ck('query_expected_set',passed,dict(dataset=name,query=q));qrows.append(dict(query_id=q,tool=name,operation_class='EXTERNAL_QUERY',native_primitives='neo_explain/neo_external_results'if name=='Neoantimon2'else'binding_filter / generate_aggregated_report actual outputs',expected=json.dumps(base_queries[q],sort_keys=True),observed=json.dumps(answer,sort_keys=True),verdict='MATCH'if passed else'MISMATCH',evidence=';'.join(sorted({r['source_locator']for r in rows})),query_code='query_saved.py:queries/apply_policy; normalization'+(' via public API artifacts'if name=='Neoantimon2'else' from native CLI artifacts')))
 # Verify entire plain report rows survived the native filter, including missing values; not just expected counts.
 before=table(P/'inputs/R09.CONTROL_INPUT.all_epitopes.tsv');after=table(E/'pvac/R09.binding.filtered.tsv');ck('binding_output_exact_input_rows_values',before==after)
 for mode in ['default','affinity-explicit']:
  raw=json.loads((E/f'pvac/{mode}/R09.aggregated.metrics.json').read_text());groups={k:x for k,x in raw.items()if isinstance(x,dict)and'good_binders'in x};ck('aggregate_Event_and_transcript_counts',len(groups)==2 and sorted(sum(x['transcript_counts'])for x in groups.values())==[1,2],mode);ck('count_limit_and_inclusion_settings',raw['aggregate_inclusion_binding_threshold']==5000 and raw['aggregate_inclusion_count_limit']==(15 if mode=='default'else 1),mode)
  for ev in groups.values():
   ck('RNA_scalar_stays_missing',ev['RNA VAF']=='NA'and ev['gene_expr']=='NA',mode)
   for ts in ev['good_binders'].values():
    ck('transcript_expression_stays_missing',all(x=='NA'for x in ts['transcript_expr']),mode)
    for pep in ts['peptides'].values():ck('WT_and_percentiles_stay_missing',all(x=='NA'for x in pep['ic50s_WT']+pep['combined_percentiles_MT']+pep['binding_percentiles_MT']),mode)
 # Baseline/new-process identity is a file-based query comparison, not a fake pVAC project API.
 baseline=json.loads(Path(v.baseline).read_text())if v.baseline else None
 if baseline:
  for name in datasets:
   ok=sigs[name]==baseline['normalized_signatures'][name]and results[name]==baseline['query_answers'][name];ck('fresh_process_same_queries_and_relations',ok,name);qrows.append(dict(query_id='Q5',tool=name,operation_class='EXTERNAL_QUERY',native_primitives='new R process already reopened through neo_project/public APIs'if name=='Neoantimon2'else'new Python process reads saved native TSV/JSON; no pvac project-reopen API claimed',expected='same as first process',observed='same'if ok else'different',verdict='MATCH'if ok else'MISMATCH',evidence='FRESH_READ_RECEIPT.json; existing Neo separate_process_reopen.json',query_code='query_saved.py:--reread --baseline'))
 for fname,rows in [('NORMALIZED_RELATIONS.tsv',allrows),('QUERY_RESULTS.tsv',qrows)]:
  with(O/fname).open('w')as f:
   wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
 # py3.10/isolated pvac import is optional in first host-only pass and explicit for fresh pass.
 version=None
 if v.reread:
  import importlib.metadata as md
  import pvactools.lib.aggregate_all_epitopes as module
  version={'pvactools':md.version('pvactools'),'actual_module':module.__file__};assert version['pvactools']=='7.1.5'
 receipt={'status':'SUPPORTED'if all(x['PASS']for x in checks)else'FAILED','at':datetime.datetime.now().astimezone().isoformat(),'pid':os.getpid(),'Python':sys.version,'process_kind':'FRESH_SEPARATE_READER'if v.reread else'FIRST_EXTERNAL_QUERY','actual_pvac_import':version,'query_answers':results,'normalized_signatures':sigs,'checks':checks,'check_count':len(checks),'failures':[x for x in checks if not x['PASS']],'scope':'EXTERNAL_QUERY on genuine native outputs; no predictor inference, no candidate status injection or native counterfactual selection'};(O/'FRESH_READ_RECEIPT.json'if v.reread else O/'FIRST_QUERY_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k]for k in ['status','process_kind','pid','check_count','failures']},indent=2));sys.exit(0 if receipt['status']=='SUPPORTED'else 1)
if __name__=='__main__':main()
