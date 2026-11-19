"""Independent saved-output checks; no Core or predictor import. Never change expectations."""
from pathlib import Path
from collections import Counter,defaultdict
import json,csv,sys,hashlib,io
P=Path(__file__).resolve().parents[1];W=Path(sys.argv[1]).resolve();V=W/'verification';V.mkdir(exist_ok=True);cases={c['case_id']:c for c in json.loads((P/'inputs/cases.json').read_text())};expected=json.loads((P/'expected/closure_expected.json').read_text());checks=[];summ=[];decision=[];originrows=[];ctxrows=[];winrows=[];ranks=[];drows=[];refs=[]
def arr(x):return []if x is None else x if isinstance(x,list)else[x]
def ck(name,ok,case,detail=None):checks.append(dict(case_id=case,check=name,PASS=bool(ok),detail=detail))
def write(name,rows):
 if not rows:return
 with(V/name).open('w')as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows({k:json.dumps(v,separators=(',',':'))if isinstance(v,(list,dict))else v for k,v in r.items()}for r in rows)
def fields(e):return(e['contig'],e['pos'],e['ref'],e['alt'])
for cpath in sorted((W/'run/observed').glob('*/receipt.json')):
 cid=cpath.parent.name;c=cases[cid];ex=expected[cid];D=cpath.parent
 def rd(n):
  p=D/n;b=p.read_bytes();refs.append(dict(case_id=cid,file=str(p.relative_to(W)),bytes=len(b),sha256=hashlib.sha256(b).hexdigest()));return json.loads(b)
 receipt=rd('receipt.json');ck('completed_exit0',receipt['status']=='COMPLETED'and receipt['exit_status']==0,cid)
 if receipt['status']!='COMPLETED':summ.append(dict(case_id=cid,execution=receipt['status'],Events=None,Origins=None,sequences=None,pMHC=None,links=None));continue
 runtime=receipt['runtime'];ck('actual_loaded_package_isolated',Path(runtime['package_path'])==W/'library/Neoantimon2',cid,runtime['package_path']);ck('package_version',runtime['package_version']=='0.3.0.9000',cid)
 events=rd('events.json');os=rd('origins.json');cs=rd('native_results.json');ts=rd('translation_resolution.json');g=rd('public_graph.json');ci={x['candidate_id']:x for x in cs};oi={x['origin_id']:x for x in os};ei={x['event_id']:x for x in events};links={(l['candidate_id'],l['origin_id'])for l in arr(g['links'])};elinks={(x['candidate_id'],o)for x in cs for o in arr(x['origin_ids'])}
 ck('all_links_exact',links==elinks and len(links)==len(arr(g['links'])),cid);ck('all_linked_sequence_identity',all(a in ci and b in oi and ci[a]['peptide_sequence']==oi[b]['peptide_sequence']for a,b in links),cid);ck('all_causative_Events_exist',all(e in ei for o in os for e in arr(o['source_event_ids'])),cid)
 ck('original_files_unchanged_by_import',rd('source_file_preservation.json')['unchanged']and receipt['native_unchanged_by_external_import'],cid)
 reopen=rd('separate_process_reopen.json');ck('separate_process_comparisons',reopen['status']=='PASS'and all(reopen['comparisons'].values()),cid);ck('separate_PID',reopen['R_process_pid']!=runtime['R_process_pid'],cid);ck('separate_actual_package_path',Path(reopen['Core_package_path'])==W/'library/Neoantimon2',cid)
 summ.append(dict(case_id=cid,execution='COMPLETED',Events=len(events),Origins=len(os),sequences=len({o['peptide_sequence']for o in os}),pMHC=len(cs),links=len(links)))
 selected=set();audits={};ep={};pol=None
 if c.get('native'):
  payload=rd('actual_native_selection.json');ck('saved_selection_found',len(payload)==1,cid)
  pol=payload[0]['value']['normalized_policy']
  for s in payload[0]['value']['selections']:
   for item in s['selected_candidates']:selected.update(arr(item['candidate_ids']))
   for a in s['candidate_audit']:audits[a['candidate_id']]=a
  elig=rd('somatic_eligibility.json');ep={x['candidate_id']:x for x in elig['candidate_eligibility']};raw=rd('somatic_evidence.json');eve={e['event_id']:e for e in raw['event_evidence']};policy=elig['policy'];paths={x['candidate_id']:x['origin_paths']for x in raw['candidate_links']}
  ck('synthetic_counts_are_input_values',all(eve[e['event_id']]['tumor']['alt_reads']==c['dna_counts'][next(q['event_id']for q in c['events']if fields(q)==fields(e))]and eve[e['event_id']]['normal']['alt_reads']==0 for e in events),cid)
  targets=[x for x in cs if x['peptide_sequence']in c['targets']and x['hla_class']=='I'];target_origins=set()
  for x in targets:
   preds=[];per=[]
   for p in paths[x['candidate_id']]:
    o=oi[p['origin_id']];target_origins.add(o['origin_id']);vals=[eve[e]for e in arr(p['somatic_event_ids'])];hard=any(v['tumor']['alt_reads']<policy['tumor_min_alt_reads']or v['normal']['alt_reads']>policy['normal_max_alt_reads']for v in vals);condition=p.get('phase_condition',{});cond=condition.get('status')=='CONDITIONAL';pred='ineligible'if hard else'not_evaluated'if cond else'eligible';obs=next(z['status']for z in ep[x['candidate_id']]['origin_results']if z['origin_id']==o['origin_id']);ck('Origin_DNA_and_phase',obs==pred,cid,dict(origin=o['origin_id'],expected=pred,observed=obs));preds.append(pred)
    rr=dict(case_id=cid,candidate_id=x['candidate_id'],peptide=x['peptide_sequence'],HLA=x['mhc_key'],origin_id=o['origin_id'],transcript=o['transcript_id'],start=o['peptide_start_aa'],end=o['peptide_end_aa'],source_event_ids=arr(o['source_event_ids']),source_event_tuples=[fields(ei[e])for e in arr(o['source_event_ids'])],synthetic_DNA=[dict(event_id=v['event_id'],tumor_alt=v['tumor']['alt_reads'],normal_alt=v['normal']['alt_reads'])for v in vals],phase_state=o['phase_status'],phase_condition=condition,expected_Origin_status=pred,observed_Origin_status=obs);originrows.append(rr);per.append(rr)
   agg='eligible'if'eligible'in preds else'ineligible'if all(z=='ineligible'for z in preds)else'not_evaluated';actual=ep[x['candidate_id']]['eligibility_status'];ck('independent_Origin_aggregation',agg==actual==x['technical_eligibility'],cid)
   a=audits[x['candidate_id']];ck('actual_selection_membership',((x['candidate_id']in selected)==(a['status']=='selected')),cid)
   decision.append(dict(case_id=cid,candidate_id=x['candidate_id'],sequence=x['peptide_sequence'],HLA=x['mhc_key'],Origin_status_counts=dict(Counter(preds)),independent_aggregate=agg,actual_native_status=actual,actual_selected=x['candidate_id']in selected,selection_audit_status=a['status'],selection_reasons=a['reasons'],native_prediction_status=x['prediction_status'],Native_Grade=x['final_grade'],normalized_native_policy=pol))
   if cid.startswith('CLOSURE-A'):
    ck('A_predeclared_Origin_composition',dict(Counter(preds))=={k:v for k,v in ex['expected_origin_counts'].items()if v},cid);ck('A_expected_status',actual==ex['expected_candidate_status'],cid)
    for rr in per:drows.append(dict(case_id=cid,candidate_id=x['candidate_id'],sequence=x['peptide_sequence'],HLA=x['mhc_key'],all_Origins_native_status=actual,all_Origins_actual_selected=x['candidate_id']in selected,representative_origin_id=rr['origin_id'],representative_transcript=rr['transcript'],representative_support=rr['synthetic_DNA'],representative_phase=rr['phase_condition'],single_Origin_external_query_status=rr['observed_Origin_status'],single_Origin_external_query_would_select=rr['observed_Origin_status']=='eligible',information_lost=[q['origin_id']for q in per if q['origin_id']!=rr['origin_id']],query_scope='External representation/query contrast only; not a Core injection or another tool measurement'))
   elif cid.startswith('CLOSURE-B'):ck('B_expected_native_status',actual==ex['expected_target_status'],cid)
  ck('predeclared_target_selected_count',len(selected)==ex['expected_selected_pMHC'],cid,dict(expected=ex['expected_selected_pMHC'],actual=len(selected)))
  if cid.startswith('CLOSURE-A'):ck('A_target_HLA_and_origin_counts',{x['mhc_key']for x in targets}==set(ex['HLA'])and len(target_origins)==3,cid)
 if not cid.startswith('CLOSURE-A'):
  unknown=cid=='CLOSURE-B-unknown'
  hypotheses=(expected['CLOSURE-B-cis']['hypotheses']+expected['CLOSURE-B-trans']['hypotheses'])if unknown else ex['hypotheses']
  wanted={(tx,w['start'],w['end'],w['sequence'])for tx in ex['transcripts']for hh in hypotheses for w in hh['windows']};obs={(o['transcript_id']+'.'+o['transcript_version'],o['peptide_start_aa'],o['peptide_end_aa'],o['peptide_sequence'])for o in os}
  ck('windows_allowed_hypotheses'if unknown else'exact_complete_expected_windows',obs<=wanted if unknown else obs==wanted,cid,dict(expected_or_possible=len(wanted),observed=len(obs),missing=[list(x)for x in sorted(wanted-obs)],extra=[list(x)for x in sorted(obs-wanted)]))
  if not unknown:ck('exact_expected_Origin_count',len(os)==len(wanted),cid)
  observed_molecules=set()
  for t in ts:
   evtuple={fields(ei[e])for e in arr(t['somatic_event_ids'])};match=[hh for hh in hypotheses if {fields(e)for e in c['events']if e['event_id']in hh['causal_event_labels']}==evtuple];tr=t['translation'];ok=bool(match)and any(t['patient_normal_transcript']==hh['normal_RNA']and tr['patient_normal_protein']==hh['normal']['protein']and t['tumor_transcript']==hh['mutant_RNA']and tr['mutant_protein']==hh['mutant']['protein']and tr['mutant_stop_nt']==hh['mutant']['first_stop_nt']for hh in match);ck('protein_RNA_stop_causal_Events',ok,cid,t['transcript_id']);observed_molecules.add(tr.get('mutant_protein'))
   ctxrows.append(dict(case_id=cid,transcript=t['transcript_id'],normal_RNA=t.get('patient_normal_transcript'),mutant_RNA=t.get('tumor_transcript'),observed_protein=tr.get('mutant_protein'),first_stop_nt=tr.get('mutant_stop_nt'),source_event_ids=arr(t['somatic_event_ids']),source_event_tuples=sorted(evtuple),phase_state=t['phase_state'],phase_context=t['phase_context'],matches_independent_hypothesis=ok,context_assumption='Unknown conditional hypotheses, not an asserted phase'if unknown else'Explicit phase / fixed annotation'))
  if not unknown:ck('all_declared_proteins_observed',observed_molecules=={hh['mutant']['protein']for hh in hypotheses},cid)
  for o in os:
   evt={fields(ei[e])for e in arr(o['source_event_ids'])};match=[hh for hh in hypotheses if {fields(e)for e in c['events']if e['event_id']in hh['causal_event_labels']}==evt];ck('window_causal_Events',any(any((w['start'],w['end'],w['sequence'])==(o['peptide_start_aa'],o['peptide_end_aa'],o['peptide_sequence'])for w in hh['windows'])for hh in match),cid,o['origin_id'])
   winrows.append(dict(case_id=cid,origin_id=o['origin_id'],transcript=o['transcript_id']+'.'+o['transcript_version'],start=o['peptide_start_aa'],end=o['peptide_end_aa'],length=o['peptide_length'],sequence=o['peptide_sequence'],source_event_ids=arr(o['source_event_ids']),source_event_tuples=sorted(evt),phase_state=o['phase_status']))
 if cs:
  rawpred=list(csv.DictReader(io.StringIO((D/'controlled_predictions.tsv').read_text()),delimiter='\t'));idx={(q['peptide_sequence'],q['mhc_key']):q for q in rawpred};er=rd('external_results.json');groups=defaultdict(list)
  for x in er:
   q=idx.get((x['peptide_sequence'],x['mhc_key']));ev=arr(x['external_evidence']);ok=(not ev and x['external_score']is None and x['external_rank']is None)if q is None else(len(ev)==1 and ev[0]['score']==q['score']and ev[0]['status']==q['status']and x['external_score']==(float(q['score'])if q['status']=='evaluated'else None));ck('controlled_import_values_missingness',ok,cid,x['candidate_id']);
   if q:ck('controlled_model_provenance_and_Origins',len(ev)==1 and all(ev[0][f]==q[f]for f in ['provider','predictor_version','model_version','run_id','task','metric','metric_kind','metric_semantics','score_direction','unit'])and set(arr(ev[0]['origin_ids']))==set(arr(ci[x['candidate_id']]['origin_ids'])),cid,x['candidate_id'])
   ck('native_fields_preserved',all(x.get(k)==ci[x['candidate_id']].get(k)for k in ['technical_eligibility','prediction_status','final_grade','rank_within_sample_class','origin_ids']),cid,x['candidate_id'])
   if x['external_status']=='evaluated':groups[(x['sample_id'],x['hla_class'])].append(x)
   ranks.append(dict(case_id=cid,candidate_id=x['candidate_id'],peptide=x['peptide_sequence'],HLA=x['mhc_key'],raw_controlled_score=q['score']if q else None,external_status=x['external_status'],external_rank=x['external_rank'],native_status=x.get('technical_eligibility'),actual_native_selected=x['candidate_id']in selected,channel='CONTROLLED_SYNTHETIC_IMPORT_NOT_REAL_PREDICTOR'))
  for group in groups.values():ck('external_rank_independent_direction_and_ties',all(x['external_rank']==i for i,x in enumerate(sorted(group,key=lambda x:(-x['external_score'],x['candidate_id'])),1)),cid)
for name,rows in [('case_summary.tsv',summ),('native_target_decisions.tsv',decision),('Origin_support_decisions.tsv',originrows),('molecular_contexts.tsv',ctxrows),('all_Origin_windows.tsv',winrows),('numeric_rank_vs_native_selection.tsv',ranks),('D_all_representatives.tsv',drows),('raw_output_references.tsv',refs)]:write(name,rows)
status={k:('NOT_APPLICABLE_WITH_REASON'if not any(x['case_id'].startswith('CLOSURE-'+k)for x in checks)else'SUPPORTED'if all(x['PASS']for x in checks if x['case_id'].startswith('CLOSURE-'+k))else'FAILED')for k in ['A','B','C']};result={'component_status':status,'independent_checks':len(checks),'failures':[x for x in checks if not x['PASS']],'checks':checks,'scope':'New bounded synthetic executions; no historical counts replaced; external scores controlled, not model inference'};(V/'INDEPENDENT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'component_status':status,'checks':len(checks),'failures':result['failures']},indent=2));sys.exit(0 if not result['failures']else 1)
