"""Read-only comparator of frozen independent oracle and new observed records; no Core imports or expectation generation."""
from pathlib import Path
import json,csv,hashlib,argparse,collections,io

def arr(x):return [] if x is None else x if isinstance(x,list) else [x]
def load(p):return json.loads(p.read_text())
def main(C,W):
 E=load(C/'expected/independent_expected.json');cases=load(C/'inputs/cases.json');checks=[];summaries=[];proteinrows=[];winrows=[];rankrows=[]
 def ck(cid,name,ok,expected=None,observed=None):checks.append({'case_id':cid,'check':name,'PASS':bool(ok),'expected':expected,'observed':observed})
 def read(D,n,default=None):return load(D/n) if (D/n).exists() else default
 for c in cases:
  cid=c['case_id'];e=E[cid];D=W/'observed'/cid;rec=read(D,'receipt.json',{})
  ck(cid,'actual_R_receipt_present',bool(rec))
  if cid=='S03':
   ck(cid,'expected_nonzero_rejection',rec.get('exit_status',0)!=0,1,rec.get('exit_status'))
   ck(cid,'declared_Sec_validation_reason',e['expected_error_substring'] in rec.get('error',''),e['expected_error_substring'],rec.get('error'))
   summaries.append({'case_id':cid,'classification':'EXPECTED_REJECTION' if all(x['PASS'] for x in checks if x['case_id']==cid) else 'FAIL','exit':rec.get('exit_status')});continue
  ck(cid,'public_run_completed_exit0',rec.get('status')=='COMPLETED' and rec.get('exit_status')==0,'COMPLETED exit0',rec.get('error',rec.get('status')))
  if rec.get('status')!='COMPLETED':summaries.append({'case_id':cid,'classification':'FAIL','exit':rec.get('exit_status')});continue
  contexts=arr(read(D,'translation_resolution.json',[]));origins=arr(read(D,'origins.json',[]));native=arr(read(D,'native.json',[]));events=arr(read(D,'events.json',[]));event_index={x['event_id']:x for x in events}
  expected_events={(x['contig'],x['pos'],x['ref'],x['alt']):x['label'] for x in c['events']}
  actual_events={(x['contig'],x['pos'],x['ref'],x['alt']) for x in events};ck(cid,'exact_source_variant_identity',actual_events==set(expected_events),[list(x) for x in expected_events],[list(x) for x in actual_events])
  wanted_causes=set(e['causative_event_labels'])
  def causes(ids):return {expected_events.get((event_index[i]['contig'],event_index[i]['pos'],event_index[i]['ref'],event_index[i]['alt']),'UNDECLARED') for i in arr(ids) if i in event_index}
  if cid=='T04':
   ck(cid,'nonempty_unresolved_molecular_context',len(contexts)>0)
   ck(cid,'conservative_unknown_extension_no_protein',all(x.get('translation',{}).get('mutant_protein') is None for x in contexts))
   ck(cid,'frozen_unknown_phase_reason',any(x.get('translation',{}).get('reason')==e['reason'] for x in contexts),e['reason'],[x.get('translation',{}).get('reason') for x in contexts])
   ck(cid,'no_asserted_unknown_phase_peptide_or_candidate',not origins and not native,0,len(origins));summaries.append({'case_id':cid,'classification':'PASS' if all(x['PASS'] for x in checks if x['case_id']==cid) else 'FAIL','exit':0});continue
  ck(cid,'exact_transcript_context_count',len(contexts)==len(c['transcripts']),len(c['transcripts']),len(contexts))
  for x in contexts:
   tx=x['transcript_id'];t=x.get('translation',{});ck(cid,tx+':expected_transcript',tx in c['transcripts'])
   for field,exp in [('patient_normal_protein',e['normal']['protein']),('mutant_protein',e['tumor']['protein']),('mutant_stop_aa',e['tumor']['stop_aa']),('mutant_stop_nt',e['tumor']['stop_nt']),('mutant_stop_origin',e['tumor']['stop_origin'])]:ck(cid,tx+':'+field,t.get(field)==exp,exp,t.get(field))
   ck(cid,tx+':independent_spliced_DNA_edit',x.get('tumor_transcript')==e['tumor_DNA'],e['tumor_DNA'],x.get('tumor_transcript'))
   ck(cid,tx+':exon_order_and_genomic_coordinate_map',t.get('spec',{}).get('exon_coordinates')==e['transcript_genomic_map'],e['transcript_genomic_map'],t.get('spec',{}).get('exon_coordinates'))
   ck(cid,tx+':causative_event_set',causes(x.get('somatic_event_ids'))==wanted_causes,sorted(wanted_causes),sorted(causes(x.get('somatic_event_ids'))))
   ck(cid,tx+':annotated_Sec_ids',arr(t.get('spec',{}).get('translation_exception_ids'))==e['sec_ids'],e['sec_ids'],t.get('spec',{}).get('translation_exception_ids'))
   ck(cid,tx+':actually_used_Sec_ids',arr(t.get('sec_sites_used'))==e['sec_ids'],e['sec_ids'],t.get('sec_sites_used'))
   proteinrows.append({'case_id':cid,'transcript':tx,'expected_protein':e['tumor']['protein'],'observed_protein':t.get('mutant_protein'),'expected_stop_aa':e['tumor']['stop_aa'],'observed_stop_aa':t.get('mutant_stop_aa'),'causative_Event_labels':sorted(causes(x.get('somatic_event_ids')))})
  wanted={(tx,w['start'],w['end'],w['peptide']) for tx in c['transcripts'] for w in e['windows']};obs=set()
  for o in origins:
   tx=o['transcript_id'];version=o.get('transcript_version');tx=tx+'.'+str(version) if version is not None and not tx.endswith('.'+str(version)) else tx
   k=(tx,o['peptide_start_aa'],o['peptide_end_aa'],o['peptide_sequence']);obs.add(k);ck(cid,'Origin:'+o['origin_id']+':exact_window',k in wanted)
   ck(cid,'Origin:'+o['origin_id']+':causative_Event_links',causes(o.get('source_event_ids'))==wanted_causes,sorted(wanted_causes),sorted(causes(o.get('source_event_ids'))))
   winrows.append({'case_id':cid,'origin_id':o['origin_id'],'transcript':tx,'start':o['peptide_start_aa'],'end':o['peptide_end_aa'],'peptide':o['peptide_sequence'],'expected':k in wanted,'event_labels':sorted(causes(o.get('source_event_ids')))})
  ck(cid,'complete_window_set_no_missing_extra',obs==wanted,{'count':len(wanted),'missing':[list(x) for x in sorted(wanted-obs)]},{'count':len(obs),'extra':[list(x) for x in sorted(obs-wanted)]})
  ck(cid,'all_molecular_Origins_retained_no_duplicate_ids',len(origins)==len(wanted)==len({x['origin_id'] for x in origins}),len(wanted),len(origins))
  expected_pairs={(w['peptide'],hla) for w in e['windows'] if w['length'] in c['class_i_lengths'] for hla in c['hla']};actual_pairs={(x['peptide_sequence'],x['mhc_key']) for x in native};ck(cid,'complete_unique_sequence_HLA_candidate_pairs',actual_pairs==expected_pairs,len(expected_pairs),len(actual_pairs))
  graph=read(D,'graph.json',{});by_seq=collections.defaultdict(set)
  for o in origins:by_seq[o['peptide_sequence']].add(o['origin_id'])
  for x in native:ck(cid,'Candidate:'+x['candidate_id']+':all_Origin_links',set(arr(x['origin_ids']))==by_seq[x['peptide_sequence']],sorted(by_seq[x['peptide_sequence']]),x['origin_ids'])
  expected_links={(x['candidate_id'],o) for x in native for o in by_seq[x['peptide_sequence']]};actual_links={(x['candidate_id'],x['origin_id']) for x in graph.get('links',[])};ck(cid,'public_graph_complete_candidate_Origin_links',expected_links==actual_links,len(expected_links),len(actual_links))
  peptides=arr(read(D,'peptide_provenance.json',[]))
  for p in peptides:
   pv=p.get('provenance',{});ck(cid,'Peptide:'+p['peptide_id']+':Sec_provenance',arr(pv.get('annotated_translation_exception_ids'))==e['sec_ids'] and arr(pv.get('sec_sites_used'))==e['sec_ids'])
  if c['public_external_chain']:
   ck(cid,'explanation_public_API_recorded',(D/'native_explain.json').is_file())
   external=arr(read(D,'external.json',[]));qi=sorted(actual_pairs);controls={}
   for i,key in enumerate(qi,1):controls[key]=('NOT_PREDICTED',None) if i%7==0 else ('unsupported',None) if i%11==0 else ('failed',None) if i%13==0 else ('evaluated',(i%3)/2)
   ck(cid,'all_candidates_retained_in_external_projection',len(external)==len(native) and {x['candidate_id'] for x in external}=={x['candidate_id'] for x in native})
   ready=[]
   for x in external:
    status,score=controls.get((x['peptide_sequence'],x['mhc_key']),('UNKNOWN',None));ck(cid,x['candidate_id']+':control_status_and_score',x['external_status']==status and x.get('external_score')==score,[status,score],[x['external_status'],x.get('external_score')])
    ck(cid,x['candidate_id']+':external_keeps_all_Origins',set(arr(x['origin_ids']))==by_seq[x['peptide_sequence']])
    if status=='evaluated':ready.append(x)
    else:ck(cid,x['candidate_id']+':unscored_not_ranked',x.get('external_rank') is None)
   for rank,x in enumerate(sorted(ready,key=lambda x:(-x['external_score'],x['candidate_id'])),1):ck(cid,x['candidate_id']+':independent_rank_and_tie',x['external_rank']==rank,rank,x['external_rank']);rankrows.append({'case_id':cid,'candidate_id':x['candidate_id'],'peptide':x['peptide_sequence'],'HLA':x['mhc_key'],'score':x['external_score'],'expected_rank':rank,'actual_rank':x['external_rank']})
   ck(cid,'actual_missing_zero_and_tied_values_present',any(s[0]=='NOT_PREDICTED' for s in controls.values()) and any(s==('evaluated',0) for s in controls.values()) and len(ready)>3)
   ck(cid,'native_records_and_old_files_unchanged',rec.get('native_unchanged_by_import')==True and read(D,'native_file_invariance.json',{}).get('unchanged')==True)
   reader=read(D,'separate_reader.json',{});ck(cid,'actual_separate_R_reader_all_public_views',reader.get('status')=='PASS' and reader.get('R_process_pid')!=rec['runtime']['R_process_pid'] and bool(reader.get('comparisons')) and all(reader['comparisons'].values()))
  summaries.append({'case_id':cid,'classification':'PASS' if all(x['PASS'] for x in checks if x['case_id']==cid) else 'FAIL','exit':rec.get('exit_status'),'expected_Origins':len(wanted),'actual_Origins':len(origins),'expected_pMHC':len(expected_pairs),'actual_pMHC':len(native)})
 V=W/'verification';V.mkdir(exist_ok=True)
 def dump(name,rows):
  fields=list(dict.fromkeys(k for x in rows for k in x));f=io.StringIO();w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for row in rows:w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()})
  (V/name).write_text(f.getvalue())
 dump('CASE_SUMMARY.tsv',summaries);dump('PROTEIN_COMPARISON.tsv',proteinrows);dump('WINDOW_ORIGIN_COMPARISON.tsv',winrows);dump('RANK_COMPARISON.tsv',rankrows);dump('ALL_CHECKS.tsv',checks)
 result={'status':'PASS' if len(summaries)==11 and all(x['PASS'] for x in checks) else 'PARTIAL','new_study':'MR11_20261116_A','case_count':len(summaries),'check_count':len(checks),'case_summary':summaries,'failures':[x for x in checks if not x['PASS']],'checks':checks,'scope':'New11fixed-source synthetic conditions; NOT historical18/2004 or prior6/946 replay'};(V/'INDEPENDENT_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','case_count','check_count','case_summary','failures']},indent=2));return 0 if result['status']=='PASS' else 1
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(main(a.root.resolve(),a.out.resolve()))
