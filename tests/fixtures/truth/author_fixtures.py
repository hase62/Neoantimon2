"""Independent fixture authoring/reference hash implementation (Python stdlib only).
No Neoantimon2 R function is called. Biological truth below is hand specified.
Schema is used ONLY for mechanical null/empty field padding; never to infer truth.
Run from repo root only when deliberately reviewing a new fixture version.
"""
import json, hashlib, unicodedata, struct, pathlib, copy
root=pathlib.Path(__file__).resolve().parents[3]
c=json.loads((root/'inst/schema/contract.json').read_text())
out=pathlib.Path(__file__).resolve().parent

def canonical(x):
 def net(t,p):return t+str(len(p.encode('utf-8')))+':'+p
 if x is None:return 'N0:'
 if isinstance(x,bool):return net('B','1' if x else '0')
 if isinstance(x,int):return net('I',str(x))
 if isinstance(x,float):return net('D',struct.pack('>d',0.0 if x==0 else x).hex())
 if isinstance(x,str):return net('S',unicodedata.normalize('NFC',x))
 if isinstance(x,list):return net('A',''.join(canonical(v) for v in x))
 d={unicodedata.normalize('NFC',k):v for k,v in x.items()}
 return net('O',''.join(canonical(k)+canonical(d[k]) for k in sorted(d,key=lambda s:s.encode())))

def ID(kind,*args):return kind+':'+hashlib.sha256(canonical(['Neoantimon2','1.0.0',kind]+list(args)).encode()).hexdigest()
def blank(table,record=False):
 return {f['name']:([] if f['type']=='array' else {} if f['type']=='object' else None) for f in c['records' if record else 'tables'][table]['fields'] if f['required']}
def context(mode):
 ct=blank('analysis_context',True)
 ct.update(context_id='ctx-somatic' if mode=='somatic' else 'ctx-donor-recipient',mode=mode,source_subject_id='patient-1' if mode=='somatic' else 'donor-1',source_subject_role='tumor' if mode=='somatic' else 'donor',self_subject_id='patient-1' if mode=='somatic' else 'recipient-1',self_subject_role='normal' if mode=='somatic' else 'recipient',presenter_subject_id='patient-1' if mode=='somatic' else 'recipient-1',presenter_subject_role='patient' if mode=='somatic' else 'recipient',sample_bindings=['sample-1'],direction='within_patient' if mode=='somatic' else 'donor_to_recipient',scenario_id='scenario-1',missing_reasons={'organ':'NOT_PROVIDED','tissue':'NOT_PROVIDED','timepoint':'NOT_PROVIDED'})
 return ct

def snapshot(mode='somatic'):
 ct=context(mode);m=blank('run_manifest');m.update(schema_version='1.0.0',id_schema_version='1.0.0',analysis_id='synthetic-analysis',package_version='0.0.1',fixture_kind='synthetic',analysis_contexts=[ct],allo_policy=None,output_schema_versions={t:'1.0.0' for t in ['events','peptide_origins','neoantigen_results','multisample_summary']})
 if mode!='somatic':
  p=blank('allo_policy',True);p.update(policy_status='NOT_RUN',reason='T002 contract fixture; no allo analysis performed');m['allo_policy']=p
 return dict(run_manifest=m,events=[],peptide_origins=[],neoantigen_results=[],multisample_summary=[])

def shared(s):
 ct=s['run_manifest']['analysis_contexts'][0]
 return dict(schema_version='1.0.0',analysis_id='synthetic-analysis',analysis_context=ct['context_id'],patient_id=ct['self_subject_id'],sample_id='sample-1')
def event(s,label='event-A',alt_count=5):
 ct=s['run_manifest']['analysis_contexts'][0];e=blank('events');e.update(shared(s))
 facts={'namespace':'synthetic-hand-defined-v1','identity':{'label':label,'contig':'synthetic-contig','pos':10,'ref':'A','alt':'T'}}
 ek=ID('event_key','synthetic-assembly',facts)
 e.update(canonical_event=facts,genome_assembly='synthetic-assembly',event_key=ek,event_id=ID('event_id',ct['context_id'],ct['self_subject_id'],'sample-1',ek),event_type='SNV',event_role='somatic',event_detection_status='CALLED',reconstruction_status='NOT_RUN',candidate_generation_status='NOT_RUN',tumor_alt_count=alt_count,tumor_count_unit='fragment',tumor_assessment_status='EVALUATED',tumor_vaf=0.25,rna_assessment_status='NOT_PROVIDED')
 if ct['mode']!='somatic':e.update(event_type='ALLO_HLA_SOURCE',event_role='source',event_detection_status='SOURCE_PROVIDED',tumor_alt_count=None,tumor_vaf=None,tumor_count_unit=None,tumor_assessment_status='NOT_APPLICABLE')
 s['events'].append(e);return e

def origin(s,e,ligand,label='protein-A'):
 ct=s['run_manifest']['analysis_contexts'][0];o=blank('peptide_origins');o.update(shared(s))
 pk=ID('peptide_key',ligand);sp=ID('sample_peptide_id',ct['context_id'],ct['self_subject_id'],'sample-1',pk)
 antigen=ID('source_antigen_key',ct['source_subject_id'],'synthetic-allele-A','synthetic-checksum-'+label,'scenario-1') if ct['mode']!='somatic' else None
 ii=dict(transcript_id=label,transcript_version='1',reconstructed_sequence_checksum='synthetic-checksum-'+label,peptide_start_aa=1,peptide_end_aa=len(ligand),source_event_keys=[e['event_key']],phase_namespace=None,phase_set=None,haplotype_id=None,source_subject_id=ct['source_subject_id'],source_antigen_key=antigen,scenario_id='scenario-1')
 o.update(origin_id=ID('origin_id',sp,ii),sample_peptide_id=sp,peptide_key=pk,peptide_sequence=ligand,ligand_sequence=ligand,origin_type=ct['mode'],sequence_id=label,sequence_context_id=label+'-context',peptide_length=len(ligand),peptide_start_aa=1,peptide_end_aa=len(ligand),transcript_id=label,transcript_version='1',source_event_ids=[e['event_id']],origin_validity='VALID',origin_identity=ii,source_antigen_key=antigen)
 if ct['mode']!='somatic':o.update(source_protein_id=label,source_allele='synthetic-allele-A',source_locus='synthetic-locus',source_sequence_ref='synthetic://'+label,source_domain='extracellular')
 s['peptide_origins'].append(o);return o

def result(s,origins,mhc='UNSPECIFIED',cls='UNSPECIFIED'):
 ct=s['run_manifest']['analysis_contexts'][0];o=origins[0];r=blank('neoantigen_results');r.update(shared(s));lig=o['ligand_sequence']
 ck=ID('candidate_key','NCBITaxon:9606',cls,lig,mhc)
 cid=ID('candidate_id',ct['context_id'],'sample-1',ck)
 r.update(candidate_id=cid,candidate_key=ck,sample_peptide_id=o['sample_peptide_id'],peptide_key=o['peptide_key'],ligand_sequence=lig,peptide_sequence=lig,peptide_length=len(lig),n_origins=len(origins),origin_ids=[z['origin_id'] for z in origins],mhc_key=mhc,hla_class=cls,species='NCBITaxon:9606',prediction_status='NOT_RUN',grading_status='NOT_GRADED',final_grade='NOT_GRADED',grade_before_warnings='NOT_GRADED',not_evaluated_codes=['NOT_RUN'] if ct['mode']=='somatic' else ['ALLO_GRADING_NOT_IMPLEMENTED'],default_include=False)
 s['neoantigen_results'].append(r);return r

def assignment(s,r,core,start,label,primary=False):
 z=blank('epitope_assignment',True);ep=ID('epitope_key','NCBITaxon:9606','II',core,r['mhc_key'])
 z.update(assignment_id=label,candidate_id=r['candidate_id'],prediction_ref='synthetic-prediction',request_id='synthetic-request',predictor_version='hand-authored-no-predictor',model_checksum='synthetic-no-model',binding_core=core,binding_core_start=start,binding_core_length=len(core),binding_orientation='forward',position_mapping=list(range(start,start+len(core))),epitope_key=ep,status='EVALUATED',is_primary=primary,selection_policy_id='synthetic-primary-only')
 r['epitope_assignments_json'].append(z)
 s['run_manifest']['predictions']=[dict(prediction_ref='synthetic-prediction',request_id='synthetic-request',predictor_version='hand-authored-no-predictor',model_checksum='synthetic-no-model')]
 if primary:
  for n in ['binding_core','binding_core_start','binding_core_length','binding_orientation','epitope_key']:r[n]=z[n]
  r['primary_prediction_ref']='synthetic-prediction'
 return z

cases={}
s=snapshot();e=event(s);o=origin(s,e,'ACDEFGHIK');result(s,[o]);cases['simple_somatic']=s
s=snapshot();a1=event(s);a2=event(s,'event-B',1);o1=origin(s,a1,'ACDEFGHIK');o2=origin(s,a2,'ACDEFGHIK','protein-B');o1['gene_tpm']=1.0;o2['gene_tpm']=20.0;a2.update(rna_alt_count=9,rna_count_unit='read',rna_assessment_status='EVALUATED');r=result(s,[o1,o2]);r.update(grading_origin_id=o1['origin_id'],tumor_alt_count=5,gene_tpm=1.0);cases['multiple_origins']=s
s=snapshot();e=event(s)
for i,lig in enumerate(['AAACDEFGHIKLMMM','TTTACDEFGHIKLNNN']):
 o=origin(s,e,lig,'protein-'+str(i));r=result(s,[o],'HLA-DRA*01:01~HLA-DRB1*01:01','II');assignment(s,r,'ACDEFGHIK',3+i,'assignment-'+str(i),True)
cases['class_ii_shared_core']=s
s=snapshot();e=event(s);o=origin(s,e,'AAACDEFGHIKLMMM');r=result(s,[o],'HLA-DRA*01:01~HLA-DRB1*01:01','II');assignment(s,r,'ACDEFGHIK',3,'register-1');assignment(s,r,'CDEFGHIKL',4,'register-2');r['not_evaluated_codes'].append('AMBIGUOUS_CORE');cases['multiple_registers']=s
s=snapshot('allo_hla');e=event(s);o=origin(s,e,'AAACDEFGHIKLMMM');r=result(s,[o],'HLA-DRA*01:01~HLA-DRB1*01:01','II');cases['subject_roles']=s
s=snapshot('allo_hla');e=event(s);o1=origin(s,e,'AAACDEFGHIKLMMM','source-HLA-A');o2=origin(s,e,'AAACDEFGHIKLMMM','source-HLA-B');r=result(s,[o1,o2],'HLA-DRA*01:01~HLA-DRB1*01:01','II')
ev=blank('evidence',True);ev.update(evidence_id='synthetic-B-evidence',analysis_context=r['analysis_context'],target_type='source_antigen',target_id=o1['source_antigen_key'],evidence_type='structural_b_cell',assessment_status='NOT_RUN',reason='NOT_RUN',applicability='UNKNOWN');r['evidence_records_json']=[ev]
link=blank('linked_recognition',True);link.update(link_id='synthetic-link',analysis_context=r['analysis_context'],source_antigen_key=o1['source_antigen_key'],b_cell_evidence_ids=[ev['evidence_id']],t_cell_origin_ids=[o1['origin_id']],candidate_ids=[r['candidate_id']],relationship_status='NOT_EVALUATED',method='hand-defined-relationship');r['linked_recognition_json']=[link];cases['linked_antigens']=s
s=copy.deepcopy(cases['simple_somatic']);s['events'][0].update(rna_alt_count=0,rna_count_unit='fragment',rna_assessment_status='EVALUATED');r=s['neoantigen_results'][0]
for label,status,value,finding in [('missing','NOT_RUN',None,None),('zero','EVALUATED',0.0,'NEGATIVE')]:
 ev=blank('evidence',True);ev.update(evidence_id=label,analysis_context=r['analysis_context'],target_type='candidate',target_id=r['candidate_id'],evidence_type='source_antigen_expression',assessment_status=status,value=value,finding=finding,reason='NOT_RUN' if status=='NOT_RUN' else None,unit='TPM',metric_id='synthetic-expression',method='hand-authored',source_ref='synthetic-truth',applicability='APPLICABLE',details={'scope':'synthetic source only'});r['evidence_records_json'].append(ev)
cases['zero_vs_missing']=s
for name,s in cases.items():(out/(name+'.json')).write_text(json.dumps(s,indent=2)+'\n')
vectors=[]
for value in [None,[],{},'','é','e\u0301','a:b|c',False,True,0,0.0,-0.0,1.5,{'z':None,'a':'é'},['a',None,'']]:
 payload=canonical(value);vectors.append(dict(value=value,canonical=payload,sha256=hashlib.sha256(payload.encode()).hexdigest()))
(out/'canonical_vectors.json').write_text(json.dumps(vectors,ensure_ascii=False,indent=2)+'\n')
(out/'expected.json').write_text(json.dumps({'fixture_kind':'synthetic','truth_source':'Hand-authored biological relationships, independent Python hashlib ID oracle; no R production calls, no public data, no real prediction.', 'cases':{'simple_somatic':{'events':1,'origins':1,'candidates':1},'multiple_origins':{'events':2,'origins':2,'candidates':1},'class_ii_shared_core':{'events':1,'origins':2,'candidates':2,'epitopes':1},'multiple_registers':{'events':1,'origins':1,'candidates':1,'registers':2,'epitopes':2},'subject_roles':{'events':1,'origins':1,'candidates':1},'zero_vs_missing':{'events':1,'origins':1,'candidates':1},'linked_antigens':{'events':1,'origins':2,'candidates':1}}},indent=2)+'\n')

# Same molecular candidate and sample, distinct contexts; somatic + allo in one run.
s=copy.deepcopy(cases['simple_somatic'])
second=snapshot();second['run_manifest']['analysis_contexts'][0]['context_id']='ctx-somatic-second'
e=event(second);o=origin(second,e,'ACDEFGHIK');result(second,[o])
allo=copy.deepcopy(cases['subject_roles'])
for other in [second,allo]:
 s['run_manifest']['analysis_contexts']+=other['run_manifest']['analysis_contexts']
 for table in ['events','peptide_origins','neoantigen_results']:s[table]+=other[table]
s['run_manifest']['allo_policy']=allo['run_manifest']['allo_policy']
for ct in s['run_manifest']['analysis_contexts'][:2]:
 z=blank('multisample_summary');z.update(schema_version='1.0.0',analysis_id='synthetic-analysis',analysis_context=ct['context_id'],patient_id='patient-1',peptide_key=ID('peptide_key','ACDEFGHIK'),peptide_sequence='ACDEFGHIK',peptide_length=9,n_samples_total=1,n_samples_detected=1,n_samples_peptide_assessed=1,n_samples_not_assessed=0,samples_detected=['sample-1'])
 s['multisample_summary'].append(z)
(out/'multi_context.json').write_text(json.dumps(s,indent=2)+'\n')
# Literal ordered inputs are independent of the R implementation and schema field lists.
inputs=[('candidate_key',dict(species='NCBITaxon:9606',hla_class='II',ligand_sequence='AAACDEFGHIKLMMM',mhc_key='HLA-DRA*01:01~HLA-DRB1*01:01')),
 ('epitope_key',dict(species='NCBITaxon:9606',hla_class='II',binding_core='ACDEFGHIK',mhc_key='HLA-DRA*01:01~HLA-DRB1*01:01'))]
ck=ID('candidate_key',*inputs[0][1].values())
inputs += [('candidate_id',dict(context_id=ctx,sample_id='sample-1',candidate_key=ck)) for ctx in ['ctx-somatic','ctx-somatic-second','ctx-donor-recipient']]
(out/'identity_vectors.json').write_text(json.dumps([dict(kind=k,values=v,expected=ID(k,*v.values())) for k,v in inputs],indent=2)+'\n')
