"""Independent DNA edit / fixed-exon splicing / standard-code expectations. No Core import."""
from pathlib import Path
import json,hashlib,datetime,re
P=Path(__file__).resolve().parents[1]
AA={}
for aa,cs in {'F':'TTT TTC','L':'TTA TTG CTT CTC CTA CTG','I':'ATT ATC ATA','M':'ATG','V':'GTT GTC GTA GTG','S':'TCT TCC TCA TCG AGT AGC','P':'CCT CCC CCA CCG','T':'ACT ACC ACA ACG','A':'GCT GCC GCA GCG','Y':'TAT TAC','*':'TAA TAG TGA','H':'CAT CAC','Q':'CAA CAG','N':'AAT AAC','K':'AAA AAG','D':'GAT GAC','E':'GAA GAG','C':'TGT TGC','W':'TGG','R':'CGT CGC CGA CGG AGA AGG','G':'GGT GGC GGA GGG'}.items():
 for co in cs.split():AA[co]=aa
assert len(AA)==64
def fa(p):
 out={};key=None
 for l in p.read_text().splitlines():
  if l.startswith('>'):key=l[1:].split('|')[0];out[key]=''
  else:out[key]+=l.strip()
 return out
def protein(s):
 out=''
 for i in range(0,len(s)-2,3):
  a=AA[s[i:i+3]]
  if a=='*':return dict(protein=out,first_stop_nt=i+1,stop_codon=s[i:i+3])
  out+=a
 return dict(protein=out,first_stop_nt=None,stop_codon=None)
def txafter(g,coords,ev):
 # Each original genome base carries its coordinate; editing never invokes Core.
 chars={i:g[i-1]for i in coords}
 for e in ev:
  pos=e['pos'];ref=e['ref'];alt=e['alt'];assert g[pos-1:pos-1+len(ref)]==ref
  assert len(alt)<=len(ref)
  for j in range(len(ref)):
   if pos+j in chars:chars[pos+j]=alt[j]if j<len(alt)else''
 return ''.join(chars[i]for i in coords)
def windows(p,lengths,changed):
 return [dict(start=s,end=s+l-1,length=l,sequence=p[s-1:s+l-1])for l in lengths for s in range(1,len(p)-l+2)if any(s<=k<s+l for k in changed)]
cases=json.loads((P/'inputs/cases.json').read_text());exp={};boundary={}
for c in cases:
 if c['case_id'].startswith('CLOSURE-A'):
  exp[c['case_id']]={'target':'MKACDFGH','HLA':['HLA-A*02:01','HLA-B*07:02'],'expected_origin_counts':{'eligible':2 if c['expected_selected']else 0,'ineligible':0 if c['expected_selected']else 2,'not_evaluated':1},'expected_candidate_status':c['expected_status'],'expected_selected_pMHC':c['expected_selected'],'DNA_counts':c['dna_counts'],'native_policy':{'allowlist_sequences':c['targets'],'max_peptides':2,'fill_to_capacity':True},'clinical_scope':'Synthetic DNA support; full-context conditionality; not biological negativity'};continue
 genome=fa(P/c['reference']/'genome.fa')['chr1'];coords=list(range(1,18))+list(range(28,116));rna=''.join(genome[i-1]for i in coords);assert rna==fa(P/c['reference']/'transcripts.fa')['AUD1.1']
 normal=protein(rna);haps=[]
 if 'CLOSURE-B' in c['case_id']:
  modes=[c['events']]if c['case_id'].endswith('-cis')else[[e]for e in c['events']]
  for events in modes:
   mt=txafter(genome,coords,events);tr=protein(mt);haps.append(dict(causal_event_labels=[e['event_id']for e in events],normal_RNA=rna,mutant_RNA=mt,normal=normal,mutant=tr,windows=windows(tr['protein'],[8,15],[2])))
  # Unknown expectations describe trans-compatible singleton hypotheses, not mandatory exhaustive output.
  exp[c['case_id']]={'observational_rule':'For unknown, inspect generated hypotheses and require conditional native gate; exhaustive cis/trans output NOT required'if c['case_id'].endswith('unknown')else'exact declared molecular haplotypes','hypotheses':haps,'hypothesis_semantics':'possible singleton contexts; not certain trans'if c['case_id'].endswith('unknown')else'declared phase','expected_target_status':c['expected_target_status'],'expected_selected_pMHC':c['expected_selected'],'native_policy':{'allowlist_sequences':c['targets'],'max_peptides':6,'fill_to_capacity':True},'transcripts':['AUD1.1','AUD2.1']}
 else:
  assert genome[16:18]=='TG'and genome[26:28]=='GT'and rna[16:18]=='TT';mt=txafter(genome,coords,c['events']);assert mt==rna[:17]+rna[18:]==rna[:16]+'T'+rna[18:]
  tr=protein(mt);assert tr['protein']=='MEACDLVILN';assert tr['first_stop_nt']==31
  h=dict(causal_event_labels=['e27_GT_G'],normal_RNA=rna,mutant_RNA=mt,normal=normal,mutant=tr,windows=windows(tr['protein'],c['lengths'],list(range(6,len(tr['protein'])+1))))
  exp[c['case_id']]={'hypotheses':[h],'transcripts':['AUD1.1','AUD2.1'],'expected_Origins':2*len(h['windows']),'assumption':'fixed exons1-17/28-115; deletion of genomic28; splice alterations not modeled'}
  boundary={'old_REF':'TT','actual_genome17_18':genome[16:18],'old_REF_valid':False,'tx_nt17_genome':coords[16],'tx_nt18_genome':coords[17],'valid_variant':c['events'][0],'actual_REF27_28':genome[26:28],'fixed_splicing_RNA_equivalent_to_old_edit':True,'repeat_ambiguity':'Deleting either of adjacent transcript TT yields same RNA; genomic deletion28 with intronic anchor27 is one valid distinct genomic representation, not unique recovery of historical mutation','expected':h}
assert exp['CLOSURE-B-cis']['hypotheses'][0]['mutant']['protein']!=exp['CLOSURE-B-trans']['hypotheses'][0]['mutant']['protein']
for name,data in [('closure_expected.json',exp),('boundary_representation.json',boundary)]:
 p=P/'expected'/name
 if p.exists():raise RuntimeError('Expectations already exist; do not overwrite')
 p.write_text(json.dumps(data,indent=2)+'\n')
# Timestamp/hash evidence before any Core load or run.
files={p.relative_to(P).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for folder in ['inputs','expected','scripts','sources']for p in(P/folder).rglob('*')if p.is_file()}
(P/'FREEZE_BEFORE_CORE.json').write_text(json.dumps({'frozen_at':datetime.datetime.now().astimezone().isoformat(),'Core_used_for_expectations':False,'files':files},indent=2)+'\n')
print('Frozen independent expectations and inputs before Core execution')
