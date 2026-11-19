"""Independently authored MR11 inputs/oracle. Standard DNA operations only; no Neoantimon2 or observed-data access."""
from pathlib import Path
import json,hashlib,itertools,datetime,argparse
CODONS=dict(zip((''.join(x) for x in itertools.product('TCAG',repeat=3)), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
ENC={'M':'ATG','A':'GCC','V':'GTG','K':'AAG','L':'CTG','P':'CCA','D':'GAC','N':'AAT','Q':'CAG','E':'GAA','G':'GGC','H':'CAT','I':'ATC','R':'CGC','S':'TCC','T':'ACC','F':'TTC','C':'TGT','Y':'TAC','W':'TGG','U':'TGA'}
RC=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(RC)[::-1]
def write(p,s):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s)
def js(p,x):write(p,json.dumps(x,indent=2)+'\n')
def fasta(p,rows):
 text='';fai=''
 for name,seq in rows:
  header='>'+name+'\n';offset=len(text.encode())+len(header.encode());text+=header+seq+'\n';fai+=f'{name}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq)+1}\n'
 write(p,text);return fai

def translate(dna,sec_positions,origin=None):
 if origin is None:origin=list(range(1,len(dna)+1))
 protein='';trace=[];used=[]
 for start in range(0,len(dna)-2,3):
  codon=dna[start:start+3];old=origin[start:start+3];a=CODONS[codon]
  # Explicit Sec recoding only if the original designated codon is intact/in-frame and still TGA.
  if codon=='TGA' and any(old==list(range((p-1)*3+1,p*3+1)) for p in sec_positions):a='U';used += [p for p in sec_positions if old==list(range((p-1)*3+1,p*3+1))]
  trace.append({'aa_position':len(protein)+1,'codon':codon,'residue':a,'original_tx_positions':old,'mutant_tx_start':start+1})
  if a=='*':return {'protein':protein,'stop_aa':len(protein)+1,'stop_nt':start+1,'stop_origin':old,'trace':trace,'sec_used':used}
  protein+=a
 return {'protein':protein,'stop_aa':None,'stop_nt':None,'stop_origin':None,'trace':trace,'sec_used':used,'terminal_residual_nt':len(dna)%3}

def splice(genome,exons,strand):
 coords=[x for a,b in exons for x in range(a,b+1)]
 if strand=='-':coords=coords[::-1]
 dna=''.join(genome[x-1] for x in coords)
 if strand=='-':dna=dna.translate(RC)
 return dna,coords

def edit(dna,coords,events,strand):
 edits=[]
 for e in events:
  inds=[coords.index(g) for g in range(e['pos'],e['pos']+len(e['ref']))];start=min(inds);old=e['ref'] if strand=='+' else rc(e['ref']);new=e['alt'] if strand=='+' else rc(e['alt']);assert dna[start:start+len(old)]==old
  edits.append((start,old,new,e['label']))
 mapping=list(range(1,len(dna)+1))
 for start,old,new,label in sorted(edits,reverse=True):
  keep=0
  while keep<min(len(old),len(new)) and old[keep]==new[keep]:keep+=1
  suffix=0
  while suffix<min(len(old),len(new))-keep and old[-suffix-1]==new[-suffix-1]:suffix+=1
  oldmap=mapping[start:start+len(old)]
  newmap=oldmap if len(old)==len(new) else oldmap[:keep]+[None]*(len(new)-keep-suffix)+(oldmap[len(old)-suffix:] if suffix else [])
  dna=dna[:start]+new+dna[start+len(old):];mapping=mapping[:start]+newmap+mapping[start+len(old):]
 return dna,mapping,edits

def make(root):
 assert not (root/'inputs').exists();cases=[];oracle={};lengths=[8,9,10,15];base='MAVKLPDNQEGHIRSTFCYWASDGLPK';coding=''.join(ENC[a] for a in base);intr='GTACGTTTCCAG';ids=['M01','M02','M03','M04','T01','T02','T03','T04','S01','S02','S03']
 for cid in ids:
  secgroup=cid.startswith('S');dup=cid in ['M01','M02'];txs=['MRSEC1.1'] if secgroup else (['MRP1.1','MRP2.1'] if dup else ['MRP1.1']);strand='-' if cid=='M04' else '+'
  prot=base if not secgroup else base[:9]+'U'+base[10:];cod=''.join(ENC[a] for a in prot)
  transcript=cod+'TAA'+('GCTGGCA CCTAA'.replace(' ','') if cid.startswith('T') else 'GGGTAACTAACTAACC C'.replace(' ',''))
  # T extension encodes lost-stop Q, A, G, T, then the next TAA.
  if cid.startswith('T'):transcript=cod+'TAA'+'GCT'+'GGC'+'ACC'+'TAA'+'CCC'
  if strand=='+':genome=transcript[:18]+intr+transcript[18:];exons=[[1,18],[31,len(genome)]]
  else:genome=rc(transcript[:18]+intr+transcript[18:]);exons=[[1,len(transcript)-18],[len(transcript)-18+13,len(genome)]]
  tx,coords=splice(genome,exons,strand);assert tx==transcript
  normal_sec=[10] if cid in ['S01','S03'] else [];normal=translate(tx,normal_sec)
  cds_len=27 if cid=='S02' else len(cod);stop_positions=coords[cds_len:cds_len+3];assert len(stop_positions)==3
  # Reference protein/coding annotations deliberately follow ordinary TGA termination in S02; same DNA, different explicit Sec/ORF declaration.
  events=[]
  def snv(nt,new,label,gt='0|1',ps='901'):
   pos=coords[nt-1];ref=genome[pos-1];alt=new if strand=='+' else rc(new);assert ref!=alt;return {'label':label,'contig':'chr1','pos':pos,'ref':ref,'alt':alt,'GT':gt,'PS':ps}
  if cid in ['M01','M04']:events=[snv(25,'G','coding-SNV')];interval=(9,9)
  elif cid=='M02':
   nt=27;g=coords[nt-1];assert coords[nt+2]==g+3;events=[{'label':'codon-deletion','contig':'chr1','pos':g,'ref':genome[g-1:g+3],'alt':genome[g-1],'GT':'0|1','PS':'901'}];interval=(10,9)
  elif cid=='M03':
   g=coords[23];events=[{'label':'frame-insertion','contig':'chr1','pos':g,'ref':genome[g-1],'alt':genome[g-1]+'C','GT':'0|1','PS':'901'}];interval=(9,None)
  elif cid.startswith('T'):
   stop=snv(len(cod)+1,'C','stop-loss');events=[stop]
   if cid!='T01':events.append(snv(len(cod)+4,'A','downstream-substitution','1|0' if cid=='T03' else '0|1'))
   if cid=='T04':
    for e in events:e['GT']='0/1';e['PS']='.'
   interval=(len(prot)+1,None)
  else:events=[snv(13,'G','pre-Sec-SNV')];interval=(5,5)
  for e in events:assert genome[e['pos']-1:e['pos']-1+len(e['ref'])]==e['ref']
  applied=events if cid!='T03' else events[:1]
  mdna,mapping,edits=edit(tx,coords,applied,strand);mut=translate(mdna,normal_sec,mapping);a,b=interval;b=len(mut['protein']) if b is None else b
  windows=[]
  for length in lengths:
   for start in range(1,len(mut['protein'])-length+2):
    end=start+length-1;hit=(start<=b and end>=a) if a<=b else (start<a and end>=a)
    if hit:windows.append({'start':start,'end':end,'length':length,'peptide':mut['protein'][start-1:end]})
  if cid not in ['T04','S03']:assert windows and mut['stop_aa'] is not None,(cid,mut)
  ref=root/'inputs'/cid/'reference';ref.mkdir(parents=True);fai=fasta(ref/'genome.fa',[('chr1',genome)]);write(ref/'genome.fa.fai',fai)
  txrows=[];prows=[];gtf=['##description: independently authored MR11 synthetic reference; NOT real GRCh38/GENCODE data'];mane=[]
  def segments(positions):
   out=[]
   for p in sorted(positions):
    if out and p==out[-1][1]+1:out[-1][1]=p
    else:out.append([p,p])
   return out
  for i,t in enumerate(txs):
   p='P'+t;attr=f'gene_id "MRG.1"; transcript_id "{t}"; gene_type "protein_coding"; transcript_type "protein_coding"; protein_id "{p}"; gene_name "MR_SYNTHETIC";'
   lines=[]
   for feat,segs in [('transcript',[[1,len(genome)]]),('exon',exons),('CDS',segments(coords[:cds_len])),('start_codon',segments(coords[:3])),('stop_codon',segments(stop_positions))]:
    for l,r in segs:lines.append(f'chr1\tMR11_synthetic\t{feat}\t{l}\t{r}\t.\t{strand}\t'+('0' if feat in ['CDS','start_codon','stop_codon'] else '.')+'\t'+attr)
   gtf+=lines
   if i==0:mane=lines
   txrows.append((f'{t}|MRG.1|-|-|synthetic|synthetic|{len(tx)}',tx));prows.append((f'{p}|{t}|MRG.1|-|-|synthetic|synthetic|{len(normal["protein"])}',normal['protein']))
  txrows.append(('MRDUMMY.1|MRDUMG.1|-|-|catalog_only|synthetic|12','ATGTGACCCTAA'));prows.append(('PMRDUMMY.1|MRDUMMY.1|MRDUMG.1|-|-|catalog_only|synthetic|3','MUP'))
  fasta(ref/'transcripts.fa',txrows);fasta(ref/'coding.fa',txrows);fasta(ref/'proteins.fa',prows);write(ref/'annotation.gtf','\n'.join(gtf)+'\n');write(ref/'mane.gtf','\n'.join(mane)+'\n')
  header='#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n'
  write(ref/'mane.tsv',header+f'GeneID:MR\tMRG.1\tHGNC:MR\tMR_SYNTHETIC\tsynthetic\tNM_MR.1\tNP_MR.1\t{txs[0]}\tP{txs[0]}\tMANE Select\tNC_000001.11\t1\t{len(genome)}\t{strand}\n')
  sec='MRDUMMY.1\t2\n'
  if cid in ['S01','S03']:sec+='MRSEC1.1\t'+('999' if cid=='S03' else '10')+'\n'
  write(ref/'sec.tsv',sec)
  lines=['##fileformat=VCFv4.2',f'##contig=<ID=chr1,length={len(genome)}>','##FORMAT=<ID=GT,Number=1,Type=String,Description="Declared synthetic genotype">','##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Declared synthetic phase set">','#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMR_SAMPLE']
  for e in sorted(events,key=lambda x:x['pos']):lines.append(f'chr1\t{e["pos"]}\t{e["label"]}\t{e["ref"]}\t{e["alt"]}\t.\tPASS\t.\tGT:PS\t{e["GT"]}:{e["PS"]}')
  write(root/'inputs'/cid/'somatic.vcf','\n'.join(lines)+'\n')
  c={'case_id':cid,'reference':f'inputs/{cid}/reference','vcf':f'inputs/{cid}/somatic.vcf','transcripts':txs,'strand':strand,'events':events,'exons_genomic':exons,'lengths':lengths,'coding_length_nt':cds_len,'class_i_lengths':[8,9,10],'class_ii_catalog_length':15,'hla':['HLA-A*02:01','HLA-B*07:02'],'public_external_chain':cid.startswith('M')};cases.append(c)
  expectation={'case_id':cid,'condition':'NEW_MR11_SYNTHETIC_NOT_HISTORICAL18','reference_genomic_DNA':genome,'spliced_reference_DNA':tx,'transcript_genomic_map':coords,'normal':normal,'tumor_DNA':mdna,'tumor_original_nt_map':mapping,'tumor':mut,'causative_event_labels':[e['label'] for e in applied],'altered_interval_mutant_aa':[a,b],'transcripts':txs,'windows':windows,'sec_ids':[f'GENCODE_SEC:{txs[0]}:10'] if cid=='S01' else [],'expected_exit':1 if cid=='S03' else 0}
  if cid=='T04':expectation.update({'expected_unresolved':True,'reason':'EXTENSION_PHASE_UNRESOLVED','windows':[],'tumor':None,'tumor_DNA':None,'hypotheses_only_not_required_to_enumerate':[translate(edit(tx,coords,events[:1],strand)[0],[]),mut]})
  if cid=='S03':expectation.update({'expected_rejection_stage':'reference resource validation','expected_error_substring':'Sec/protein relationship mismatch','windows':[],'tumor':None})
  oracle[cid]=expectation
 js(root/'inputs/cases.json',cases);js(root/'expected/independent_expected.json',oracle)
 js(root/'expected/GENERATION_RECORD.json',{'generated_at':datetime.datetime.now().astimezone().isoformat(),'generator':'scripts/generate_independent.py','generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'method':'Independent nucleotide projection/splicing/editing, standard64codon table+explicit Sec, first-stop translation, predeclared mutation/junction/tail window geometry; no Core or observed data access','case_ids':ids,'new_study':'MR11_20261116_A','exact_comparison_no_sequence_tolerance':True})
 print(json.dumps({cid:{'normal':v['normal']['protein'],'mutant':v['tumor']['protein'] if v['tumor'] else None,'stop':v['tumor']['stop_aa'] if v['tumor'] else None,'windows_per_transcript':len(v['windows'])} for cid,v in oracle.items()},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);make(p.parse_args().root)
