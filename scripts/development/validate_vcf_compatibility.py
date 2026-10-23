#!/usr/bin/env python3
"""Opt-in, manifest-pinned VCF component validation; never runs peptide stages."""
import argparse
import collections
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('neo_public_acquisition', REPO / 'scripts/data/download_public_data.py')
ACQ = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ACQ)
MANIFEST = REPO / 'scripts/data/vcf_compatibility_manifest.json'
CATEGORIES = ('SUPPORTED_PASS','EXPECTED_REJECTION','UNEXPECTED_ERROR','INCORRECT_INTERPRETATION','UNRESOLVED')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def vcf_header_facts(path):
    """Record all declared identities without loading a population panel's rows."""
    opener = gzip.open if str(path).endswith('.gz') else open
    version = None
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('##fileformat='):
                version = line.strip().split('=', 1)[1]
            elif line.startswith('#CHROM'):
                samples = line.rstrip('\r\n').split('\t')[9:]
                return dict(vcf_fileformat=version, declared_sample_count=len(samples),
                            declared_sample_ids=samples)
    raise ValueError('Verified source has no VCF column header: ' + str(path))


def command(args, output=None):
    if output:
        with Path(output).open('wb') as stream:
            subprocess.run([str(x) for x in args], stdout=stream, check=True)
        return None
    return subprocess.check_output([str(x) for x in args], text=True)


def verify_sources(root, acquire=False, large_downloads=False):
    entries = json.loads(MANIFEST.read_text())['datasets']
    paths = {}; records = []
    for entry in entries:
        path = (root / entry['local_relative_path']).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Manifest path escapes data root')
        if not path.exists():
            if not acquire:
                raise ValueError('Missing verified source; use --acquire: ' + str(path))
            if entry['bytes'] > 100_000_000 and not large_downloads:
                raise ValueError('Large acquisition requires --allow-large-downloads: ' + entry['id'])
            if entry.get('derived_from'):
                parent = paths[entry['derived_from']]
                path.parent.mkdir(parents=True, exist_ok=True)
                if entry['id'] == 'genome':
                    ACQ.decompress_fasta(parent, path)
                elif entry['id'] == 'genome-fai':
                    ACQ.ensure_fasta(parent)
            else:
                def validate(candidate, record=entry):
                    ACQ.validate(candidate, record['kind'])
                    if candidate.stat().st_size != record['bytes'] or sha256(candidate) != record['sha256']:
                        raise ValueError('Downloaded source checksum/size differs from pinned manifest: ' + record['id'])
                ACQ.acquire(entry['url'], path, validator=validate)
        if path.stat().st_size != entry['bytes'] or sha256(path) != entry['sha256']:
            raise ValueError('Cached source checksum/size mismatch; retained unchanged: ' + entry['id'])
        paths[entry['id']] = path
        record = dict(entry, verified_sha256=entry['sha256'])
        if entry['kind'] == 'vcf':
            record.update(vcf_header_facts(path))
        records.append(record)
    return paths, records


def read_vcf(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    metadata = []; header = None; records = []; definitions = {'INFO':{},'FORMAT':{}}
    with opener(path, 'rt') as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.rstrip('\r\n')
            if line.startswith('##'):
                metadata.append(line)
                match = re.match(r'##(INFO|FORMAT)=<(.*)>$', line)
                if match:
                    id_match = re.search(r'(?:^|,)ID=([^,>]+)', match[2])
                    number = re.search(r'(?:^|,)Number=([^,>]+)', match[2])
                    if id_match and number: definitions[match[1]][id_match[1]] = number[1]
            elif line.startswith('#CHROM'):
                header = line.split('\t')
            elif line and not line.startswith('#'):
                records.append((line_number, line.split('\t')))
    if header is None: raise ValueError('Independent inspection: missing #CHROM header')
    return metadata, header, records, definitions


def write_vcf(path, metadata, header, records):
    with Path(path).open('w') as stream:
        for line in metadata: stream.write(line + '\n')
        stream.write('\t'.join(header) + '\n')
        for fields in records: stream.write('\t'.join(fields) + '\n')


def prepare(paths, root):
    derived = root / 'public_datasets/vcf-compatibility-v1/derived'
    derived.mkdir(parents=True, exist_ok=True)
    full_ref = paths['genome']
    # Explicit reference derivation, never an implicit parser contig repair.
    numeric_ref = derived / 'GRCh38-chr22-as-22.fa'
    temp = derived / 'chr22-extract.fa'
    command(['samtools','faidx',full_ref,'chr22'], temp)
    staged_ref = derived / 'chr22-numeric-staged.fa'
    with temp.open('rb') as src, staged_ref.open('wb') as dst:
        src.readline(); dst.write(b'>22\n'); shutil.copyfileobj(src,dst)
    os.replace(staged_ref, numeric_ref)
    temp.unlink()
    ACQ.ensure_fasta(numeric_ref)
    ACQ.ensure_fasta(paths['strelka-demo'])
    reference_hashes = {str(full_ref):sha256(full_ref), str(numeric_ref):sha256(numeric_ref),
                        str(paths['strelka-demo']):sha256(paths['strelka-demo'])}
    cases = []
    def add(id, path, reference=full_ref, sample='HCC1395', role='somatic_vcf', assembly='GRCh38.p14',
            expected='SUPPORTED_PASS', error=None, parent=None, operation=None, equivalent=None, large=False):
        cases.append(dict(id=id,path=str(path),reference=str(reference),sample_id=sample,patient_id='P1',role=role,
                          assembly=assembly,expected=expected,error=error,parent=parent,preprocessing=operation,
                          equivalent_to=equivalent,large=large,sha256=sha256(path),bytes=Path(path).stat().st_size,
                          reference_sha256=reference_hashes[str(reference)],reference_fai_sha256=sha256(str(reference)+'.fai')))
    add('seqc2-quick',paths['seqc2-quick'],parent='seqc2-quick',operation='none; original site-only VCF')
    meta,header,records,_ = read_vcf(paths['seqc2-full'])
    subset = [r[1] for r in records[:256]]
    plain = derived / 'seqc2-first256.vcf'; write_vcf(plain,meta,header,subset)
    compressed = derived / 'seqc2-first256.vcf.gz'
    with plain.open('rb') as src, compressed.open('wb') as dst:
        with gzip.GzipFile(filename='',mode='wb',fileobj=dst,mtime=0) as zipped: shutil.copyfileobj(src,zipped)
    add('seqc2-subset',plain,parent='seqc2-full',operation='first 256 records, original header, no allele/header edits')
    add('seqc2-gzip',compressed,parent='seqc2-full',operation='gzip of seqc2-subset',equivalent='seqc2-subset')
    reverse = derived / 'seqc2-reversed.vcf'; write_vcf(reverse,meta,header,list(reversed(subset)))
    add('seqc2-reversed',reverse,parent='seqc2-full',operation='reverse the same 256 rows',equivalent='seqc2-subset')
    duplicate = derived / 'seqc2-duplicate.vcf'; write_vcf(duplicate,meta,header,subset+[subset[0]])
    add('seqc2-duplicate',duplicate,parent='seqc2-full',operation='same 256 rows plus one explicit duplicate',equivalent='seqc2-subset')
    for name in ['strelka-snvs','strelka-indels']:
        original=paths[name]
        add(name+'-multisample',original,reference=paths['strelka-demo'],sample='NORMAL',assembly='STRELKA_DEMO',
            expected='EXPECTED_REJECTION',error='Multi-sample VCF',parent=name,operation='none; NORMAL,TUMOR header retained')
        output=derived/(name+'-TUMOR.vcf')
        command(['bcftools','view','--no-version','--no-update','-s','TUMOR','-Ov','-o',output,original])
        add(name+'-TUMOR',output,reference=paths['strelka-demo'],sample='TUMOR',assembly='STRELKA_DEMO',parent=name,
            operation='bcftools view --no-version --no-update -s TUMOR -Ov; no header/allele repair')
    add('mutect2-mito-cardinality',paths['mutect2-mito'],sample='NA12878',expected='EXPECTED_REJECTION',
        error='FORMAT MMQ Number=A cardinality',parent='mutect2-mito',operation='none; upstream MMQ declares Number=A but has ref,alt values')
    add('mutect2-na12878-multisample',paths['mutect2-na12878'],sample='SM-612V3',expected='EXPECTED_REJECTION',
        error='Multi-sample VCF',parent='mutect2-na12878',operation='none; both original samples retained')
    # Indexed public panel extraction is explicit, independent preprocessing.
    igsr = derived/'1000g-NA12878-region.vcf'
    command(['bcftools','view','--no-version','--no-update','-s','NA12878','-r','22:20000000-20050000','-Ov','-o',igsr,paths['1000g-panel']])
    add('1000g-selected',igsr,reference=numeric_ref,sample='NA12878',role='germline_vcf',parent='1000g-panel',
        operation='bcftools single-sample/region extraction; reference chr22 sequence explicitly selected and named 22')
    add('1000g-contig-mismatch',igsr,sample='NA12878',role='germline_vcf',expected='EXPECTED_REJECTION',error='Contig absent',
        parent='1000g-panel',operation='same VCF against unmodified chr-prefixed primary assembly; no alias repair')
    add('1000g-sample-mismatch',igsr,reference=numeric_ref,sample='WRONG',role='germline_vcf',expected='EXPECTED_REJECTION',
        error='sample does not match',parent='1000g-panel',operation='same VCF with deliberately incorrect authoritative sample binding')
    panel=derived/'1000g-panel-first2.vcf'
    # Retain all 2,548 sample columns, but only two complete source records.
    with gzip.open(paths['1000g-panel'],'rt') as src, panel.open('w') as dst:
        count=0
        for line in src:
            dst.write(line)
            if not line.startswith('#'): count+=1
            if count==2:break
    add('1000g-multisample',panel,reference=numeric_ref,sample='HG00096',role='germline_vcf',expected='EXPECTED_REJECTION',
        error='Multi-sample VCF',parent='1000g-panel',operation='first two records; original multi-sample header retained')
    bad_ref=derived/'seqc2-wrong-ref.vcf'; row=list(subset[0]); row[3]=next(b for b in 'ACGT' if b!=row[3] and b!=row[4])
    write_vcf(bad_ref,meta,header,[row])
    add('seqc2-ref-mismatch',bad_ref,expected='EXPECTED_REJECTION',error='REF mismatch',parent='seqc2-full',operation='negative control: explicitly mutate one REF base')
    bad_build=derived/'seqc2-wrong-build.vcf'
    write_vcf(bad_build,[m if not m.startswith('##reference=') else '##reference=GRCh37' for m in meta],header,subset[:1])
    add('seqc2-build-mismatch',bad_build,expected='EXPECTED_REJECTION',error='contradicts selected reference',parent='seqc2-full',operation='negative control: explicitly replace reference declaration with GRCh37')
    corrupt=derived/'seqc2-corrupt.vcf.gz'; contents=bytearray(compressed.read_bytes()); contents[-8:]=b'\x00'*8; corrupt.write_bytes(contents)
    add('seqc2-corrupt-gzip',corrupt,expected='EXPECTED_REJECTION',error='error reading|compressed|compression|CRC',parent='seqc2-full',operation='negative control: corrupt gzip CRC/ISIZE trailer')
    truncated=derived/'seqc2-truncated.vcf'; write_vcf(truncated,meta,header,subset[:-1]+[subset[-1][:5]])
    add('seqc2-truncated-record',truncated,expected='EXPECTED_REJECTION',error='Malformed VCF record',parent='seqc2-full',operation='negative control: truncate final row after ALT')
    missing_header=derived/'seqc2-missing-header.vcf'
    missing_header.write_text('\t'.join(subset[0])+'\n')
    add('seqc2-missing-header',missing_header,expected='EXPECTED_REJECTION',error='VCF data before column header',
        parent='seqc2-full',operation='negative control: one original row with all headers deliberately removed')
    add('seqc2-full-large',paths['seqc2-full'],parent='seqc2-full',operation='none; all 41,072 original records',large=True)
    return cases


def source_facts(case):
    meta,header,records,defs = read_vcf(case['path'])
    expected=[]; ref_filtered=0; raw_alts=0; malformed=0; kinds=collections.Counter(); cardinality=[]
    for line,fields in records:
        if len(fields)!=len(header):malformed+=1;continue
        chrom,pos,id,ref,alt=fields[:5];alts=alt.split(','); raw_alts+=len(alts)
        fmt=dict(zip(fields[8].split(':'),fields[9].split(':'))) if len(fields)==10 else {}
        called=list(range(1,len(alts)+1))
        if case['role']=='germline_vcf' and 'GT' in fmt:
            tokens=re.split(r'[/|]',fmt['GT'])
            if any(not t.isdigit() for t in tokens): raise ValueError('Unresolved germline oracle GT')
            called=sorted({int(t) for t in tokens if int(t)>0})
            if not called:ref_filtered+=1
        info={item.split('=',1)[0]:item.split('=',1)[1] if '=' in item else None for item in fields[7].split(';') if item!='.'}
        for category, values in [('INFO',info),('FORMAT',fmt)]:
            for key,value in values.items():
                number=defs[category].get(key)
                if number in ('A','R') and value not in (None,'.'):
                    wanted=len(alts)+(number=='R')
                    if len(value.split(','))!=wanted:cardinality.append(dict(line=line,field=category+'/'+key,number=number,value=value,expected=wanted))
        for index in called:
            value=alts[index-1]
            kinds['SNV' if len(ref)==len(value)==1 else 'MNV' if len(ref)==len(value) else 'indel']+=1
            expected.append(dict(line=line,index=index,fields=fields,alt=value,info=info,format=fmt))
    summary=dict(input_records=len(records),input_alt_alleles=raw_alts,selected_alt_alleles=len(expected),
                 hom_ref_filtered_records=ref_filtered,malformed_records=malformed,variant_types=dict(kinds),
                 sample_count=max(0,len(header)-9),declared_samples=header[9:13],cardinality_errors=cardinality[:10])
    return summary,expected,defs


def independent_oracle(case, directory, expected):
    oracle=directory/'oracle-input.vcf'; fai=Path(case['reference']+'.fai')
    contigs=[line.split('\t')[:2] for line in fai.read_text().splitlines()]
    meta=['##fileformat=VCFv4.2']+['##contig=<ID='+c+',length='+n+'>' for c,n in contigs]
    rows=[]
    for i,e in enumerate(expected):
        f=e['fields'];rows.append([f[0],f[1],'oracle_'+str(i),f[3],e['alt'],'.','PASS','.'])
    write_vcf(oracle,meta,['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'],rows)
    norm=directory/'oracle-normalized.vcf'
    probe=subprocess.run(['bcftools','norm','--no-version','-f',case['reference'],'-Ov','-o',str(norm),str(oracle)],capture_output=True,text=True)
    (directory/'bcftools-norm.log').write_text(probe.stderr)
    if probe.returncode: raise ValueError('Independent normalization failed: '+probe.stderr[:400])
    normalized={}
    for _,fields in read_vcf(norm)[2]:
        index=int(fields[2].split('_')[1]); e=expected[index]
        normalized[(e['line'],e['index'])]=(fields[0],int(fields[1]),fields[3],fields[4])
    return normalized


def projections(values,definitions,index,alt_count):
    out={}
    for key,value in values.items():
        number=definitions.get(key)
        if number not in ('A','R') or value is None:continue
        wanted=alt_count+(number=='R')
        if value=='.' and wanted>1:continue
        items=value.split(',')
        if len(items)!=wanted:raise ValueError('Invalid oracle cardinality '+key)
        out[key]=items[index-1] if number=='A' else {'ref':items[0],'alt':items[index]}
    return out


def check_semantics(case,result,expected,defs,normalized):
    errors=[]; actual={}; ids=[]
    if result['raw_event_count']!=len(expected):errors.append('selected ALT expansion count differs')
    expected_unique=set(normalized.values())
    if result['unique_event_count']!=len(expected_unique):errors.append('unique normalized event count differs')
    for row in result['rows']:
        ids.append(row['event_id'])
        locus=(row['contig'],row['pos'],row['ref'],row['alt'])
        if row['sample_id']!=case['sample_id'] or row['patient_id']!=case['patient_id']:errors.append('identity mismatch')
        for source in row['sources']:
            key=(source['line'],source['alt_index'])
            if key in actual:errors.append('duplicate source occurrence')
            actual[key]=(locus,source)
    if set(actual)!=set(normalized):errors.append('missing/extra source ALT occurrences')
    for e in expected:
        key=(e['line'],e['index'])
        if key not in actual:continue
        locus,source=actual[key];fields=e['fields'];count=len(fields[4].split(','))
        if locus!=normalized[key]:errors.append('normalized REF/ALT/position mismatch at '+str(key))
        if source['sample_id']!=case['sample_id']:errors.append('source sample mismatch')
        if source['raw_info']!=fields[7]:errors.append('raw INFO changed')
        info=projections(e['info'],defs['INFO'],e['index'],count)
        fmt=projections(e['format'],defs['FORMAT'],e['index'],count)
        if (source['projected_info'] or {})!=info:errors.append('ALT-specific INFO mismatch at '+str(key))
        if (source['projected_format'] or {})!=fmt:errors.append('ALT-specific FORMAT mismatch at '+str(key))
        if len(fields)==10:
            if source['raw_format']!={'keys':fields[8],'sample':fields[9]}:errors.append('raw FORMAT changed')
            gt=e['format'].get('GT');phase=source['phase']
            if gt and re.fullmatch(r'\d+[|/]\d+',gt):
                alleles=[int(x) for x in re.split(r'[|/]',gt)]
                side='BOTH' if alleles.count(e['index'])==2 else 'LEFT' if alleles[0]==e['index'] else 'RIGHT' if alleles[1]==e['index'] else 'ABSENT'
                status='INTERPRETABLE' if '|' in gt else 'UNPHASED'
                if not any(a.get('status')==status and a.get('genotype_state')==side for a in phase):errors.append('GT/phase interpretation mismatch')
    if ids!=sorted(ids):errors.append('event ordering is not deterministic ID order')
    if len(ids)!=len(set(ids)):errors.append('duplicate occurrence IDs after merge')
    return sorted(set(errors))


def evaluate(case,run_root,timeout):
    directory=run_root/case['id'];directory.mkdir(parents=True,exist_ok=True)
    record={k:v for k,v in case.items() if k not in ('path','reference')}
    record.update(path=case['path'],reference=case['reference'],scope='VCF component import, merge and per-source genotype/phase assertions; no peptide stages')
    expected=[];defs={};normalized=None
    try:
        facts,expected,defs=source_facts(case); record.update(facts)
        if case['expected']=='SUPPORTED_PASS':normalized=independent_oracle(case,directory,expected)
    except (ValueError,EOFError,OSError,IndexError) as error:
        record['oracle_diagnostic']=str(error)
        if case['expected']=='SUPPORTED_PASS':record['outcome']='UNRESOLVED';return record
    request=directory/'request.json';request.write_text(json.dumps(case,indent=2)+'\n')
    output=directory/'import.json';time_flag='-l' if platform.system()=='Darwin' else '-v'
    cmd=['/usr/bin/time',time_flag,'Rscript',str(REPO/'scripts/development/vcf_compatibility_import.R'),str(REPO),str(request),str(output)]
    started=time.monotonic()
    try:
        process=subprocess.run(cmd,cwd=REPO,capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        record.update(outcome='UNRESOLVED',resource_limit='R import exceeded configured timeout',timeout_seconds=timeout)
        return record
    record.update(wall_seconds=round(time.monotonic()-started,3),command=cmd,process_exit_code=process.returncode)
    (directory/'process.log').write_text(process.stdout+'\n'+process.stderr)
    rss=re.search(r'(\d+)\s+maximum resident set size',process.stderr) if platform.system()=='Darwin' else re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',process.stderr)
    record['peak_rss_bytes']=int(rss[1])*(1 if platform.system()=='Darwin' else 1024) if rss else None
    if process.returncode or not output.exists():record.update(outcome='UNEXPECTED_ERROR',error='R process failed without component result');return record
    result=json.loads(output.read_text());record['warnings']=result.get('warnings',[])
    if not result['success']:
        record['error']=result['error'];record['emitted_events']=0
        record['outcome']='EXPECTED_REJECTION' if case['expected']=='EXPECTED_REJECTION' and re.search(case['error'],result['error']) else 'UNEXPECTED_ERROR'
        return record
    record.update(emitted_events=result['unique_event_count'],accepted_alt_occurrences=result['raw_event_count'])
    if case['expected']=='EXPECTED_REJECTION':record.update(outcome='INCORRECT_INTERPRETATION',errors=['Unsupported/invalid input accepted']);return record
    errors=check_semantics(case,result,expected,defs,normalized)
    record.update(outcome='INCORRECT_INTERPRETATION' if errors else 'SUPPORTED_PASS',errors=errors)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',required=True,type=Path)
    parser.add_argument('--run-root',required=True,type=Path)
    parser.add_argument('--acquire',action='store_true',help='Reuse existing acquisition validation/atomic downloads for missing pinned inputs')
    parser.add_argument('--allow-large-downloads',action='store_true')
    parser.add_argument('--large',action='store_true',help='Include all 41,072 SEQC2 records in a separate timed import')
    parser.add_argument('--case',action='append',help='Only these prepared case identifiers (repeatable)')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--timeout',type=int,default=600,help='Seconds per isolated R import; timeouts are resource-limited UNRESOLVED')
    args=parser.parse_args();root=args.data_root.expanduser().resolve();run=args.run_root.expanduser().resolve()
    run.mkdir(parents=True,exist_ok=True);(run/'tmp').mkdir(exist_ok=True)
    os.environ['TMPDIR']=str(run/'tmp')
    ACQ.require('bcftools','samtools','Rscript','gzip')
    paths,verified=verify_sources(root,args.acquire,args.allow_large_downloads)
    ACQ.atomic_json(run/'verified-sources.json',verified)
    cases=prepare(paths,root)
    chosen=[c for c in cases if (args.large or not c['large']) and (not args.case or c['id'] in args.case)]
    if args.case and set(args.case)-{c['id'] for c in chosen}:raise ValueError('Unknown/unselected --case; large cases also require --large')
    ACQ.atomic_json(run/'prepared-cases.json',cases)
    ACQ.atomic_json(root/'acquisition/vcf-compatibility.json',dict(
        schema_version='vcf-compatibility-1',manifest_sha256=sha256(MANIFEST),
        acquisition_script='scripts/development/validate_vcf_compatibility.py',
        source_artifacts=verified,derived_cases=cases))
    versions={tool:command([tool,'--version']).splitlines()[0] for tool in ['bcftools','samtools']}
    ACQ.atomic_json(run/'environment.json',dict(versions=versions,python=sys.version,platform=platform.platform(),manifest_sha256=sha256(MANIFEST)))
    if args.prepare_only:return
    results=[]
    for case in chosen:
        print('[case]',case['id'],flush=True)
        record=evaluate(case,run,args.timeout)
        if record['outcome']=='SUPPORTED_PASS' and case.get('equivalent_to'):
            prior=next((r for r in results if r['id']==case['equivalent_to']),None)
            if prior and prior['outcome']=='SUPPORTED_PASS':
                def signature(id):
                    rows=json.loads((run/id/'import.json').read_text())['rows']
                    return [(r['event_id'],r['event_key'],r['contig'],r['pos'],r['ref'],r['alt']) for r in rows]
                if signature(case['id'])!=signature(prior['id']):record.update(outcome='INCORRECT_INTERPRETATION',errors=['Equivalent input changed ordered semantic/ID output'])
        results.append(record)
        ACQ.atomic_json(run/'compatibility-results.json',results)
        print(record['outcome'],case['id'],'records=',record.get('input_records'),'events=',record.get('emitted_events'),flush=True)
    counts={category:sum(r['outcome']==category for r in results) for category in CATEGORIES}
    ACQ.atomic_json(run/'compatibility-summary.json',dict(case_denominator=len(results),counts=counts,notes='Cases are file-level operations, not independent biological samples. Rejected files emit no partial events. Do not add repeated/equivalent source record denominators as independent evidence.'))
    print(json.dumps(counts,indent=2))
    if counts['UNEXPECTED_ERROR'] or counts['INCORRECT_INTERPRETATION'] or counts['UNRESOLVED']:return 1
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,OSError,subprocess.CalledProcessError) as error:sys.exit('[error] '+str(error))
