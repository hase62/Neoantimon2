#!/usr/bin/env python3
"""Public input acquisition; standalone script, independent of the R installation.

Compatibility uses declared assembly/release/accession/scope, not mandatory
historic byte hashes. Completed inputs are validated before reuse or promotion.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

IMGTHLA_COMMIT = '5b915f27f7f620361cf83cb626eeac8a03c0247c'
IMGTHLA_RELEASE = '3.65.0'
RESOURCES = ('genome', 'gencode', 'mane', 'imgthla', '1000g', 'hcc1395',
             'processed-rna', 'longread', 'full-rna', 'full-truth')
PROFILES = {'minimum': ('genome', 'gencode', 'mane'), 'full': RESOURCES}


def run(*args, capture=False):
    result = subprocess.run([str(a) for a in args], check=True,
                            stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
                            text=capture)
    return result.stdout if capture else None


def require(*commands):
    missing = [c for c in commands if not shutil.which(c)]
    if missing:
        raise RuntimeError('Required command(s): ' + ', '.join(missing))


def validate(path, kind=None):
    """Readability and format checks also reject nonempty truncated downloads."""
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        raise ValueError('Missing/empty artifact: ' + str(p))
    name = p.name.lower()
    kind = kind or ('bai' if name.endswith('.bai') else 'tbi' if name.endswith('.tbi')
                    else 'fai' if name.endswith('.fai') else 'vcf' if '.vcf' in name
                    else 'bam' if name.endswith('.bam') else 'fasta' if '.fa' in name
                    else 'gtf' if '.gtf' in name else 'bed' if name.endswith('.bed') else 'text')
    # Avoid reading a large artifact into RAM. Identify gzip with two bytes.
    with p.open('rb') as stream:
        magic = stream.read(4)
    compressed = magic[:2] == b'\x1f\x8b'
    if name.endswith('.gz') or compressed or kind == 'tbi':
        require('gzip')
        run('gzip', '-t', p)
    if kind == 'vcf':
        require('bcftools')
        run('bcftools', 'view', '-Ou', p)
        header = run('bcftools', 'view', '-h', p, capture=True)
        if '#CHROM\tPOS\tID\tREF\tALT' not in header:
            raise ValueError('Invalid VCF header')
    elif kind == 'bam':
        require('samtools')
        run('samtools', 'quickcheck', '-v', p)
        run('samtools', 'view', '-H', p)
    elif kind == 'bai':
        if magic != b'BAI\x01' or p.stat().st_size < 8:
            raise ValueError('Invalid BAM index')
        with p.open('rb') as stream:
            stream.read(4)
            if struct.unpack('<i', stream.read(4))[0] < 1:
                raise ValueError('Empty BAM index')
    elif kind == 'tbi':
        with gzip.open(p, 'rb') as stream:
            if stream.read(4) != b'TBI\x01' or len(stream.read(4)) != 4:
                raise ValueError('Invalid tabix index')
    elif kind == 'fastq':
        opener = gzip.open if compressed else open
        records = 0
        with opener(p, 'rt') as stream:
            while True:
                name = stream.readline()
                if not name: break
                sequence, plus, quality = stream.readline(), stream.readline(), stream.readline()
                if not name.startswith('@') or not plus.startswith('+') or not quality or len(sequence.rstrip()) != len(quality.rstrip()):
                    raise ValueError('Invalid/incomplete FASTQ')
                records += 1
        if not records: raise ValueError('Empty FASTQ')
    else:
        opener = gzip.open if compressed else open
        with opener(p, 'rt') as stream:
            lines, headers = [], []
            for _ in range(200):
                line = stream.readline()
                if not line:
                    break
                if line.startswith('#'): headers.append(line.rstrip('\n'))
                if line.strip() and not line.startswith('#'):
                    lines.append(line.rstrip('\n'))
        if not lines:
            raise ValueError('No readable content: ' + str(p))
        if kind == 'fasta' and (not lines[0].startswith('>') or
                                not any(line and not line.startswith('>') for line in lines[1:])):
            raise ValueError('Invalid FASTA')
        if kind in ('gtf', 'gencode-gtf') and any(len(line.split('\t')) != 9 for line in lines):
            raise ValueError('Invalid GTF')
        if kind == 'gencode-gtf' and not any('GRCh38' in line and re.search(r'version 50\b', line) for line in headers):
            raise ValueError('Expected GENCODE 50 / GRCh38 annotation header absent')
        if kind == 'mane-summary' and not any('NCBI_GeneID' in line and 'MANE_status' in line and 'GRCh38_chr' in line for line in headers + lines):
            raise ValueError('Expected MANE/GRCh38 summary columns absent')
        if kind == 'sec' and any(len(line.split('\t')) != 2 or not line.split('\t')[0].startswith('ENST') or not line.split('\t')[1].isdigit() for line in lines):
            raise ValueError('Invalid GENCODE Sec metadata')
        if kind == 'bed'  and any(len(line.split('\t')) < 3 for line in lines):
            raise ValueError('Invalid BED')
        if kind == 'fai':
            for line in lines:
                fields = line.split('\t')
                if len(fields) < 5 or any(int(x) < 0 for x in fields[1:5]) or int(fields[1]) == 0:
                    raise ValueError('Invalid FASTA index')
    return p


def fasta_scope(path, expected, exact=False):
    validate(path, 'fasta')
    with Path(path).open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    lengths = {}; current = None
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('>'):
                current = line[1:].split()[0]
                if current in lengths: raise ValueError('Duplicate FASTA contig')
                lengths[current] = 0
            elif line.strip():
                if current is None: raise ValueError('FASTA sequence before header')
                lengths[current] += len(line.strip())
    if any(lengths.get(name) != length for name, length in expected.items()) or (exact and set(lengths) != set(expected)):
        raise ValueError('FASTA does not match selected assembly/contig scope')
    return Path(path)


def acquire(url, output, kind=None, validator=None):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    check = validator or (lambda p: validate(p, kind))
    if output.exists():
        try:
            check(output)
            print('[validated reuse]', output, flush=True)
            return output
        except (ValueError, OSError, subprocess.CalledProcessError) as error:
            print('[invalid existing artifact; acquiring replacement]', output, str(error), flush=True)
    require('curl')
    # Each attempt owns fresh staging state. Never resume/reuse an unverified
    # previous destination or partial file. The old destination is not removed.
    with tempfile.TemporaryDirectory(prefix='.neo-download-', dir=output.parent) as directory:
        stage = Path(directory) / output.name
        run('curl', '--location', '--fail', '--retry', '3', '--retry-delay', '2',
            '--output', stage, url)
        check(stage)
        os.replace(stage, output)
    print('[completed]', output, flush=True)
    return output


def acquire_pair(url, output, suffix, pair_validator):
    """Validate a data/index pair in staging before promoting either member."""
    output = Path(output)
    index = Path(str(output) + suffix)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and index.exists():
        try:
            pair_validator(output)
            print('[validated pair reuse]', output, flush=True)
            return output
        except (ValueError, OSError, subprocess.CalledProcessError):
            pass
    with tempfile.TemporaryDirectory(prefix='.neo-pair-', dir=output.parent) as directory:
        stage = Path(directory) / output.name
        stage_index = Path(str(stage) + suffix)
        reused = set()
        for prior, candidate, endpoint in ((output, stage, url), (index, stage_index, url + suffix)):
            if prior.exists():
                try:
                    validate(prior)
                    os.link(prior, candidate)
                    reused.add(candidate)
                except (ValueError, OSError, subprocess.CalledProcessError):
                    acquire(endpoint, candidate)
            else:
                acquire(endpoint, candidate)
        try:
            pair_validator(stage)
        except (ValueError, OSError, subprocess.CalledProcessError):
            if stage_index not in reused:
                raise
            # This is only our temporary hard link, not the existing input.
            stage_index.unlink()
            reused.remove(stage_index)
            acquire(url + suffix, stage_index)
            pair_validator(stage)
        if stage_index not in reused: os.replace(stage_index, index)
        if stage not in reused: os.replace(stage, output)
    return output


def vcf_pair(path, region=None, sample=None):
    p = validate(path, 'vcf')
    validate(str(p) + '.tbi', 'tbi')
    count = int(run('bcftools', 'index', '-n', p, capture=True).strip())
    if count < 1:
        raise ValueError('Expected nonzero VCF records')
    if sample and sample not in run('bcftools', 'query', '-l', p, capture=True).splitlines():
        raise ValueError('Expected sample absent: ' + sample)
    if region:
        text = run('bcftools', 'view', '-H', '-r', region, p, capture=True)
        if not text.strip():
            raise ValueError('Expected variants absent from region: ' + region)
    else:
        contigs = [line.split('\t')[0] for line in run('bcftools', 'index', '-s', p, capture=True).splitlines()]
        if not contigs:
            raise ValueError('No indexed contigs')
        run('bcftools', 'view', '-H', '-r', contigs[0] + ':1-1', p)
    return count


def bam_pair(path):
    p = validate(path, 'bam')
    validate(str(p) + '.bai', 'bai')
    header = run('samtools', 'view', '-H', p, capture=True)
    contigs = [part[3:] for line in header.splitlines() if line.startswith('@SQ')
               for part in line.split('\t') if part.startswith('SN:')]
    if not contigs:
        raise ValueError('BAM reference header absent')
    run('samtools', 'view', '-c', p, contigs[0] + ':1-1')


def index_bam(path):
    """Supported samtools index syntax; promote only a readable new index."""
    p = validate(path, 'bam')
    with tempfile.TemporaryDirectory(prefix='.neo-index-', dir=p.parent) as directory:
        stage = Path(directory) / (p.name + '.bai')
        run('samtools', 'index', p, stage)
        validate(stage, 'bai')
        os.replace(stage, str(p) + '.bai')
    bam_pair(p)


def ensure_fasta(path):
    p = validate(path, 'fasta')
    require('samtools')
    with tempfile.TemporaryDirectory(prefix='.neo-fasta-', dir=p.parent) as directory:
        stage = Path(directory) / (p.name + '.fai')
        run('samtools', 'faidx', '-o', stage, p)
        validate(stage, 'fai')
        target = Path(str(p) + '.fai')
        if not target.exists() or target.read_bytes() != stage.read_bytes():
            os.replace(stage, target)
    return p


def decompress_fasta(source, output):
    output = Path(output)
    validate(source, 'fasta')
    matches = output.is_file()
    if matches:
        with gzip.open(source, 'rb') as src, output.open('rb') as dst:
            while True:
                a, b = src.read(1024 * 1024), dst.read(1024 * 1024)
                if a != b:
                    matches = False
                    break
                if not a: break
    if not matches:
        with tempfile.TemporaryDirectory(prefix='.neo-fasta-', dir=output.parent) as directory:
            stage = Path(directory) / output.name
            with gzip.open(source, 'rb') as src, stage.open('wb') as dst:
                shutil.copyfileobj(src, dst)
            validate(stage, 'fasta')
            os.replace(stage, output)
    return ensure_fasta(output)


def extract_1000g(panel, destination):
    require('bcftools')
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    vcf_pair(panel, sample='NA12878')
    contigs = [line.split('\t')[0] for line in run('bcftools', 'index', '-s', panel, capture=True).splitlines()]
    selected = [c for c in contigs if c in ('22', 'chr22')]
    if len(selected) != 1:
        raise ValueError('Expected exactly one chr22/22 source contig')
    region = selected[0] + ':20000000-25000000'
    for name, arguments in [('NA12878.chr22.GRCh38.phased.vcf.gz', []),
                            ('NA12878.chr22.20-25Mb.GRCh38.phased.vcf.gz', ['-r', region])]:
        out = destination / name
        if out.exists():
            try:
                vcf_pair(out, region=region, sample='NA12878')
                continue
            except (ValueError, OSError, subprocess.CalledProcessError):
                pass
        with tempfile.TemporaryDirectory(prefix='.neo-1000g-', dir=destination) as directory:
            stage = Path(directory) / name
            run('bcftools', 'view', '-s', 'NA12878', *arguments, '-Oz', '-o', stage, panel)
            run('bcftools', 'index', '-t', stage)
            vcf_pair(stage, region=region, sample='NA12878')
            # Keep the source contig name. Acquisition reconciles the query;
            # Neoantimon2's exact reference-contig semantics remain unchanged.
            os.replace(str(stage) + '.tbi', str(out) + '.tbi')
            os.replace(stage, out)
    return {'source_contig': selected[0], 'subset_region': region,
            'subset_records': vcf_pair(destination / 'NA12878.chr22.20-25Mb.GRCh38.phased.vcf.gz', region=region)}


def acquire_imgthla(destination, commit=IMGTHLA_COMMIT, release=IMGTHLA_RELEASE):
    require('git')
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    def check(p):
        if run('git', '-C', p, 'rev-parse', 'HEAD', capture=True).strip() != commit:
            raise ValueError('Existing IMGTHLA checkout differs from selected commit; choose another destination')
        if run('git', '-C', p, 'status', '--porcelain', capture=True).strip():
            raise ValueError('IMGTHLA checkout contains local changes')
        if release not in (p / 'release_version.txt').read_text():
            raise ValueError('IMGTHLA release mismatch')
    if destination.exists():
        check(destination)
    else:
        with tempfile.TemporaryDirectory(prefix='.neo-imgthla-', dir=destination.parent) as directory:
            stage = Path(directory) / 'IMGTHLA'
            run('git', 'clone', '--no-checkout', 'https://github.com/ANHIG/IMGTHLA.git', stage)
            run('git', '-C', stage, 'checkout', '--detach', commit)
            check(stage)
            os.replace(stage, destination)
    return destination


def acquire_sra(accession, directory):
    require('prefetch', 'vdb-validate')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    out = directory / accession / (accession + '.sra')
    if out.exists():
        try:
            run('vdb-validate', out)
            return out
        except subprocess.CalledProcessError:
            print('[invalid existing SRA; staging replacement]', out, flush=True)
    with tempfile.TemporaryDirectory(prefix='.neo-sra-', dir=directory) as stage:
        run('prefetch', '--max-size', 'u', '--output-directory', stage, accession)
        candidate = Path(stage) / accession / (accession + '.sra')
        if not candidate.is_file():
            raise ValueError('Expected accession file absent')
        run('vdb-validate', candidate)
        out.parent.mkdir(parents=True, exist_ok=True)
        os.replace(candidate, out)
    return out


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False, prefix='.neo-metadata-') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stage = Path(stream.name)
    os.replace(stage, path)


def acquire_resource(resource, root, observed_sha256=False):
    root = Path(root)
    files, sources = [], []
    def get(url, relative, kind=None):
        p = acquire(url, root / relative, kind, validator=(lambda p: fasta_scope(p, {'chr1':248956422,'chr22':50818468,'chrM':16569})) if resource == 'genome' else (lambda p: fasta_scope(p, {'chr1':248956422}, exact=True)) if relative.endswith('GRCh38_no_alts_chr1.fasta') else None)
        files.append(p)
        sources.append(url)
        return p
    def get_pair(url, relative, suffix, checker):
        p = acquire_pair(url, root / relative, suffix, checker)
        files.extend([p, Path(str(p) + suffix)])
        sources.extend([url, url + suffix])
        return p
    scope = {}
    release = None
    license_url = None
    if resource in ('genome', 'gencode'):
        base = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/'
        release = 'GENCODE 50 / GRCh38.p14'
        license_url = 'https://www.gencodegenes.org/pages/data_access.html'
        if resource == 'genome':
            p = get(base + 'GRCh38.primary_assembly.genome.fa.gz', 'reference/genome/GRCh38.primary_assembly.genome.fa.gz', 'fasta')
            files.extend([decompress_fasta(p, root / 'reference/genome/GRCh38.primary_assembly.genome.fa'), root / 'reference/genome/GRCh38.primary_assembly.genome.fa.fai'])
            scope = {'assembly': 'GRCh38.p14', 'scope': 'primary_assembly'}
        else:
            for name in ('primary_assembly.annotation.gtf', 'transcripts.fa', 'pc_transcripts.fa', 'pc_translations.fa', 'metadata.Selenocysteine'):
                filename = 'gencode.v50.' + name + '.gz'
                get(base + filename, 'reference/gencode/' + filename, 'gencode-gtf' if name.endswith('.gtf') else 'sec' if name=='metadata.Selenocysteine' else None)
            scope = {'annotation_release': '50', 'annotation_scope': 'primary_assembly', 'catalogs': 'all'}
    elif resource == 'mane':
        license_url = 'https://www.ncbi.nlm.nih.gov/home/about/policies/'
        release = 'MANE 1.5 / GRCh38'
        base = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/'
        for filename in ('MANE.GRCh38.v1.5.summary.txt.gz', 'MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz'):
            get(base + filename, 'reference/mane/' + filename, 'mane-summary' if filename.endswith('summary.txt.gz') else None)
    elif resource == 'imgthla':
        release = IMGTHLA_RELEASE
        p = acquire_imgthla(root / 'reference/hla/IMGTHLA')
        sources = ['https://github.com/ANHIG/IMGTHLA.git']
        scope = {'git_commit': IMGTHLA_COMMIT, 'format': 'original upstream data; not normalized allo manifest'}
        license_url = 'https://www.ebi.ac.uk/ipd/imgt/hla/licence/'
        atomic_json(root / 'reference/hla/IMGTHLA.acquisition.json', {'commit': IMGTHLA_COMMIT, 'release': release})
    elif resource == '1000g':
        release = '20190312_biallelic_SNV_and_INDEL / v2a_27022019 / GRCh38'
        base = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/'
        filename = 'ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
        panel = get_pair(base + filename, 'public_datasets/1000genomes/' + filename, '.tbi', lambda p: vcf_pair(p, sample='NA12878'))
        scope = extract_1000g(panel, root / 'samples/germline_vcf')
    elif resource in ('hcc1395', 'full-truth'):
        base = 'https://storage.googleapis.com/deepvariant/deepsomatic-case-studies/'
        if resource == 'hcc1395':
            base += 'quick-start/'
            directory = 'public_datasets/hcc1395_seqc2/'
            names = [s + ext for s in ('S1395_WGS_ilm_normal.bwa.dedup.chr1.quickstart', 'S1395_WGS_ilm_tumor.bwa.dedup.chr1.quickstart') for ext in ('.bam', '.bam.bai')]
            names += ['GRCh38_no_alts_chr1.fasta', 'GRCh38_no_alts_chr1.fasta.fai', 'SEQC2_truth.chr1.quick_start.vcf.gz', 'SEQC2_truth.chr1.quick_start.vcf.gz.tbi', 'SEQC2_truth.chr1.quick_start.bed']
            release = 'DeepSomatic quick-start named HCC1395/SEQC2 chr1 artifacts'
            scope = {'assembly': 'GRCh38', 'contig': 'chr1', 'scope': 'quick-start'}
        else:
            base += 'SEQC2-S1395-truth/'
            directory = 'public_datasets/hcc1395_seqc2/full_truth/'
            names = ['high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz', 'high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz.tbi', 'High-Confidence_Regions_v1.2.bed']
            release = 'SEQC2 VCF 1.2.1 / BED 1.2'
            scope = {'assembly': 'GRCh38', 'scope': 'full truth'}
        for filename in names:
            if filename.endswith(('.bam.bai', '.vcf.gz.tbi')): continue
            if filename.endswith('.bam'):
                get_pair(base + filename, directory + filename, '.bai', bam_pair)
            elif filename.endswith('.vcf.gz'):
                get_pair(base + filename, directory + filename, '.tbi', vcf_pair)
            else:
                get(base + filename, directory + filename)
        for p in list(files):
            if p.name.endswith('.bam'):
                bam_pair(p)
            elif p.name.endswith('.vcf.gz'):
                vcf_pair(p)
            elif p.name.endswith('.fasta'):
                ensure_fasta(p)
    elif resource == 'processed-rna':
        release = 'HCC1395 HISAT tutorial named tumor/normal BAM artifacts'
        scope = {'source': 'genomedata RNA-seq tutorial', 'release_published': False}
        base = 'https://genomedata.org/rnaseq-tutorial/results/alignments/hisat/'
        for name in ('HCC1395_tumor.bam', 'HCC1395_normal.bam'):
            get_pair(base + name, 'public_datasets/hcc1395_rna/' + name, '.bai', bam_pair)
        license_url = 'https://genomedata.org/rnaseq-tutorial/'
    else:
        license_url = 'https://www.ncbi.nlm.nih.gov/home/about/policies/'
        accession = 'SRR23456709' if resource == 'longread' else 'SRR37849525'
        directory = root / ('public_datasets/hcc1395_longread' if resource == 'longread' else 'public_datasets/hcc1395_full_rna/sra')
        files = [acquire_sra(accession, directory)]
        release = accession
        sources = ['https://www.ncbi.nlm.nih.gov/sra/' + accession]
        scope = {'accession': accession, 'format': 'original SRA; FASTQ/alignment generation separate'}
    artifacts = []
    for p in files:
        entry = {'local_relative_path': str(p.relative_to(root)), 'bytes': p.stat().st_size}
        if observed_sha256:
            digest = hashlib.sha256()
            with p.open('rb') as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(block)
            entry['observed_sha256'] = digest.hexdigest()
        artifacts.append(entry)
    versions = {}
    commands = ('samtools','bcftools') if resource in ('1000g','hcc1395','processed-rna','full-truth') else ('samtools',) if resource=='genome' else ('git',) if resource=='imgthla' else ('prefetch','vdb-validate') if resource in ('longread','full-rna') else ()
    for command in commands:
        if shutil.which(command):
            probe = subprocess.run([command,'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            versions[command] = probe.stdout.splitlines()[0] if probe.returncode==0 and probe.stdout.splitlines() else 'version probe unavailable'
    atomic_json(root / 'acquisition' / (resource + '.json'), {
        'schema_version': '1.0.0', 'resource_id': resource, 'release': release,
        'source_urls': sources, 'semantic_scope': scope, 'license_url': license_url,
        'license_status': 'consult upstream; no redistribution grant inferred', 'tool_versions': versions,
        'acquisition_mode': 'public_download', 'redistribution': 'consult selected upstream terms; not granted by this script',
        'artifacts': artifacts, 'acquisition_script_version': 'operational-1',
        'checksum_policy': 'observed optional; semantic/release and format validation required'})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dest', default=os.environ.get('NEOANTIMON2_DATA_ROOT'), help='chosen input/data directory (or NEOANTIMON2_DATA_ROOT)')
    parser.add_argument('--profile', choices=PROFILES)
    parser.add_argument('--resource', help='comma-separated selection: ' + ','.join(RESOURCES))
    parser.add_argument('--record-sha256', action='store_true', help='optional observed file hashes')
    # Shared operations for the shell workflows; these do not need --dest.
    parser.add_argument('--download', metavar='PUBLIC_URL')
    parser.add_argument('--output')
    parser.add_argument('--kind')
    parser.add_argument('--index-bam')
    args = parser.parse_args(argv)
    if args.index_bam:
        index_bam(args.index_bam)
        return
    if args.download:
        if not args.output:
            parser.error('--download requires --output')
        acquire(args.download, args.output, args.kind)
        return
    if not args.dest or not (args.profile or args.resource):
        parser.error('choose --dest and --profile minimum/full or --resource selection')
    selected = list(PROFILES.get(args.profile, ()))
    if args.resource:
        selected.extend(args.resource.split(','))
    if any(s not in RESOURCES for s in selected):
        parser.error('unknown resource selection')
    root = Path(args.dest).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    for resource in dict.fromkeys(selected):
        print('[resource]', resource, flush=True)
        acquire_resource(resource, root, args.record_sha256)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        sys.exit('[error] ' + str(error))
