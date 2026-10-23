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
RESOURCE_CATALOG_SHA256 = 'f2693c4edd0176376b7b1e078ffac8cedad5ecd31efe8076fe608f89cd9589e1'
_CATALOG = None
_CATALOG_PATH = None
_CHECK_ONLY = False
_NO_OVERWRITE = False


def resource_catalog():
    """One catalog for installed R, this utility and standalone downloads."""
    global _CATALOG
    if _CATALOG is not None:
        return _CATALOG
    candidates = ([Path(_CATALOG_PATH)] if _CATALOG_PATH else
                  [Path(__file__).resolve().parents[2] / 'inst/setup/resource-catalog.json',
                   Path(__file__).resolve().parent / 'resource-catalog.json'])
    local = next((p for p in candidates if p.is_file()), None)
    if local:
        payload = local.read_bytes()
    elif _CATALOG_PATH:
        raise ValueError('Resource catalog not found: ' + str(_CATALOG_PATH))
    else:
        url = 'https://raw.githubusercontent.com/hase62/Neoantimon2/main/inst/setup/resource-catalog.json'
        require('curl')
        payload = run('curl', '--location', '--fail', '--retry', '3', url, capture=True).encode('utf-8')
    if not _CATALOG_PATH and hashlib.sha256(payload).hexdigest() != RESOURCE_CATALOG_SHA256:
        raise ValueError('Resource catalog checksum mismatch; obtain the matching catalog/utility')
    catalog = json.loads(payload)
    if (catalog.get('schema_version') != 'neo-resource-catalog-1' or
            set(catalog['resources']) != set(RESOURCES) or
            {k: tuple(v) for k, v in catalog['profiles'].items()} != PROFILES):
        raise ValueError('Resource catalog selection contract mismatch')
    for spec in catalog['resources'].values():
        for item in spec['files']:
            path = Path(item['path'])
            if path.is_absolute() or '..' in path.parts:
                raise ValueError('Catalog path escapes destination')
    _CATALOG = catalog
    return catalog


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
    if _CHECK_ONLY and not output.exists():
        raise ValueError('Missing resource: ' + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    check = validator or (lambda p: validate(p, kind))
    if output.exists():
        try:
            check(output)
            print('[validated reuse]', output, flush=True)
            return output
        except (ValueError, OSError, subprocess.CalledProcessError) as error:
            if _CHECK_ONLY or _NO_OVERWRITE:
                raise ValueError('Invalid existing resource retained: ' + str(output)) from error
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
    if _CHECK_ONLY:
        pair_validator(output)
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and index.exists():
        try:
            pair_validator(output)
            print('[validated pair reuse]', output, flush=True)
            return output
        except (ValueError, OSError, subprocess.CalledProcessError):
            if _NO_OVERWRITE:
                raise ValueError('Invalid existing data/index pair retained: ' + str(output))
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
                    if _NO_OVERWRITE:
                        raise ValueError('Invalid existing pair member retained: ' + str(prior))
                    acquire(endpoint, candidate)
            else:
                acquire(endpoint, candidate)
        try:
            pair_validator(stage)
        except (ValueError, OSError, subprocess.CalledProcessError):
            if stage_index not in reused:
                raise
            if _NO_OVERWRITE:
                raise ValueError('Conflicting existing index retained: ' + str(index))
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
    if _CHECK_ONLY or (_NO_OVERWRITE and Path(str(p) + '.bai').exists()):
        bam_pair(p)
        return
    with tempfile.TemporaryDirectory(prefix='.neo-index-', dir=p.parent) as directory:
        stage = Path(directory) / (p.name + '.bai')
        run('samtools', 'index', p, stage)
        validate(stage, 'bai')
        os.replace(stage, str(p) + '.bai')
    bam_pair(p)


def ensure_fasta(path):
    p = validate(path, 'fasta')
    require('samtools')
    with tempfile.TemporaryDirectory(prefix='.neo-fasta-', dir=None if _CHECK_ONLY else p.parent) as directory:
        stage = Path(directory) / (p.name + '.fai')
        run('samtools', 'faidx', '-o', stage, p)
        validate(stage, 'fai')
        target = Path(str(p) + '.fai')
        if not target.exists() or target.read_bytes() != stage.read_bytes():
            if _CHECK_ONLY or (_NO_OVERWRITE and target.exists()):
                def geometry(file):
                    rows = [line.split('\t') for line in file.read_text().splitlines()]
                    return [(f[0], tuple(int(v) for v in f[1:])) for f in rows]
                if not target.exists() or geometry(target) != geometry(stage):
                    raise ValueError('Missing/inconsistent FASTA index retained: ' + str(target))
                return p
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
        if _CHECK_ONLY or (_NO_OVERWRITE and output.exists()):
            raise ValueError('Missing/inconsistent decompressed FASTA retained: ' + str(output))
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
    if _CHECK_ONLY and not destination.is_dir():
        raise ValueError('Missing subset directory: ' + str(destination))
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
                if _CHECK_ONLY or _NO_OVERWRITE:
                    raise ValueError('Invalid existing subset retained: ' + str(out))
        if _CHECK_ONLY:
            raise ValueError('Missing subset: ' + str(out))
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


def acquire_imgthla(destination, commit=IMGTHLA_COMMIT, release=IMGTHLA_RELEASE,
                    url='https://github.com/ANHIG/IMGTHLA.git'):
    require('git')
    destination = Path(destination)
    if _CHECK_ONLY and not destination.exists():
        raise ValueError('Missing IMGTHLA checkout: ' + str(destination))
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
            run('git', 'clone', '--no-checkout', url, stage)
            run('git', '-C', stage, 'checkout', '--detach', commit)
            check(stage)
            os.replace(stage, destination)
    return destination


def acquire_sra(accession, directory):
    require('prefetch', 'vdb-validate')
    directory = Path(directory)
    out = directory / accession / (accession + '.sra')
    if _CHECK_ONLY and not out.exists():
        raise ValueError('Missing SRA: ' + str(out))
    directory.mkdir(parents=True, exist_ok=True)
    if out.exists():
        try:
            run('vdb-validate', out)
            return out
        except subprocess.CalledProcessError:
            if _CHECK_ONLY or _NO_OVERWRITE:
                raise ValueError('Invalid existing SRA retained: ' + str(out))
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
    if _NO_OVERWRITE and path.exists():
        # Preserve earlier acquisition provenance under the opt-in R policy.
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False, prefix='.neo-metadata-') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stage = Path(stream.name)
    os.replace(stage, path)


def acquire_resource(resource, root, observed_sha256=False):
    root = Path(root)
    catalog = resource_catalog()
    spec = catalog['resources'][resource]
    files, sources = [], []
    downloads = [f for f in spec['files'] if f['action'] == 'download']
    def get(item):
        kind = 'fasta' if item['kind'] in ('transcript', 'protein') else item['kind']
        checker = (lambda p: fasta_scope(p, catalog['reference']['genome_contigs'])) if resource == 'genome' else (
            (lambda p: fasta_scope(p, {'chr1': 248956422}, exact=True))
            if item['path'].endswith('GRCh38_no_alts_chr1.fasta') else None)
        path = acquire(item['url'], root / item['path'], kind, validator=checker)
        files.append(path); sources.append(item['url'])
        return path
    scope = dict(spec['semantic_scope'])
    if resource == 'imgthla':
        checkout = spec['checkout']
        acquire_imgthla(root / checkout['path'], checkout['commit'], spec['release'], checkout['url'])
        sources = [checkout['url']]
        if not _CHECK_ONLY:
            atomic_json(root / 'reference/hla/IMGTHLA.acquisition.json',
                        {'commit': checkout['commit'], 'release': spec['release']})
    elif resource in ('longread', 'full-rna'):
        item = spec['files'][0]
        accession = spec['accession']
        files = [acquire_sra(accession, (root / item['path']).parent.parent)]
        sources = ['https://www.ncbi.nlm.nih.gov/sra/' + accession]
    else:
        for item in downloads:
            if item['kind'] in ('tbi', 'bai'):
                continue
            suffix = '.tbi' if item['kind'] == 'vcf' else '.bai' if item['kind'] == 'bam' else None
            if suffix:
                checker = bam_pair if suffix == '.bai' else (
                    (lambda p: vcf_pair(p, sample='NA12878')) if resource == '1000g' else vcf_pair)
                path = acquire_pair(item['url'], root / item['path'], suffix, checker)
                files.extend([path, Path(str(path) + suffix)])
                sources.extend([item['url'], item['url'] + suffix])
            else:
                get(item)
        if resource == 'genome':
            target = next(f for f in spec['files'] if f['action'] == 'gunzip')
            path = decompress_fasta(root / target['parent'], root / target['path'])
            files.extend([path, Path(str(path) + '.fai')])
        elif resource == '1000g':
            scope = extract_1000g(root / downloads[0]['path'], root / 'samples/germline_vcf')
        elif resource == 'hcc1395':
            for path in files:
                if path.name.endswith('.fasta'):
                    ensure_fasta(path)
    artifacts = []
    for path in files:
        entry = {'local_relative_path': str(path.relative_to(root)), 'bytes': path.stat().st_size}
        if observed_sha256:
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(block)
            entry['observed_sha256'] = digest.hexdigest()
        artifacts.append(entry)
    versions = {}
    commands = ('samtools','bcftools') if resource in ('1000g','hcc1395','processed-rna','full-truth') else ('samtools',) if resource=='genome' else ('git',) if resource=='imgthla' else ('prefetch','vdb-validate') if resource in ('longread','full-rna') else ()
    for name in commands:
        if shutil.which(name):
            probe = subprocess.run([name,'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            versions[name] = probe.stdout.splitlines()[0] if probe.returncode==0 and probe.stdout.splitlines() else 'version probe unavailable'
    result = {
        'schema_version': '1.0.0', 'resource_id': resource, 'release': spec['release'],
        'source_urls': sources, 'semantic_scope': scope, 'license_url': spec['license_url'],
        'license_status': 'consult upstream; no redistribution grant inferred', 'tool_versions': versions,
        'acquisition_mode': 'public_download', 'redistribution': 'consult selected upstream terms; not granted by this script',
        'artifacts': artifacts, 'acquisition_script_version': 'operational-1',
        'checksum_policy': 'observed optional; semantic/release and format validation required'}
    if not _CHECK_ONLY:
        atomic_json(root / 'acquisition' / (resource + '.json'), result)
    return result


def main(argv=None):
    global _CATALOG_PATH, _CATALOG, _CHECK_ONLY, _NO_OVERWRITE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dest', default=os.environ.get('NEOANTIMON2_DATA_ROOT'), help='chosen input/data directory (or NEOANTIMON2_DATA_ROOT)')
    parser.add_argument('--profile', choices=PROFILES)
    parser.add_argument('--resource', help='comma-separated selection: ' + ','.join(RESOURCES))
    parser.add_argument('--record-sha256', action='store_true', help='optional observed file hashes')
    parser.add_argument('--catalog', help='explicit matching resource catalog (normally located automatically)')
    parser.add_argument('--check-only', action='store_true', help='validate selected existing resources without acquisition or promotion')
    parser.add_argument('--no-overwrite', action='store_true', help='retain invalid existing files rather than replacing them')
    # Shared operations for the shell workflows; these do not need --dest.
    parser.add_argument('--download', metavar='PUBLIC_URL')
    parser.add_argument('--output')
    parser.add_argument('--kind')
    parser.add_argument('--index-bam')
    args = parser.parse_args(argv)
    _CATALOG_PATH, _CATALOG = args.catalog, None
    _CHECK_ONLY, _NO_OVERWRITE = args.check_only, args.no_overwrite
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
    if not _CHECK_ONLY:
        root.mkdir(parents=True, exist_ok=True)
    for resource in dict.fromkeys(selected):
        print('[resource]', resource, flush=True)
        acquire_resource(resource, root, args.record_sha256)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        sys.exit('[error] ' + str(error))
