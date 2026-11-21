#!/usr/bin/env python3
"""Small, fixed-model MHCflurry BA software verification, not calibration validity.

Uses an explicitly supplied existing installation/model. No automatic downloads.
Original public-reference model identity cannot be inferred from this experiment.
"""
import argparse
import csv
import datetime
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys

for key in ['TF_NUM_INTEROP_THREADS', 'TF_NUM_INTRAOP_THREADS', 'OMP_NUM_THREADS']:
    os.environ[key] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


def sha(path):
    with Path(path).open('rb') as stream:
        # Python 3.10 is supported by the fixed MHCflurry 2.0.6 environment.
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
        return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ['model-dir', 'output-dir', 'converter']:
        parser.add_argument('--'+arg, required=True)
    args = parser.parse_args()
    target = Path(args.output_dir); target.mkdir(parents=True, exist_ok=False)
    model = Path(args.model_dir).resolve()
    records = [{'path': str(p.relative_to(model)), 'sha256': sha(p)}
               for p in sorted(model.rglob('*')) if p.is_file()]
    identity = hashlib.sha256(json.dumps(records, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if identity != '4bf7e97baaaa53b4d90ec3cbc9f2710fc3ca1da77a7df6363544142b7e09b06d':
        raise ValueError('Fixed reference model files do not match the recorded model identity')
    converter_sha = sha(args.converter)
    if converter_sha != 'cfd9573ec3a48d1f3bd29e663b51e33eff2b19191d677efba224dd8b636ef5b5':
        raise ValueError('Pinned MHCflurry BA converter checksum mismatch')
    version = importlib.metadata.version('mhcflurry')
    if version != '2.0.6':
        raise ValueError('This reproducible verification pins MHCflurry 2.0.6')
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    peptides = [''.join(alphabet[b % 20] for b in hashlib.sha256(
        ('Neoantimon2-reference-software-v1:'+str(i)).encode()).digest()[:9]) for i in range(256)]
    queries = ['ACDEFGHIK', 'LMNPQRSTV']
    (target/'reference-peptides.txt').write_text('\n'.join(peptides)+'\n')
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load(str(model))
    values = predictor.predict(peptides + queries, allele='HLA-A*02:01')
    native_percentiles = predictor.percentile_ranks(values[256:], allele='HLA-A*02:01')
    reference = target/'scores.tsv'
    reference.write_text('score\n'+''.join(format(float(v),'.17g')+'\n' for v in values[:256]))
    original_reference = target/'reference-native.csv'
    with original_reference.open('w', newline='') as stream:
        writer = csv.writer(stream);writer.writerow(['peptide','allele','mhcflurry_affinity'])
        for peptide,value in zip(peptides,values[:256]):
            writer.writerow([peptide,'HLA-A*02:01',format(float(value),'.17g')])
    # Preserve native output before conversion and the original converter output.
    raw = target/'queries.csv'
    with raw.open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['peptide','allele','mhcflurry_affinity','mhcflurry_affinity_percentile'])
        for peptide,value,percentile in zip(queries, values[256:],native_percentiles):
            writer.writerow([peptide,'HLA-A*02:01',format(float(value),'.17g'),format(float(percentile),'.17g')])
    spec = importlib.util.spec_from_file_location('pinned_converter', args.converter)
    converter = importlib.util.module_from_spec(spec); spec.loader.exec_module(converter)
    normalized = converter.convert(raw.read_text(), converter.SPEC, {'HLA-A*02:01'})
    def table(path, rows):
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=converter.COLUMNS, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    table(target/'converter-original.tsv', normalized)
    # An explicit receipt from THIS run, not an inferred historical model version.
    for row in normalized:
        row['model_version'] = 'weights-sha256:'+identity
    table(target/'predictions.tsv', normalized)
    record = dict(reference_id='fixed-mhcflurry-ba-software-v1',provider='mhcflurry-ba',
        predictor_version=version,model_version='weights-sha256:'+identity,task='mhc1_binding',
        metric='affinity_nM',unit='nM',score_direction='lower_is_better',mhc_key='HLA-A*02:01',
        peptide_length=9,species='NCBITaxon:9606',hla_class='I',reference_set=
        '256 deterministic SHA256-derived AA20 9mers; software test only, not a representative biological population',
        source='https://github.com/openvax/mhcflurry/releases/download/pre-2.0/models_class1_pan.selected.20200610.tar.bz2',
        source_version='e9263e842da394abd9b41e2418c58ce8e6662fb8 / models_class1_pan.selected.20200610',
        original_file='reference-native.csv',original_sha256=sha(original_reference),
        reference_peptides_file='reference-peptides.txt',reference_peptides_sha256=sha(target/'reference-peptides.txt'),
        scores_file='scores.tsv',scores_sha256=sha(reference),
        count=256,license='MHCflurry Apache-2.0; review upstream model/data terms',model_provenance='receipt.json')
    (target/'manifest.json').write_text(json.dumps({'schema_version':'reference-percentiles-1.0.0','distributions':[record]},indent=2)+'\n')
    receipt = dict(purpose='same-fixed-model software verification, not validation of pVACtools model equivalence or biological calibration',
        completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        script_sha256=sha(__file__),execution_settings={key:os.environ[key] for key in
            ['TF_NUM_INTEROP_THREADS','TF_NUM_INTRAOP_THREADS','OMP_NUM_THREADS','CUDA_VISIBLE_DEVICES','TF_CPP_MIN_LOG_LEVEL']},
        model_files=records,model_sha256=identity,predictor_version=version,python=sys.version,platform=platform.platform(),
        packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()
                  if d.metadata['Name'].lower() in ['tensorflow-macos','keras','h5py','numpy','mhcflurry']},
        converter_sha256=converter_sha,provider='mhcflurry-ba',model_binding='Only model_version field added to converter-original.tsv from observed weights hash',
        reference_count=256,query_count=2,outputs={p.name:sha(p) for p in sorted(target.iterdir()) if p.is_file()},
        expected_percentiles=[100*sum(float(r)<float(q) for r in values[:256])/256+
                              50*sum(float(r)==float(q) for r in values[:256])/256 for q in values[256:]])
    (target/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','model_identity':identity,'query_scores':[float(v) for v in values[256:]],
                      'expected_percentiles':receipt['expected_percentiles']}),flush=True)


if __name__ == '__main__':
    main()
