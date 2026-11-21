#!/usr/bin/env python3
"""Explicit offline HDF5 -> one-column local reference. Requires optional h5py.

Never downloads data, infers a model version, mixes lengths, or runs in Core.
The caller must provide a verified model identity and its primary record.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['hdf5', 'expected-sha256', 'dataset', 'output-dir', 'reference-id',
                 'provider', 'predictor-version', 'model-version', 'model-provenance',
                 'task', 'metric', 'unit', 'score-direction', 'hla', 'reference-set',
                 'source', 'source-version', 'license']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    if args.model_version.strip().lower() in ('', 'unknown', 'unversioned', 'na'):
        parser.error('Model weights identity must be verified; HDF5 algorithm names are insufficient')
    if sha(args.hdf5) != args.expected_sha256:
        parser.error('HDF5 checksum mismatch')
    import h5py
    import numpy as np
    with h5py.File(args.hdf5, 'r') as stream:
        dataset = stream[args.dataset]
        if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 1:
            parser.error('Select exactly one one-dimensional length-specific dataset')
        values = dataset[...]
    length = args.dataset.rsplit('/', 1)[-1]
    if not length.endswith('mer') or not length[:-3].isdigit():
        parser.error('Dataset must explicitly identify peptide length, e.g. MHCflurry/9mer')
    if not len(values) or not np.isfinite(values).all():
        parser.error('Empty or non-finite reference scores')
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=False)
    scores = target / 'scores.tsv'
    scores.write_text('score\n' + ''.join(format(float(v), '.17g') + '\n' for v in values))
    record = {k: getattr(args, k) for k in ['reference_id', 'provider', 'predictor_version',
        'model_version', 'task', 'metric', 'unit', 'score_direction', 'reference_set',
        'source', 'source_version', 'license']}
    record.update(species='NCBITaxon:9606', hla_class='I', mhc_key=args.hla,
        peptide_length=int(length[:-3]), original_sha256=sha(args.hdf5), scores_file='scores.tsv',
        scores_sha256=sha(scores), count=len(values), model_provenance=args.model_provenance,
        source_dataset=args.dataset)
    (target / 'manifest.json').write_text(json.dumps({
        'schema_version': 'reference-percentiles-1.0.0', 'distributions': [record]}, indent=2)+'\n')


if __name__ == '__main__':
    main()
