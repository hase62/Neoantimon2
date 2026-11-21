#!/usr/bin/env python3
"""Compare independent counts with the exact pinned pVACtools 7.1.5 methods.

Inputs are explicit freshly verified local files, never a moving branch download.
Only the two unmodified methods are compiled from the supplied upstream source.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import h5py
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hdf5', required=True); p.add_argument('--source', required=True)
    p.add_argument('--output-dir', required=True); args=p.parse_args()
    source=Path(args.source); hdf5=Path(args.hdf5)
    observed=hashlib.sha256(hdf5.read_bytes()).hexdigest()
    if observed!='2caf69fca14bf70b1e6b69a38f544708139b88243dd0b1e72d190369f5b78b54':
        raise ValueError('Pinned reference file checksum mismatch')
    source_bytes = source.read_bytes()
    source_sha = hashlib.sha256(source_bytes).hexdigest()
    if source_sha != '5af398a9e44f14c9c70942278848ad9bb50dce68d54f82329d4492ff1f0e3541':
        raise ValueError('Pinned pVACtools source checksum mismatch')
    functions=[n for n in ast.walk(ast.parse(source_bytes)) if isinstance(n,ast.FunctionDef)
               and n.name in ('calculate_normalized_percentile','_normalize_allele')]
    if len(functions)!=2: raise ValueError('Expected exactly two upstream methods')
    cls=ast.ClassDef(name='PinnedParser', bases=[], keywords=[], body=functions, decorator_list=[])
    namespace={'np':np,'h5py':h5py,'os':os}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[cls],type_ignores=[])),str(source),'exec'),namespace)
    parser=namespace['PinnedParser']();parser.reference_scores_path=str(hdf5.parent);parser.reference_scores={}
    target=Path(args.output_dir);target.mkdir(parents=True,exist_ok=False)
    results=[]
    with h5py.File(hdf5,'r') as f:
        for method,score,reversed_,per_length,mixed in [('MHCflurry',500,False,5.068,2.017),('BigMHC_EL',0.5,True,0.617,0.267)]:
            values=f[method+'/9mer'][...]
            if len(values)!=100000 or not np.isfinite(values).all(): raise ValueError('Invalid distribution')
            better=int(np.count_nonzero(values>score if reversed_ else values<score));ties=int(np.count_nonzero(values==score))
            independent=100*(better+0.5*ties)/len(values)
            pvac=parser.calculate_normalized_percentile('HLA-A*02:01',9,score,method,is_reversed=reversed_,mode='per_length')
            pooled=parser.calculate_normalized_percentile('HLA-A*02:01',9,score,method,is_reversed=reversed_)
            if round(independent,3)!=pvac or pvac!=per_length or pooled!=mixed: raise ValueError('Checkpoint mismatch')
            file=target/(method+'-9mer.tsv');file.write_text('score\n'+''.join(format(float(v),'.17g')+'\n' for v in values))
            results.append(dict(method=method,score=score,count=len(values),strictly_better=better,ties=ties,
                independent_unrounded=independent,pvac_per_length=pvac,pvac_default_mixed=pooled,
                score_file_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),dataset_attributes=dict(f[method+'/9mer'].attrs)))
    (target/'receipt.json').write_text(json.dumps(dict(status='PASS',reference_sha256=observed,
        pvac_commit='871f2b3d979d5d1988af7a2d3ef2d52e179d1524',reference_commit='24fdb2aed1c575d1421eae909a473706edda84dc',
        source_sha256=source_sha,model_identity='UNKNOWN: computation check only, not applied to historical scores',
        methods=results),indent=2)+'\n')
    print(json.dumps(results),flush=True)


if __name__=='__main__':main()
