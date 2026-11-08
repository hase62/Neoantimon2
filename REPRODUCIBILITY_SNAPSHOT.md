# Neoantimon2 Reproducibility Snapshot v1

A versioned research reproducibility snapshot with supporting data and verification scripts.

- Tag: `reproducibility-v1` (20 NAMESPACE exports; S3 methods/imports excluded).
- Scientific baseline: `19152001c964bc2078c6d9b163365866f6a8b9d7`.
- This tag's docs-only child retains all baseline non-documentation bytes,
  including R source, NAMESPACE, DESCRIPTION, tests, examples and schemas.
- Release: https://github.com/hase62/Neoantimon2/releases/tag/reproducibility-v1
- Named asset: `Neoantimon2_Reproducibility_v1.zip` and its `.sha256` companion.

Download the named asset and extract into a fresh folder. From its root:

```sh
python3 -I scripts/verify_package.py
python3 -I scripts/verify_snapshot.py
```

The eight saved input TSVs reaggregate to 15 unique mutant sequences, 75 exact
peptide–HLA candidates, 36 Origins, 180 candidate–Origin links and 300 external
measurements. AGRN AAPCRVLPQG: 5 Origins, 5 HLA candidates, 25 links, 20 measurements.
These are preserved technical results, not independently measured immunogenicity.
The archive checksums and scientific source hashes are included in provenance.

Replay does not run molecular reconstruction, model inference, oracle cases,
pVACseq or full HCC1395 processing. The original 15-sequence selector executable
remains unrecovered; saved membership is not a complete VCF-to-pilot replay.
Historical qualifications, including the original driver/controller termination,
remain unchanged. No submission destination or manuscript DOI is asserted.

Core/Converter retain PolyForm Noncommercial and notices. Independent archive
helpers are MIT; author-owned package explanations are CC BY 4.0, with explicit
scope and attribution. Scientific data and third-party terms are unchanged.
No model weights, vendor binaries or bulk references are redistributed.

For latest software use main (22 APIs); for the released baseline use v0.3.0
(16 APIs). These refs are distinct despite the shared development version number
on this snapshot and main. Record the exact commit in every reproduction.
