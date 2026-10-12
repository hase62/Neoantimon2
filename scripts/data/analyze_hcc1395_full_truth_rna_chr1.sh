#!/usr/bin/env bash
set -euo pipefail

source "$(dirname -- "${BASH_SOURCE[0]}")/workflow_paths.sh"
neo_workflow_args "$@"
DATA_ROOT="$ROOT"

SOURCE="${DATA_ROOT}/public_datasets/hcc1395_seqc2/full_truth"
RNA_SOURCE="${RUN_ROOT}/development/HCC1395_full_rna/sra_chr1_alignment"

RUN="${RUN_ROOT}/development/HCC1395_dna_rna_chr1_full_truth"
INTERMEDIATE="${RUN}/intermediate"
RESULTS="${RUN}/results"
LOGS="${RUN}/logs"

mkdir -p "$INTERMEDIATE" "$RESULTS" "$LOGS"

BAM="${BAM_INPUT:-${RNA_SOURCE}/HCC1395.SRR37849525.chr1.sorted.bam}"
CHR1_FA="${FASTA_INPUT:-${DATA_ROOT}/public_datasets/hcc1395_seqc2/GRCh38_no_alts_chr1.fasta}"

BASE="https://storage.googleapis.com/deepvariant/deepsomatic-case-studies/SEQC2-S1395-truth"

FULL_VCF="${TRUTH_INPUT:-${SOURCE}/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
FULL_TBI="${FULL_VCF}.tbi"
FULL_BED="${SOURCE}/High-Confidence_Regions_v1.2.bed"

CHR1_VCF="${INTERMEDIATE}/SEQC2_high-confidence.chr1.vcf.gz"
CHR1_NORM="${INTERMEDIATE}/SEQC2_high-confidence.chr1.biallelic.norm.vcf.gz"

ALLELES="${INTERMEDIATE}/SEQC2_chr1.alleles.tsv.gz"
EVIDENCE_VCF="${INTERMEDIATE}/HCC1395_RNA_evidence_at_truth_sites.vcf.gz"
EVIDENCE_RAW="${INTERMEDIATE}/HCC1395_RNA_evidence_raw.tsv"

FINAL="${RESULTS}/SEQC2_chr1_RNA_evidence.tsv"
SUMMARY="${RESULTS}/summary.txt"

for cmd in curl bcftools samtools bgzip tabix python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[error] missing command: $cmd"
        exit 1
    fi
done

neo_ensure_rna_bam

if [ ! -s "$BAM" ]; then
    echo "[error] RNA BAM missing: $BAM"
    exit 1
fi

if [ ! -s "${BAM}.bai" ]; then
    echo "[error] RNA BAM index missing: ${BAM}.bai"
    exit 1
fi

if [ ! -s "$CHR1_FA" ]; then
    echo "[error] chr1 FASTA missing: $CHR1_FA"
    exit 1
fi

echo "=== Download full SEQC2 truth set ==="

if [[ -z "$TRUTH_INPUT" ]]; then
    if [[ -z "$DATA_ROOT" ]]; then echo '[error] choose --data-root or --truth' >&2; exit 1; fi
    python3 "$SCRIPT_DIR/download_public_data.py" --dest "$DATA_ROOT" --resource full-truth
fi

echo "=== Extract chromosome 1 truth ==="

bcftools view \
    -r chr1 \
    -Oz \
    -o "$CHR1_VCF" \
    "$FULL_VCF"

bcftools index -f -t "$CHR1_VCF"

CHR1_ORIGINAL_N="$(
    bcftools view -H "$CHR1_VCF" | wc -l | tr -d ' '
)"

echo "Original chr1 truth records: $CHR1_ORIGINAL_N"

echo "=== Normalize to biallelic truth records ==="

bcftools norm \
    -f "$CHR1_FA" \
    -m -any \
    -Oz \
    -o "$CHR1_NORM" \
    "$CHR1_VCF"

bcftools index -f -t "$CHR1_NORM"

CHR1_NORM_N="$(
    bcftools view -H "$CHR1_NORM" | wc -l | tr -d ' '
)"

echo "Normalized biallelic truth records: $CHR1_NORM_N"

echo "=== Build constrained allele target ==="

bcftools query \
    -f '%CHROM\t%POS\t%REF,%ALT\n' \
    "$CHR1_NORM" \
    | bgzip -c \
    > "$ALLELES"

tabix -f -s1 -b2 -e2 "$ALLELES"

echo "=== Measure RNA evidence at known truth alleles ==="

bcftools mpileup \
    -f "$CHR1_FA" \
    -R "$CHR1_NORM" \
    -d 1000000 \
    -a FORMAT/AD,FORMAT/DP \
    -Ou \
    "$BAM" \
    | bcftools call \
        -m \
        -C alleles \
        -T "$ALLELES" \
        -A \
        -Oz \
        -o "$EVIDENCE_VCF"

bcftools index -f -t "$EVIDENCE_VCF"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' \
    "$EVIDENCE_VCF" \
    > "$EVIDENCE_RAW"

echo "=== Build complete truth-site evidence table ==="

python3 - "$CHR1_NORM" "$EVIDENCE_RAW" "$FINAL" <<'PY'
import csv
import subprocess
import sys

truth_vcf = sys.argv[1]
evidence_tsv = sys.argv[2]
out_tsv = sys.argv[3]

truth_text = subprocess.check_output(
    [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
        truth_vcf
    ],
    text=True
)

truth = []
for line in truth_text.splitlines():
    if not line:
        continue
    chrom, pos, ref, alt = line.split("\t")
    truth.append((chrom, pos, ref, alt))

evidence = {}

with open(evidence_tsv, newline="") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if len(row) < 6:
            continue

        chrom, pos, ref, alt, dp, ad = row[:6]
        key = (chrom, pos, ref, alt)

        try:
            dp_i = int(dp) if dp not in ("", ".") else 0
        except ValueError:
            dp_i = 0

        ref_count = 0
        alt_count = 0

        if ad not in ("", "."):
            parts = ad.split(",")

            if len(parts) >= 1 and parts[0] not in ("", "."):
                try:
                    ref_count = int(parts[0])
                except ValueError:
                    ref_count = 0

            if len(parts) >= 2 and parts[1] not in ("", "."):
                try:
                    alt_count = int(parts[1])
                except ValueError:
                    alt_count = 0

        evidence[key] = (dp_i, ref_count, alt_count)

with open(out_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
        "rna_dp",
        "rna_ref_count",
        "rna_alt_count",
        "rna_vaf",
        "rna_support_status"
    ])

    for chrom, pos, ref, alt in truth:
        key = (chrom, pos, ref, alt)

        dp, ref_count, alt_count = evidence.get(
            key,
            (0, 0, 0)
        )

        denom = ref_count + alt_count

        if denom > 0:
            vaf = alt_count / denom
            vaf_text = f"{vaf:.6f}"
        else:
            vaf_text = "NA"

        if len(ref) == 1 and len(alt) == 1:
            variant_type = "SNV"
        elif len(ref) == len(alt):
            variant_type = "MNV"
        elif len(ref) < len(alt):
            variant_type = "INS"
        else:
            variant_type = "DEL"

        if alt_count >= 1:
            status = "RNA_ALT_SUPPORTED"
        elif dp > 0:
            status = "RNA_COVERED_NO_ALT"
        else:
            status = "RNA_NO_COVERAGE"

        writer.writerow([
            chrom,
            pos,
            ref,
            alt,
            variant_type,
            dp,
            ref_count,
            alt_count,
            vaf_text,
            status
        ])
PY

echo "=== Summarize ==="

python3 - "$FINAL" "$SUMMARY" "$CHR1_ORIGINAL_N" "$CHR1_NORM_N" <<'PY'
import csv
import sys
from collections import Counter

path = sys.argv[1]
summary_path = sys.argv[2]
original_n = int(sys.argv[3])
normalized_n = int(sys.argv[4])

rows = []

with open(path, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)

status = Counter(r["rna_support_status"] for r in rows)
types = Counter(r["variant_type"] for r in rows)

alt_ge_1 = sum(int(r["rna_alt_count"]) >= 1 for r in rows)
alt_ge_2 = sum(int(r["rna_alt_count"]) >= 2 for r in rows)
alt_ge_3 = sum(int(r["rna_alt_count"]) >= 3 for r in rows)
alt_ge_5 = sum(int(r["rna_alt_count"]) >= 5 for r in rows)

covered = sum(int(r["rna_dp"]) > 0 for r in rows)

lines = [
    "=== HCC1395 DNA truth / RNA evidence summary ===",
    f"Original SEQC2 chr1 truth records:       {original_n}",
    f"Normalized biallelic truth records:      {normalized_n}",
    f"Truth alleles with RNA coverage >0:      {covered}",
    f"Truth alleles with RNA ALT >=1:          {alt_ge_1}",
    f"Truth alleles with RNA ALT >=2:          {alt_ge_2}",
    f"Truth alleles with RNA ALT >=3:          {alt_ge_3}",
    f"Truth alleles with RNA ALT >=5:          {alt_ge_5}",
    "",
    "RNA status:",
]

for key in [
    "RNA_ALT_SUPPORTED",
    "RNA_COVERED_NO_ALT",
    "RNA_NO_COVERAGE",
]:
    lines.append(f"  {key}: {status.get(key, 0)}")

lines.append("")
lines.append("Variant type:")

for key in sorted(types):
    lines.append(f"  {key}: {types[key]}")

text = "\n".join(lines) + "\n"

print(text, end="")

with open(summary_path, "w") as out:
    out.write(text)
PY

echo
echo "[ok] Full SEQC2 chr1 truth × HCC1395 RNA analysis complete."
echo
echo "Source truth:"
echo "  $CHR1_NORM"
echo
echo "Final evidence:"
echo "  $FINAL"
echo
echo "Summary:"
echo "  $SUMMARY"
