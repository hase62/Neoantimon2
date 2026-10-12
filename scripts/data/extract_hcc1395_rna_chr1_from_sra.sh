#!/usr/bin/env bash
set -euo pipefail

source "$(dirname -- "${BASH_SOURCE[0]}")/workflow_paths.sh"
neo_workflow_args "$@"

ACC="SRR37849525"

PUBLIC="${ROOT}/public_datasets/hcc1395_full_rna"
SRA_FILE="${SRA_INPUT:-${PUBLIC}/sra/${ACC}/${ACC}.sra}"

HCC="${ROOT}/public_datasets/hcc1395_seqc2"
TRUTH_VCF="${TRUTH_INPUT:-${HCC}/SEQC2_truth.chr1.quick_start.vcf.gz}"
CHR1_FA="${FASTA_INPUT:-${HCC}/GRCh38_no_alts_chr1.fasta}"

OUT="${RUN_ROOT}/development/HCC1395_full_rna/sra_chr1_alignment"
TMP="${OUT}/tmp"
CHECK="${OUT}/dna_rna_overlap"


mkdir -p "$OUT" "$TMP" "$CHECK"

for cmd in sam-dump samtools bcftools python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[error] command not found: $cmd"
        exit 1
    fi
done

if [ ! -s "$SRA_FILE" ]; then
    echo "[error] SRA file not found:"
    echo "  $SRA_FILE"
    exit 1
fi

if [[ "$EXTRACT_ONLY" == 0 && ! -s "$TRUTH_VCF" ]]; then
    echo "[error] truth VCF not found:"
    echo "  $TRUTH_VCF"
    exit 1
fi

if [ ! -s "$CHR1_FA" ]; then
    echo "[error] chr1 FASTA not found:"
    echo "  $CHR1_FA"
    exit 1
fi

echo "=== Inspect SRA alignment reference ==="

HEADER="${OUT}/${ACC}.header.sam"

# Dump a tiny region only to recover the SAM header.
# Try common GRCh38 chr1 names.
REFNAME=""

for candidate in chr1 1 NC_000001.11; do
    echo "[probe] ${candidate}"

    if sam-dump \
        --header \
        --primary \
        --aligned-region "${candidate}:1-100000" \
        "$SRA_FILE" \
        > "$HEADER" 2>/dev/null
    then
        if grep -q '^@SQ' "$HEADER"; then
            if grep -Fq "SN:${candidate}" "$HEADER"; then
                REFNAME="$candidate"
                break
            fi
        fi
    fi
done

# Fallback: identify chromosome 1 by GRCh38 length.
if [ -z "$REFNAME" ]; then
    REFNAME="$(
        awk -F '\t' '
        $1=="@SQ" {
            sn="";
            ln="";
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SN:/) {
                    sn=substr($i,4)
                }
                if ($i ~ /^LN:/) {
                    ln=substr($i,4)
                }
            }
            if (ln=="248956422") {
                print sn;
                exit
            }
        }' "$HEADER"
    )"
fi

if [ -z "$REFNAME" ]; then
    echo "[error] Could not identify GRCh38 chromosome 1 in SRA alignment."
    echo
    echo "Reference header:"
    grep '^@SQ' "$HEADER" | head -30 || true
    exit 1
fi

echo "[ok] chromosome 1 reference in SRA: $REFNAME"

CHR1_LEN="$(
    awk -F '\t' -v target="$REFNAME" '
    $1=="@SQ" {
        sn="";
        ln="";
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) sn=substr($i,4)
            if ($i ~ /^LN:/) ln=substr($i,4)
        }
        if (sn==target) {
            print ln;
            exit
        }
    }' "$HEADER"
)"

echo "[info] chromosome 1 length: $CHR1_LEN"

if [ "$CHR1_LEN" != "248956422" ]; then
    echo "[error] chromosome 1 length does not match GRCh38."
    echo "Expected: 248956422"
    echo "Observed: $CHR1_LEN"
    exit 1
fi

echo "=== Extract aligned RNA reads from chromosome 1 ==="

RAW_BAM="${OUT}/HCC1395.${ACC}.chr1.raw.bam"
SORTED_BAM="${OUT}/HCC1395.${ACC}.chr1.sorted.bam"

rm -f "$RAW_BAM"

if sam-dump --help 2>&1 | grep -q -- '--rna-splicing'; then
    echo "[info] sam-dump supports --rna-splicing"

    sam-dump \
        --header \
        --primary \
        --rna-splicing \
        --aligned-region "$REFNAME" \
        "$SRA_FILE" \
        | samtools view \
            -@ "$THREADS" \
            -b \
            -o "$RAW_BAM" \
            -
else
    echo "[info] sam-dump does not expose --rna-splicing; using stored alignment"

    sam-dump \
        --header \
        --primary \
        --aligned-region "$REFNAME" \
        "$SRA_FILE" \
        | samtools view \
            -@ "$THREADS" \
            -b \
            -o "$RAW_BAM" \
            -
fi

echo "=== Normalize chr1 reference name to chr1 ==="

if [ "$REFNAME" != "chr1" ]; then
    samtools view -H "$RAW_BAM" \
      | awk -v old="$REFNAME" '
          BEGIN { OFS="\t" }
          {
              if ($1=="@SQ") {
                  for (i=2; i<=NF; i++) {
                      if ($i=="SN:" old) {
                          $i="SN:chr1"
                      }
                  }
              }
              print
          }' \
      > "${OUT}/header.chr1.sam"

    samtools reheader \
      "${OUT}/header.chr1.sam" \
      "$RAW_BAM" \
      > "${OUT}/HCC1395.${ACC}.chr1.renamed.bam"

    RAW_FOR_SORT="${OUT}/HCC1395.${ACC}.chr1.renamed.bam"
else
    RAW_FOR_SORT="$RAW_BAM"
fi

# Rebuild from the original SRA on every invocation; an old sorted BAM is
# never an input requirement or an unbound cache from a different source.
samtools sort -@ "$THREADS" -T "${TMP}/HCC1395_chr1_sort" -o "$SORTED_BAM" "$RAW_FOR_SORT"
neo_index_bam "$SORTED_BAM"

samtools quickcheck -v "$SORTED_BAM"

if [[ "$EXTRACT_ONLY" == 1 ]]; then echo "[ok] generated RNA BAM: $SORTED_BAM"; exit 0; fi

echo "=== Alignment summary ==="

samtools flagstat "$SORTED_BAM" \
    > "${CHECK}/HCC1395.${ACC}.chr1.flagstat.txt"

samtools idxstats "$SORTED_BAM" \
    > "${CHECK}/HCC1395.${ACC}.chr1.idxstats.txt"

head -12 "${CHECK}/HCC1395.${ACC}.chr1.flagstat.txt"

echo "=== Create canonical symlink ==="

echo "=== DNA truth positions with RNA coverage ==="

TRUTH_BED="${CHECK}/SEQC2_truth.chr1.sites.bed"

bcftools query \
    -f '%CHROM\t%POS0\t%POS\n' \
    "$TRUTH_VCF" \
    > "$TRUTH_BED"

samtools depth \
    -a \
    -b "$TRUTH_BED" \
    "$SORTED_BAM" \
    > "${CHECK}/SEQC2_truth.chr1.RNA_depth.tsv"

python3 - \
    "${CHECK}/SEQC2_truth.chr1.RNA_depth.tsv" \
    "$TRUTH_VCF" <<'PY'
import subprocess
import sys

depth_file = sys.argv[1]
truth_vcf = sys.argv[2]

truth = subprocess.check_output(
    ["bcftools", "view", "-H", truth_vcf],
    text=True,
)

truth_n = len([x for x in truth.splitlines() if x])

covered = 0
depths = []

with open(depth_file) as fh:
    for line in fh:
        f = line.rstrip().split("\t")
        if len(f) < 3:
            continue
        dp = int(f[2])
        depths.append(dp)
        if dp > 0:
            covered += 1

print(f"Truth variants: {truth_n}")
print(f"Truth positions with RNA coverage >0: {covered}")

if depths:
    print(f"Maximum RNA depth at truth positions: {max(depths)}")
PY

echo "=== RNA variant calling at SEQC2 truth sites ==="

RNA_VCF="${CHECK}/HCC1395.${ACC}.RNA_calls_at_truth_sites.vcf.gz"

bcftools mpileup \
    -f "$CHR1_FA" \
    -R "$TRUTH_VCF" \
    -a FORMAT/AD,FORMAT/DP \
    -Ou \
    "$SORTED_BAM" \
    | bcftools call \
        -Am \
        -Ou \
    | bcftools norm \
        -f "$CHR1_FA" \
        -m -any \
        -Oz \
        -o "$RNA_VCF"

bcftools index -f -t "$RNA_VCF"

OVERLAP_VCF="${CHECK}/SEQC2_truth_with_RNA_variant_support.vcf.gz"

bcftools isec \
    -n=2 \
    -w1 \
    -Oz \
    -o "$OVERLAP_VCF" \
    "$TRUTH_VCF" \
    "$RNA_VCF"

bcftools index -f -t "$OVERLAP_VCF"

TRUTH_N="$(
    bcftools view -H "$TRUTH_VCF" \
      | wc -l \
      | tr -d ' '
)"

RNA_N="$(
    bcftools view -H "$RNA_VCF" \
      | wc -l \
      | tr -d ' '
)"

OVERLAP_N="$(
    bcftools view -H "$OVERLAP_VCF" \
      | wc -l \
      | tr -d ' '
)"

echo
echo "=== DNA / RNA summary ==="
echo "SEQC2 truth variants:         $TRUTH_N"
echo "RNA calls at truth regions:   $RNA_N"
echo "Exact REF/ALT overlap:        $OVERLAP_N"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "$OVERLAP_VCF" \
    > "${CHECK}/SEQC2_truth_with_RNA_variant_support.tsv"

echo
echo "[ok] HCC1395 chr1 RNA extraction complete."
echo
echo "RNA BAM:"
echo "  $SORTED_BAM"
echo
echo "RNA depth:"
echo "  ${CHECK}/SEQC2_truth.chr1.RNA_depth.tsv"
echo
echo "DNA/RNA variant overlap:"
echo "  ${CHECK}/SEQC2_truth_with_RNA_variant_support.tsv"
