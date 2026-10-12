#!/usr/bin/env bash
set -euo pipefail

source "$(dirname -- "${BASH_SOURCE[0]}")/workflow_paths.sh"
neo_workflow_args "$@"

PUBLIC="${ROOT}/public_datasets/hcc1395_full_rna"
OUT="${RUN_ROOT}/development/HCC1395_dna_rna_quickstart"
CHECK="${OUT}/dna_rna_overlap"

BAM="${BAM_INPUT:-${RUN_ROOT}/development/HCC1395_full_rna/sra_chr1_alignment/HCC1395.SRR37849525.chr1.sorted.bam}"

HCC="${ROOT}/public_datasets/hcc1395_seqc2"
TRUTH_VCF="${TRUTH_INPUT:-${HCC}/SEQC2_truth.chr1.quick_start.vcf.gz}"
CHR1_FA="${FASTA_INPUT:-${HCC}/GRCh38_no_alts_chr1.fasta}"

mkdir -p "$CHECK"

neo_ensure_rna_bam

for f in "$BAM" "${BAM}.bai" "$TRUTH_VCF" "$CHR1_FA"; do
    if [ ! -s "$f" ]; then
        echo "[error] missing file: $f"
        exit 1
    fi
done

echo "=== Basic RNA BAM summary ==="

samtools flagstat "$BAM" \
  > "${CHECK}/HCC1395.chr1.flagstat.txt"

samtools idxstats "$BAM" \
  > "${CHECK}/HCC1395.chr1.idxstats.txt"

head -12 "${CHECK}/HCC1395.chr1.flagstat.txt"

echo "=== Count SEQC2 truth variants ==="

TRUTH_N="$(
    bcftools view -H "$TRUTH_VCF" \
      | wc -l \
      | tr -d ' '
)"

echo "SEQC2 truth variants: $TRUTH_N"

echo "=== RNA pileup at truth variants ==="

bcftools mpileup \
  -f "$CHR1_FA" \
  -R "$TRUTH_VCF" \
  -a FORMAT/AD,FORMAT/DP \
  -Ou \
  "$BAM" \
  | bcftools call \
      -Am \
      -Ou \
  | bcftools norm \
      -f "$CHR1_FA" \
      -m -any \
      -Oz \
      -o "${CHECK}/HCC1395.RNA_calls_at_truth_sites.vcf.gz"

bcftools index -f -t \
  "${CHECK}/HCC1395.RNA_calls_at_truth_sites.vcf.gz"

echo "=== Exact DNA truth / RNA ALT overlap ==="

bcftools isec \
  -n=2 \
  -w1 \
  -Oz \
  -o "${CHECK}/SEQC2_truth_with_RNA_ALT_support.vcf.gz" \
  "$TRUTH_VCF" \
  "${CHECK}/HCC1395.RNA_calls_at_truth_sites.vcf.gz"

bcftools index -f -t \
  "${CHECK}/SEQC2_truth_with_RNA_ALT_support.vcf.gz"

RNA_CALL_N="$(
    bcftools view -H \
      "${CHECK}/HCC1395.RNA_calls_at_truth_sites.vcf.gz" \
      | wc -l \
      | tr -d ' '
)"

OVERLAP_N="$(
    bcftools view -H \
      "${CHECK}/SEQC2_truth_with_RNA_ALT_support.vcf.gz" \
      | wc -l \
      | tr -d ' '
)"

echo "=== Build truth-site RNA depth table ==="

bcftools query \
  -f '%CHROM\t%POS0\t%POS\n' \
  "$TRUTH_VCF" \
  > "${CHECK}/SEQC2_truth_sites.bed"

samtools depth \
  -a \
  -b "${CHECK}/SEQC2_truth_sites.bed" \
  "$BAM" \
  > "${CHECK}/SEQC2_truth_RNA_depth.tsv"

COVERED_N="$(
    awk '$3 > 0 {n++} END {print n+0}' \
      "${CHECK}/SEQC2_truth_RNA_depth.tsv"
)"

echo
echo "=== FINAL SUMMARY ==="
echo "SEQC2 truth variants:             $TRUTH_N"
echo "Truth sites with RNA coverage:    $COVERED_N"
echo "RNA variant calls at truth sites: $RNA_CALL_N"
echo "Exact DNA REF/ALT RNA overlap:    $OVERLAP_N"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  "${CHECK}/SEQC2_truth_with_RNA_ALT_support.vcf.gz" \
  > "${CHECK}/SEQC2_truth_with_RNA_ALT_support.tsv"

echo
echo "[ok] Results:"
echo "  ${CHECK}/SEQC2_truth_RNA_depth.tsv"
echo "  ${CHECK}/SEQC2_truth_with_RNA_ALT_support.tsv"
