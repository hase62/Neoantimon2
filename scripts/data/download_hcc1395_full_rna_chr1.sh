#!/usr/bin/env bash
set -euo pipefail

source "$(dirname -- "${BASH_SOURCE[0]}")/workflow_paths.sh"
neo_workflow_args "$@"
ACC="SRR37849525"

PUBLIC="${RUN_ROOT}/development/HCC1395_full_rna"
SAMPLES="${ROOT}/samples"
HCC="${ROOT}/public_datasets/hcc1395_seqc2"
GENCODE="${ROOT}/reference/gencode"

SRA_DIR="${ROOT}/public_datasets/hcc1395_full_rna/sra"
FASTQ_DIR="${PUBLIC}/fastq"
STAR_INDEX="${PUBLIC}/star_chr1_index"
ALIGN_DIR="${PUBLIC}/star_chr1_alignment"
TMP_DIR="${PUBLIC}/tmp"
CHECK_DIR="${PUBLIC}/dna_rna_overlap"

CHR1_FA="${FASTA_INPUT:-${HCC}/GRCh38_no_alts_chr1.fasta}"
GTF_GZ="${GTF_INPUT:-${GENCODE}/gencode.v50.primary_assembly.annotation.gtf.gz}"
CHR1_GTF="${PUBLIC}/gencode.v50.chr1.gtf"

TRUTH_VCF="${TRUTH_INPUT:-${HCC}/SEQC2_truth.chr1.quick_start.vcf.gz}"

R1="${FASTQ_DIR}/${ACC}_1.fastq"
R2="${FASTQ_DIR}/${ACC}_2.fastq"

PREFIX="${ALIGN_DIR}/HCC1395.fullRNA.chr1."
RNA_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"

mkdir -p \
  "$FASTQ_DIR" \
  "$STAR_INDEX" \
  "$ALIGN_DIR" \
  "$TMP_DIR" \
  "$CHECK_DIR"

for cmd in prefetch fasterq-dump vdb-validate STAR samtools bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[error] required command not found: $cmd"
        exit 1
    fi
done

if [ ! -s "$CHR1_FA" ]; then
    echo "[error] chr1 FASTA not found: $CHR1_FA"
    exit 1
fi

if [ ! -s "$GTF_GZ" ]; then
    echo "[error] GENCODE GTF not found: $GTF_GZ"
    exit 1
fi

if [ ! -s "$TRUTH_VCF" ]; then
    echo "[error] SEQC2 truth VCF not found: $TRUTH_VCF"
    exit 1
fi

echo "=== HCC1395 full RNA-seq ==="
echo "SRA accession: $ACC"

SRA_FILE="${SRA_INPUT:-${SRA_DIR}/${ACC}/${ACC}.sra}"
if [[ ! -f "$SRA_FILE" ]]; then
    if [[ -z "$ROOT" || -n "$SRA_INPUT" ]]; then
        echo '[error] supply an existing --sra, or acquire full-rna with --data-root' >&2
        exit 1
    fi
    python3 "$SCRIPT_DIR/download_public_data.py" --dest "$ROOT" --resource full-rna
fi

echo "=== Validate SRA ==="

vdb-validate "$SRA_FILE"

FASTQ_STAGE="$(mktemp -d "${FASTQ_DIR}/.fasterq-XXXXXX")"
trap 'rm -rf -- "$FASTQ_STAGE"' EXIT
fasterq-dump "$SRA_FILE" --split-files --threads "$THREADS" --temp "$TMP_DIR" --outdir "$FASTQ_STAGE"
python3 - "$SCRIPT_DIR" "$FASTQ_STAGE/${ACC}_1.fastq" "$FASTQ_STAGE/${ACC}_2.fastq" <<'PYFASTQ'
import sys
sys.path.insert(0, sys.argv[1])
from download_public_data import validate
for path in sys.argv[2:]: validate(path, 'fastq')
PYFASTQ
mv "$FASTQ_STAGE/${ACC}_1.fastq" "$R1"
mv "$FASTQ_STAGE/${ACC}_2.fastq" "$R2"

# Derive from the selected original references/reads, without requiring any
# previous index, annotation subset or alignment from the run directory.
gzip -t "$GTF_GZ"
gzip -dc "$GTF_GZ" | awk 'BEGIN{FS=OFS="\t"} /^#/ {next} $1=="chr1" {print}' > "$CHR1_GTF"
[[ -s "$CHR1_GTF" ]] || { echo '[error] expected chr1 annotation absent' >&2; exit 1; }
STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" \
     --genomeFastaFiles "$CHR1_FA" --sjdbGTFfile "$CHR1_GTF" --sjdbOverhang 100 --genomeSAindexNbases 12
STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" --readFilesIn "$R1" "$R2" \
     --outFileNamePrefix "$PREFIX" --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes Standard --outSAMunmapped None --quantMode GeneCounts

echo "=== Index and validate RNA BAM ==="

neo_index_bam "$RNA_BAM"
samtools quickcheck -v "$RNA_BAM"

echo "=== RNA alignment summary ==="

samtools flagstat "$RNA_BAM" \
  > "${CHECK_DIR}/HCC1395.fullRNA.chr1.flagstat.txt"

samtools idxstats "$RNA_BAM" \
  > "${CHECK_DIR}/HCC1395.fullRNA.chr1.idxstats.txt"

cat "${CHECK_DIR}/HCC1395.fullRNA.chr1.flagstat.txt" | head -12

echo "=== Truth-site RNA coverage ==="

bcftools query \
  -f '%CHROM\t%POS0\t%POS\n' \
  "$TRUTH_VCF" \
  > "${CHECK_DIR}/SEQC2_truth.chr1.sites.bed"

samtools depth \
  -b "${CHECK_DIR}/SEQC2_truth.chr1.sites.bed" \
  "$RNA_BAM" \
  > "${CHECK_DIR}/SEQC2_truth.chr1.RNA_depth.tsv"

python3 - \
  "${CHECK_DIR}/SEQC2_truth.chr1.RNA_depth.tsv" \
  "$TRUTH_VCF" <<'PY'
import sys
import subprocess

depth_file = sys.argv[1]
vcf = sys.argv[2]

truth_n = int(
    subprocess.check_output(
        ["bcftools", "view", "-H", vcf],
        text=True
    ).count("\n")
)

covered = 0
depths = []

with open(depth_file) as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        dp = int(fields[2])
        depths.append(dp)
        if dp > 0:
            covered += 1

print(f"Truth variants: {truth_n}")
print(f"Truth positions with RNA coverage >0: {covered}")

if depths:
    print(f"Maximum RNA depth at truth positions: {max(depths)}")
PY

echo "=== Quick RNA variant calling at SEQC2 truth regions ==="

RNA_VCF="${CHECK_DIR}/HCC1395.fullRNA.chr1.calls.vcf.gz"

bcftools mpileup \
  -f "$CHR1_FA" \
  -R "$TRUTH_VCF" \
  -a FORMAT/AD,FORMAT/DP \
  -Ou \
  "$RNA_BAM" \
  | bcftools call \
      -Am \
      -Ou \
  | bcftools norm \
      -f "$CHR1_FA" \
      -m -any \
      -Oz \
      -o "$RNA_VCF"

bcftools index -f -t "$RNA_VCF"

RNA_OVERLAP="${CHECK_DIR}/SEQC2_truth_with_RNA_variant_support.vcf.gz"

bcftools isec \
  -n=2 \
  -w1 \
  -Oz \
  -o "$RNA_OVERLAP" \
  "$TRUTH_VCF" \
  "$RNA_VCF"

bcftools index -f -t "$RNA_OVERLAP"

echo "=== DNA truth / RNA variant overlap ==="

TRUTH_N="$(bcftools view -H "$TRUTH_VCF" | wc -l | tr -d ' ')"
RNA_CALL_N="$(bcftools view -H "$RNA_VCF" | wc -l | tr -d ' ')"
OVERLAP_N="$(bcftools view -H "$RNA_OVERLAP" | wc -l | tr -d ' ')"

echo "SEQC2 truth variants:           $TRUTH_N"
echo "RNA variant calls at regions:   $RNA_CALL_N"
echo "Exact REF/ALT overlap:          $OVERLAP_N"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  "$RNA_OVERLAP" \
  > "${CHECK_DIR}/SEQC2_truth_with_RNA_variant_support.tsv"

echo
echo "[ok] Full HCC1395 RNA-seq chr1 analysis complete."
echo
echo "RNA BAM:"
echo "  $RNA_BAM"
echo
echo "DNA/RNA overlap results:"
echo "  ${CHECK_DIR}/SEQC2_truth.chr1.RNA_depth.tsv"
echo "  ${CHECK_DIR}/SEQC2_truth_with_RNA_variant_support.tsv"
echo
echo "FASTQs retained for later whole-genome alignment:"
echo "  $R1"
echo "  $R2"
