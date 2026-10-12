#!/usr/bin/env bash
# Shared CLI for public-data development workflows. No user-home layout required.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${NEOANTIMON2_DATA_ROOT:-}"
RUN_ROOT="${NEOANTIMON2_RUN_ROOT:-}"
THREADS="${THREADS:-8}"
BAM_INPUT=""; FASTA_INPUT=""; GTF_INPUT=""; TRUTH_INPUT=""; SRA_INPUT=""; EXTRACT_ONLY=0
neo_workflow_args() {
  while (($#)); do
    case "$1" in
      --data-root|--run-root|--threads|--bam|--fasta|--gtf|--truth|--sra)
        if (($# < 2)); then echo "[error] missing value for $1" >&2; return 2; fi
        case "$1" in
          --data-root) ROOT="$2";; --run-root) RUN_ROOT="$2";; --threads) THREADS="$2";;
          --bam) BAM_INPUT="$2";; --fasta) FASTA_INPUT="$2";; --gtf) GTF_INPUT="$2";;
          --truth) TRUTH_INPUT="$2";; --sra) SRA_INPUT="$2";;
        esac
        shift 2;;
      --extract-only) EXTRACT_ONLY=1; shift;;
      --help|-h)
        echo "Usage: $0 --run-root DIR [--data-root DIR] [--bam FILE --fasta FILE --gtf FILE --truth FILE --sra FILE --threads N]"
        echo 'Explicit files may be anywhere. DATA_ROOT supplies convenience paths only. All derived output goes under RUN_ROOT.'
        exit 0;;
      *) echo "[error] unknown argument: $1" >&2; return 2;;
    esac
  done
  if [[ -z "$RUN_ROOT" ]]; then echo '[error] choose --run-root or NEOANTIMON2_RUN_ROOT' >&2; return 2; fi
  if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then echo '[error] threads must be positive' >&2; return 2; fi
  # Do not generate files under the input tree or the repository.
  python3 - "$RUN_ROOT" "$ROOT" "$SCRIPT_DIR/../.." <<'PY'
from pathlib import Path
import sys
run = Path(sys.argv[1]).expanduser().resolve()
for raw in sys.argv[2:]:
    if not raw:
        continue
    protected = Path(raw).expanduser().resolve()
    if run == protected or run in protected.parents or protected in run.parents:
        sys.exit('[error] run root overlaps input data or repository')
PY
  RUN_ROOT="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$RUN_ROOT")"
  if [[ -n "$ROOT" ]]; then ROOT="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$ROOT")"; fi
  mkdir -p "$RUN_ROOT"
  RUN_ROOT="$(cd -- "$RUN_ROOT" && pwd)"
  if [[ -n "$ROOT" ]]; then ROOT="$(cd -- "$ROOT" && pwd)"; fi
}
download() { python3 "$SCRIPT_DIR/download_public_data.py" --download "$1" --output "$2"; }
neo_index_bam() { python3 "$SCRIPT_DIR/download_public_data.py" --index-bam "$1"; }

neo_ensure_rna_bam() {
  if [[ -z "$BAM_INPUT" && ! -s "$BAM" ]]; then
    local options=(--run-root "$RUN_ROOT" --fasta "$CHR1_FA" --threads "$THREADS" --extract-only)
    if [[ -n "$ROOT" ]]; then options+=(--data-root "$ROOT"); fi
    if [[ -n "$SRA_INPUT" ]]; then options+=(--sra "$SRA_INPUT"); fi
    bash "$SCRIPT_DIR/extract_hcc1395_rna_chr1_from_sra.sh" "${options[@]}"
  fi
}
