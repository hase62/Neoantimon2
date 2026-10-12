#!/usr/bin/env bash
# Compatibility entry point. The standalone Python script is downloadable from GitHub.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in
  references) shift; set -- --resource genome,gencode,mane,imgthla "$@" ;;
  phased) shift; set -- --resource 1000g "$@" ;;
  somatic) shift; set -- --resource hcc1395 "$@" ;;
  processed-rna|longread) resource="$1"; shift; set -- --resource "$resource" "$@" ;;
  light) shift; set -- --resource genome,gencode,mane,imgthla,1000g,hcc1395 "$@" ;;
esac
exec python3 "$SCRIPT_DIR/download_public_data.py" "$@"
