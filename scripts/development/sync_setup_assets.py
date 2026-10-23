#!/usr/bin/env python3
"""Synchronize the installed runtime utility; resource-catalog.json is canonical."""
import argparse
import hashlib
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check', action='store_true')
args = parser.parse_args()
catalog = root / 'inst/setup/resource-catalog.json'
source = root / 'scripts/data/download_public_data.py'
runtime = root / 'inst/setup/download_public_data.py'
sha = hashlib.sha256(catalog.read_bytes()).hexdigest()
content = re.sub(r"RESOURCE_CATALOG_SHA256 = '[^']+'",
                 "RESOURCE_CATALOG_SHA256 = '" + sha + "'", source.read_text())
if args.check:
    if source.read_text() != content or not runtime.exists() or runtime.read_text() != content:
        raise SystemExit('Setup assets differ; run scripts/development/sync_setup_assets.py')
else:
    source.write_text(content)
    runtime.write_text(content)
print('Setup assets synchronized; catalog SHA256:', sha)
