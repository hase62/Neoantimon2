# Installed setup assets

`resource-catalog.json` is the canonical public resource catalog for R and
Python. It preserves the existing resource IDs, selections and reference roles.
URLs, declared releases, terms and selection-specific tools are defined here.

`download_public_data.py` is an exact generated runtime copy of the maintained
`scripts/data/download_public_data.py`, intentionally included so installed R
packages can acquire advanced resources without a development checkout.
Edit the maintained utility, then run `scripts/development/sync_setup_assets.py`.
The `--check` mode and offline acquisition tests detect unsynchronized assets.
The catalog hash is bound into the standalone helper; explicit catalogs are an
advanced trusted-input option, not automatic support for new profiles.

`check_mhcflurry.py` only inspects/imports the chosen Python runtime. It never
installs packages or models. Older MHCflurry releases need TensorFlow even when
their package metadata omits it; newer releases use PyTorch. Prediction/model
readiness requires an explicit synthetic smoke check.
