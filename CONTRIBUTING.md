# Contributing

Report reproducible defects through the project's GitHub Issues after public
publication. Include the package version/commit, operating system, R version,
minimal synthetic input, exact command, and diagnostic output. Do not upload
patient information, credentials, restricted reference data, or licensed tools.
For scientific disagreements, identify the input, expected interpretation, and
supporting evidence; do not silently redefine missingness, biological identities,
or a frozen contract to make a test pass.

Pull requests should explain the problem, change, and relevant validation.
Keep scientific source, schemas, fixtures, and reproducibility evidence intact.
Run the relevant tests and package checks described in the
[reproducibility guide](docs/user/reproducibility.md). The existing
[PolyForm Noncommercial license](LICENSE) and [notice](NOTICE) remain controlling;
contribution review does not authorize relicensing third-party material.

Submit converter definitions and OK/NG feedback to the separate
[Community Converter Registry](https://www.hase62.jp/Neoantimon2-Converters/),
not this package. Feedback is community experience, not official validation.
Identify the converter ID/version, predictor version, and a safe synthetic
reproduction when reporting a broken converter.
