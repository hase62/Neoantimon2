# Neoantimon2 instructions for Codex

For a requested analysis, read [shared instructions](docs/ai/AI_INSTRUCTIONS.md),
then follow the single authoritative [analysis playbook](docs/ai/ANALYSIS_PLAYBOOK.md).
Use the [manifest example](docs/ai/INPUT_MANIFEST.example.yaml) and
[request examples](docs/ai/EXAMPLE_REQUESTS.md) as needed. Inspect existing inputs
before asking questions; record the installed version/commit and preserve original data.
Current main has 22 APIs; frozen reproducibility-v1 has 20 and immutable v0.3.0 has 16.
Use the version supporting the requested feature; never infer it from old logs.

For repository maintenance, follow the user's requested scope. The playbook is
for analysis, not a requirement to run an analysis during documentation work.
Do not change scientific algorithms, APIs or release tags to simplify an example.
