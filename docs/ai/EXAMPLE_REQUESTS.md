# Requests for Codex and Claude Code

Both entrypoints route to the same [playbook](ANALYSIS_PLAYBOOK.md). These prompts
are intended usage examples, not evidence of a separate agent-session test.
The documented R commands and declarations are tested; availability and behavior
of an independently launched Codex/Claude Code session remain environment-specific.

**Codex — offline:** “Read AGENTS.md and run Mode A of the analysis playbook.
Keep outputs local and explain NOT_RUN and NOT_GRADED.”

**Claude Code — offline:** “Read CLAUDE.md and follow Mode A of the shared
Neoantimon2 analysis playbook. Do not install predictors or download references.”

**Either — my data:** “Read the repository agent entrypoint and analyze the
samples defined in `/absolute/path/manifest.yaml`. Check the configuration and
sample/HLA bindings, show a concise plan, and use a new run directory. Ask before
expensive execution or downloads. Keep all identifiable results local.”

**Multiple patients:** “Use Mode D for my explicit manifest. P1 has S1 and S2;
P2 has S3. Preserve separate sample occurrences and each patient's supplied HLA.
Do not enable longitudinal analysis.”

**External:** “Use Mode E with my completed project and reviewed BigMHC-IM
output/converter. Keep IM evidence and external ranking separate from EL/native
results. Do not claim predictor execution if the inputs are synthetic.”

**Explain only:** “Open `/absolute/path/completed/run` and explain the candidate
IDs I supply. Retain all origins and evidence states. Do not rerun the pipeline.”

See [the manifest template](INPUT_MANIFEST.example.yaml). Replace its placeholders;
YAML is an agent workflow declaration, not a direct Neoantimon2 package input.
The default repository-helper command checks it without executing analysis:

```sh
Rscript --vanilla docs/ai/run-analysis.R /absolute/path/manifest.yaml
```

No prompt can supply missing HLA/assembly, dependencies, permissions or reference
files. The agent must report missing essentials instead of guessing them.
