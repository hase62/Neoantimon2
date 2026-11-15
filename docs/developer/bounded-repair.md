# Bounded molecular-condition and annotation repair

[User interpretation](../user/results-and-explain.md#molecular-conditions-and-native-support) ·
[Reproduction scope](../user/reproducibility.md#current-main-and-additional-replay-materials)

These repairs are integrated into current main through
[PR #9](https://github.com/hase62/Neoantimon2/pull/9), merge
`0fa9af969e3a0e812c2b4ff5efbbc9c1822d858a`. The immutable `v0.3.0` and
`reproducibility-v1` refs are unchanged. All 22 public APIs and the identity
constructors are preserved.

## Phase conditions and native decisions

`Origin.phase_status` remains a normal/germline-background label. Somatic
`phase_context` and exported peptide `provenance.phase_condition` retain the
complete molecular condition and `CONDITIONAL`/`RESOLVED` interpretation.
`Origin.origin_warning_codes` marks conditional support with `PHASE_UNCERTAIN`.
`T009_DETERMINISTIC` is conditional sequence determinism, not proof of all phase.

An explicitly conditional Origin cannot establish native molecular eligibility;
hard DNA failures remain failures. Candidate aggregation accepts another eligible
Origin, rejects all-ineligible support, and otherwise returns `not_evaluated`.
Vaccine admission consumes that per-Origin support. External `external_rank`
orders a selected measured channel and is not a native eligibility policy.
Missing legacy metadata is not relabelled; saved original results are not rewritten.

## Bounded extension and Sec provenance

Stop-loss extension admits only reached, explicitly cis somatic, length-preserving
substitutions under a complete standard nuclear extension. It does not enumerate
unknown combinations or add trans peers. General UTR/indel and complex/germline
extension cases remain guarded. Empty Sec positions yield no annotation IDs;
annotated exceptions and actually used exceptions remain separate.

## Source and regression references

- [Phase provenance](../../R/peptide_contract.R) and [Origin warnings](../../R/peptide_records.R).
- [Per-Origin eligibility](../../R/somatic_eligibility.R), [evidence linkage](../../R/somatic_evidence.R) and [vaccine admission](../../R/vaccine_contract.R).
- [Reached cis peers](../../R/reconstruction_core.R) and [Sec IDs](../../R/translation_spec.R).
- [Bounded regression controls](../../tests/testthat/test-bounded-repair.R) and [translation regression](../../tests/testthat/test-translation-t009.R).

Fresh baseline/current-source runs have separate source/run identities. Recorded
checks are described in [versioned validation](../development/validation-evidence.md#current-molecular-repair-integration);
public source tests and the historical eight-table snapshot do not imply that
every separate author-held replay archive has been published.
