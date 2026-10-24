# v0.3.0 release-candidate readiness

Version 0.3.0 is prepared locally; it is not a published release. Publication,
release tags and release assets still require explicit owner authorization.
The immutable v0.2.0 tag and assets must remain unchanged.

## Candidate scope and identity

The candidate consolidates readable documentation/cohort examples, localized
VCF validation corrections, four R resource/predictor setup APIs and the editable
scientific workflow overview. DESCRIPTION and CITATION.cff identify 0.3.0.
All sixteen exports, scientific schema/ID versions and existing analysis contracts
remain intact. There is no multi-sample-column VCF support or clinical validation.
Candidate commit/archive identities and test counts belong in the timestamped
release-candidate report; no unapproved tag is used as its identity.

## Required review evidence

- Full source and installed regression counts, R CMD check, clean-library install,
  literal Quick Start, four cohorts, external example, setup API/help inventory.
- Public VCF supported/rejected denominators, independent import assertions,
  unchanged semantic outputs/IDs and saved-project reopening across versions.
- Fresh representative reference download, verified cache reuse/checksums,
  explicit reference compatibility, clear missing-tool and licensing diagnostics.
- Shared catalog/runtime synchronization and existing Python/shell regression.
- Updated figure/source/exports, source-to-figure mapping, visual QA and links.
- No private history/files, datasets, run artifacts or release-only credentials.
- Untested real licensed predictors/MHCflurry environments and other platforms
  documented without treating synthetic fixtures as real prediction validation.

The prior production run used about 91 CPU minutes and 3.64 GB peak RSS; host
sleep affected wall time. It is an observation, not a reconstruction-accuracy
benchmark or an automatic release blocker. Retain verified inputs/configuration,
artifacts and reopening evidence instead of repeating it without a concrete reason.

## Only after owner release approval

1. Confirm the exact reviewed candidate commit and clean tree. Finalize release
   wording in README/getting-started/setup/API documentation, NEWS, CITATION.cff
   and installed help: remove unpublished-candidate wording only when true.
   Preserve intentional v0.2.0 historical references and fixture/schema versions.
2. Rebuild/check the exact release source and record the source-archive SHA256,
   R/dependency environment and provenance. Package build timestamps may vary;
   compare source content separately and identify the approved archive explicitly.
3. Install that archive into a clean library and run the installed examples/setup
   checks independently of a source checkout. Review all untested limitations.
4. Once authorized, create the approved immutable v0.3.0 tag at that exact commit.
   Verify real installation from the new GitHub tag and downloaded archive before
   publishing approved assets. Until then, do not advertise an existing v0.3.0 tag.
5. Publish only approved release notes, R source archive, SHA256 manifest and
   validation/provenance material; distinguish the R archive from GitHub's generated
   repository archives. Reverify remote tag/commit/assets after publication.

The completed workflow overview is in `docs/assets/neoantimon2-overview.svg`;
its editable source, exports, rebuild instructions and scientific mapping are
maintained in `docs/user/workflow-overview.md`. Explicit settings control current
project or bound longitudinal sample graphs. External export requires a completed
project, not native prediction scores/grades; derived views remain separate.
