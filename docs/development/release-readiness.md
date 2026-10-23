# Consolidated release preparation

This is a local release-validation checklist, not a published release. Version
0.3.0 is provisional and needs owner approval. Keep v0.2.0 tags/assets immutable.

## Scope and current source identity

The consolidated source includes the readable README/cohort examples, VCF
validation corrections, and four R resource/predictor setup APIs. Existing
analysis contracts and scientific limitations remain in force. There is no new
multi-sample-column VCF support, predictor accuracy claim or clinical validation.
The source DESCRIPTION still reports 0.2.0; record the source commit to distinguish
it from the immutable release. The published setup baseline is
97c13211b9fb211b6a26439efc7a32a8941dadb0.

## Before a new version is approved

- Review the exact release-validation commit, full source/installed test counts,
  public VCF outcomes, clean-library install, example results and unresolved issues.
- Preserve the catalog/runtime synchronization check and both Python test suites.
- Review licensing, citation, reference identities, manual predictor installation
  requirements, and tested/untested platform boundaries.
- Confirm no private history, credentials, datasets or run artifacts enter source.
- Confirm upstream main and immutable v0.2.0 tag/release asset identities.

## After explicit version/publication authorization

1. Set the approved DESCRIPTION version. Review version assertions, CITATION.cff,
   README/getting-started/setup/API/resources documentation, `man/neo_setup.Rd`
   development-only wording and archive names. Search for `0.2.0`, `Phase 3` and
   `development build` to distinguish current availability from historical release
   references; do not replace these strings indiscriminately. Change only
   metadata that identifies the new release; retain the immutable release examples
   where useful. Keep historical fixture/provider schema versions unchanged.
2. Make a short release description covering usability, parser corrections and
   setup APIs, plus actual limitations. A lengthy historical changelog is unnecessary.
3. Rebuild from the exact candidate commit and run full regression, asset sync,
   Python tests and R CMD check --no-manual. Retain logs and exact denominators.
4. Install the candidate source archive in a clean library with explicit dependencies,
   run installed Quick Start/cohorts/external example and setup diagnostics from an
   unrelated working directory, and record SHA256 and source identity.
5. Once approved, create a new immutable tag at that exact commit. Test installation
   using the actual new GitHub tag and the downloaded release archive; compare APIs,
   packaged resources and example semantics before publishing release assets.
6. Attach the R source archive and its SHA256 manifest, concise release notes,
   validation report/summary and license/notice/citation information as appropriate.
   GitHub's generated repository archive and R CMD build archive are distinct;
   identify each explicitly. Reverify remote commit/tag/assets after publication.

No new tag, version publication or release asset is created by this checklist.
Real licensed predictor execution, a real MHCflurry environment/models and other
platforms remain separate validation areas unless actually tested and documented.

## Facts for the later workflow diagram

Inputs bind samples to patients, HLA and explicit resources. Native analysis links
events to reconstruction/translation origins, altered peptides and HLA candidates;
configured native prediction/evidence feed their supported ranking/results/explain
views. Export Bundle preserves all candidate-origin links; external prediction
imports attach by exact peptide/full-presenter and produce independent derived
results/explanations. A completed original project is required for external import.
Longitudinal, vaccine, receptor/presentation and indirect/allo evidence have their
own implemented connections; external scores do not automatically feed them.
Missing evidence, NOT_RUN, NOT_GRADED, evaluated zero and negative evidence remain
distinct. The final graphical asset belongs to Phase 5.
