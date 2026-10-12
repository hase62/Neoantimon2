peptide_test_contexts <- function(rows, transcript = NULL) {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  old <- Neoantimon2:::reconstruction_contexts(p, b, rows)
  new <- Neoantimon2:::translation_contexts(p, b, rows, old)
  if (!is.null(transcript)) new <- Filter(function(x) x$transcript_id == transcript, new)
  new
}

peptide_test_policy <- function(i = 3L, ii = 5L)
  Neoantimon2:::peptide_policy(list(class_i_lengths = i, class_ii_lengths = ii))

peptide_test_sequences <- function(context, i = 3L, ii = 3L)
  vapply(Neoantimon2:::peptide_enumerate(context, peptide_test_policy(i, ii))$peptides,
    `[[`, "", "peptide_sequence")

peptide_test_long <- function(position = 12L) {
  normal <- "MACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
  replacement <- if (substring(normal, position, position) == "A") "V" else "A"
  mutant <- normal; substring(mutant, position, position) <- replacement
  interval <- list(normal_aa_start = position, normal_aa_end = position,
    mutant_aa_start = position, mutant_aa_end = position,
    normal_sequence = substring(normal, position, position), mutant_sequence = replacement,
    left_anchor = c(position - 1L, position - 1L), right_anchor = c(position + 1L, position + 1L))
  list(reconstruction_id = "reconstruction:long", translation = list(
    translation_id = "translation:long", status = "EVALUABLE", mutant_protein = mutant,
    patient_normal_protein = normal, reference_protein = normal,
    normal_tumor_difference_intervals = list(interval),
    reference_tumor_difference_intervals = list(interval), frame_effect = "unchanged",
    protein_effect = "substitution", termination_effect = "unchanged",
    mutant_boundary = "STOP_CODON"))
}

peptide_test_run <- function(sec = FALSE) {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    if (sec) "chrS2\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1" else
      "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"), vcf)
  base <- reconstruction_test_project()
  cfg <- base$config
  cfg$peptides <- peptide_test_policy(if (sec) 8L else 3L, if (sec) 9L else 5L)
  p <- neo_project(samples = base$samples, config = cfg,
    hla = list("patient-1" = c("HLA-A*02:01", "HLA-B*07:02",
      "HLA-DRA*01:01~HLA-DRB1*04:01", "HLA-DQA1*05:01~HLA-DQB1*02:01",
      "HLA-DPA1*01:03~HLA-DPB1*04:01")),
    inputs = list(list(input_id = "peptide-fixture", input_role = "somatic_vcf",
      sample_id = "sample-1", path = vcf)))
  neo_run(p)
}
