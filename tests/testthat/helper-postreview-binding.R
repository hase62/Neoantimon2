postreview_binding_project <- function(conflicting_expression=FALSE) {
  root <- tempfile("postreview-binding-");dir.create(root)
  tool <- file.path(root,"mock-netMHCpan");calls <- file.path(root,"calls.txt")
  writeLines(c("#!/bin/sh", "if [ \"$1\" = '-h' ]; then echo 'NetMHCpan version 4.2'; exit 0; fi",
    paste("printf 'synthetic-call\\n' >>",shQuote(calls)),
    "input=''; allele=''; while [ $# -gt 0 ]; do",
    "case \"$1\" in -p) input=$2;shift 2;; -a) allele=$2;shift 2;; *) shift;; esac;done",
    "echo 'Pos HLA Peptide Aff(nM) %Rank_BA'", "n=0;while IFS= read -r p;do n=$((n+1));",
    "printf '%s %s %s 120 0.5\\n' \"$n\" \"$allele\" \"$p\";done < \"$input\""),tool)
  Sys.chmod(tool,"0755")
  vcf <- file.path(root,"tumor.vcf")
  writeLines(c("##fileformat=VCFv4.2","##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:80,20"),vcf)
  depth <- file.path(root,"depth.tsv")
  writeLines(c("contig\tpos\tref\talt\ttumor_depth","chrW\t4\tA\tG\t100"),depth)
  expression <- file.path(root,"expression.tsv")
  writeLines(c("transcript_id\tgene_id\tgene_tpm","TW.1\tGW.1\t3",
    if(conflicting_expression) "TW\tGW\t9" else character()),expression)
  ccf <- file.path(root,"ccf.tsv")
  writeLines(c("contig\tpos\tref\talt\tccf\tmethod","chrW\t4\tA\tG\t0.95\tsynthetic supplied CCF"),ccf)
  base <- reconstruction_test_project();cfg <- base$config
  cfg$peptides <- peptide_test_policy(8L,9L)
  cfg$predictors <- list(list(provider="netmhcpan",task="mhc1_binding",executable=tool,options=list()))
  cfg$somatic_evidence <- Neoantimon2:::somatic_evidence_policy(list(sources=list(
    list(kind="dna_reads",sample_id="sample-1",path=depth,reference_build="GRCh38.p14"),
    list(kind="expression",sample_id="sample-1",path=expression),
    list(kind="ccf",sample_id="sample-1",path=ccf,reference_build="GRCh38.p14"))))
  p <- neo_run(neo_project(samples=base$samples,config=cfg,
    hla=list("patient-1"=c("HLA-A*02:01","HLA-B*07:02")),
    inputs=list(list(input_id="postreview-binding",input_role="somatic_vcf",sample_id="sample-1",path=vcf))))
  list(p=p,calls=calls,tool=tool,inputs=c(vcf,depth,expression,ccf))
}

