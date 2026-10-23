# Component tests deliberately stop before transcript reconstruction/pMHC.
vcf_robust_path <- function(name) test_path('..', 'fixtures', 'vcf-compatibility', name)
vcf_robust_write <- function(rows = 'chrT\t2\t.\tA\tC\t.\tPASS\t.',
                             sample = NULL, meta = character(), header = NULL) {
  if (is.null(header)) header <- paste0('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    if (!is.null(sample)) paste0('\tFORMAT\t', sample) else '')
  path <- tempfile(fileext = '.vcf')
  writeLines(c('##fileformat=VCFv4.2', meta, header, rows), path)
  path
}
vcf_robust_parse <- function(path, role = 'somatic_vcf', sample = 'S1') {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile('vcf-runs-')),
    scientific = list(genome_assembly = 'GRCh38.p14'))
  p <- neo_project(samples = data.frame(sample_id = sample, patient_id = 'P1'), config = cfg,
    inputs = list(list(input_id = 'vcf', input_role = role, sample_id = sample, path = normalizePath(path, mustWork = TRUE))))
  fai <- read.delim(vcf_robust_path('reference.fa.fai'), header = FALSE,
    col.names = c('contig','length','offset','line_bases','line_width'))
  con <- file(vcf_robust_path('reference.fa'), 'rb'); on.exit(close(con))
  ref <- list(con = con, fai = fai, assembly = 'GRCh38.p14')
  Neoantimon2:::event_parse_vcf(p, p$inputs[[1L]], ref)
}
vcf_robust_source <- function(row) row$event_details_json$source_records[[1L]]

test_that('VCF component tests assert hand-derived boundary and normalization alleles', {
  # The final base is read independently from the committed hand-authored FASTA.
  seq <- readLines(vcf_robust_path('reference.fa'))[[2L]]
  n <- nchar(seq)
  rows <- vcf_robust_parse(vcf_robust_write(c('chrT\t1\t.\tC\tA\t.\tPASS\t.',
    paste('chrT',n,'.',substr(seq,n,n),'A','.','PASS','.',sep='\t'),
    'chrT\t3\t.\tA\tAA\t.\tPASS\t.', 'chrT\t2\t.\tAA\tA\t.\tPASS\t.',
    'chrT\t5\t.\tCG\tCT\t.\tPASS\t.', 'chrT\t5\t.\tCG\tTA\t.\tPASS\t.')))
  expect_length(rows, 6L)
  expect_identical(lapply(rows[3:6], function(x) c(x$pos,x$ref,x$alt)),
    list(c('1','C','CA'),c('1','CA','C'),c('6','G','T'),c('5','CG','TA')))
  expect_identical(vapply(rows, `[[`, '', 'event_type'), c('SNV','SNV','INS','DEL','SNV','MNV'))
  expect_true(all(vapply(rows, function(x) x$sample_id=='S1' && x$patient_id=='P1', TRUE)))
})

test_that('multiallelic Number=A/R projections preserve the selected allele, zero and missing', {
  meta <- c('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
    '##INFO=<ID=AD,Number=R,Type=Integer,Description="AD">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
    '##FORMAT=<ID=VF,Number=A,Type=Float,Description="VF">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">')
  rows <- vcf_robust_parse(vcf_robust_write('chrT\t2\tmulti\tA\tC,G\t.\tPASS\tAF=0,.;AD=10,0,.\tGT:AD:VF:PS\t1|2:10,0,.:0,.:22',
    sample='S1',meta=meta), role='germline_vcf')
  expect_identical(vapply(rows, `[[`, '', 'alt'), c('C','G'))
  for (i in seq_along(rows)) {
    s <- vcf_robust_source(rows[[i]])
    expected <- c('0','.')[i]
    expect_identical(s$selected_alt_index, as.integer(i))
    expect_identical(s$projected_info$AF, expected)
    expect_identical(s$projected_info$AD, list(ref='10',alt=expected))
    expect_identical(s$projected_format$AD, list(ref='10',alt=expected))
    expect_identical(s$projected_format$VF, expected)
    a <- Neoantimon2:::phase_source_assertions(rows[[i]],s)[[1L]]
    expect_identical(a$status,'INTERPRETABLE')
    expect_identical(a$block_label,'22')
    expect_identical(a$genotype_state,c('LEFT','RIGHT')[i])
  }
})

test_that('germline GT selection and somatic source assertions retain their separate contract', {
  for (gt in c('0/0','0/1','1/1','1/2','2|1')) {
    rows <- vcf_robust_parse(vcf_robust_write(paste0('chrT\t2\t.\tA\tC,G\t.\tPASS\t.\tGT\t',gt),sample='S1'),role='germline_vcf')
    expected <- switch(gt,'0/0'=character(),'0/1'='C','1/1'='C',c('C','G'))
    expect_identical(vapply(rows, `[[`, '', 'alt'), expected)
  }
  no_gt <- vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC,G\t.\tPASS\t.'))
  expect_length(no_gt,2L)
  for (gt in c('.','./.','0/.','1/3','99999999999999999999'))
    expect_error(vcf_robust_parse(vcf_robust_write(paste0('chrT\t2\t.\tA\tC,G\t.\tPASS\t.\tGT\t',gt),sample='S1'),role='germline_vcf'), 'germline GT|Germline GT')
  for (case in list(c('0/0','UNPHASED','ABSENT'),c('.','MISSING','UNKNOWN'),c('0/3','INVALID_INDEX','UNKNOWN'))) {
    rows <- vcf_robust_parse(vcf_robust_write(paste0('chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT\t',case[1]),sample='S1'))
    expect_length(rows,1L)
    a <- Neoantimon2:::phase_source_assertions(rows[[1]],vcf_robust_source(rows[[1]]))[[1L]]
    expect_identical(a$status,case[2]); expect_identical(a$genotype_state,case[3])
  }
})

test_that('valid omitted trailing FORMAT fields remain explicitly missing', {
  rows <- vcf_robust_parse(vcf_robust_path('short-format.vcf'))
  expect_length(rows,1L)
  s <- vcf_robust_source(rows[[1]])
  expect_identical(s$raw_format$sample,'0/1')
  expect_null(s$projected_format$AD)
  expect_identical(Neoantimon2:::phase_source_assertions(rows[[1]],s)[[1]]$status,'UNPHASED')
  # A missing entire sample field is not a numerical zero or a measured genotype.
  rows <- vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT:AD\t.',sample='S1'))
  s <- vcf_robust_source(rows[[1]])
  expect_null(s$projected_format$AD)
  expect_identical(Neoantimon2:::phase_source_assertions(rows[[1]],s)[[1]]$status,'MISSING')
})

test_that('invalid ALT delimiters and malformed A/R cardinalities are rejected without dropping assertions', {
  expect_error(vcf_robust_parse(vcf_robust_path('trailing-alt.vcf')), 'Malformed VCF ALT')
  meta <- c('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC,G\t.\tPASS\tAF=0.1',meta=meta)), 'cardinality')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC,G\t.\tPASS\t.\tAD\t10,2',meta=meta,sample='S1')), 'cardinality')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\tAF=0.1;AF=0.2',meta=meta)), 'duplicate VCF INFO')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\tDP=1;')), 'VCF INFO key')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT:GT\t0/1:1/1',sample='S1')), 'FORMAT')
})

test_that('header structure, POS spelling and authoritative sample bindings fail safely', {
  expect_error(vcf_robust_parse(vcf_robust_write(header='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT',rows='chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT')), 'VCF column header')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1',sample='OTHER')), 'sample does not match')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/1',sample='S1\tS2')), 'Multi-sample VCF')
  expect_error(vcf_robust_parse(vcf_robust_write(c('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    'chrT\t2\t.\tA\tC\t.\tPASS\t.'))), 'Duplicate VCF column header')
  for (pos in c('1.0','1e0','0','-1','2147483648'))
    expect_error(vcf_robust_parse(vcf_robust_write(paste('chrT',pos,'.','C','A','.','PASS','.',sep='\t'))),'VCF POS')
  expect_error(vcf_robust_parse(vcf_robust_write(header='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',rows='chrT\t2\t.\tA\tC\t.\tPASS\t.\tEXTRA')), 'Malformed VCF record')
})

test_that('unsupported sequence syntax and reference discrepancies are explicit rejections', {
  for (alt in c('*','<NON_REF>bad','N','.','A,,C'))
    expect_error(vcf_robust_parse(vcf_robust_write(paste('chrT',2,'.','A',alt,'.','PASS','.',sep='\t'))),'unsupported|Unsupported|Malformed')
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tC\tG\t.\tPASS\t.')),'REF mismatch')
  expect_error(vcf_robust_parse(vcf_robust_write('T\t2\t.\tA\tC\t.\tPASS\t.')),'Contig absent')
  expect_error(vcf_robust_parse(vcf_robust_write(meta='##reference=GRCh37')),'contradicts selected reference')
  expect_error(vcf_robust_parse(vcf_robust_write(meta='##contig=<ID=chrT,length=999>')),'contig length mismatch')
})

test_that('symbolic and breakend assertions are retained as unsupported reconstruction, not resolved peptides', {
  rows <- vcf_robust_parse(vcf_robust_write(c('chrT\t2\tsv\tA\t<DEL>\t.\tPASS\tEND=4',
    'chrT\t2\tbnd\tA\tA]chrT:5]\t.\tPASS\t.')))
  expect_identical(vapply(rows, `[[`, '', 'event_type'),c('SV','BND'))
  expect_true(all(vapply(rows,function(x)x$reconstruction_status=='UNSUPPORTED',TRUE)))
  expect_error(vcf_robust_parse(vcf_robust_write('chrT\t2\t.\tA\tA]absent:5]\t.\tPASS\t.')),'BND mate contig absent')
})

test_that('plain and gzip inputs preserve content and bad compression fails', {
  plain <- vcf_robust_write(); compressed <- tempfile(fileext='.vcf.gz')
  con <- gzfile(compressed,'wt'); writeLines(readLines(plain),con); close(con)
  a <- vcf_robust_parse(plain); b <- vcf_robust_parse(compressed)
  a[[1L]]$analysis_id <- b[[1L]]$analysis_id <- NULL
  expect_identical(a,b)
  mismatch <- tempfile(fileext='.vcf.gz'); file.copy(plain,mismatch)
  expect_error(vcf_robust_parse(mismatch),'Compression/extension mismatch')
  bytes <- readBin(compressed,'raw',n=file.info(compressed)$size)
  corrupt <- tempfile(fileext='.vcf.gz'); bytes[seq.int(length(bytes)-7L,length(bytes))] <- as.raw(0)
  writeBin(bytes,corrupt)
  expect_error(suppressWarnings(vcf_robust_parse(corrupt)),'compressed VCF|compressed data|compression|CRC|error reading')
})

test_that('INFO/FORMAT order, record order and duplicate occurrence provenance are stable', {
  meta <- c('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">','##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">','##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">')
  a <- vcf_robust_parse(vcf_robust_write('chrT\t2\tx\tA\tC\t.\tPASS\tAF=0.2;DP=10\tGT:PS:DP\t0|1:44:10',sample='S1',meta=meta))
  b <- vcf_robust_parse(vcf_robust_write('chrT\t2\tx\tA\tC\t.\tPASS\tDP=10;AF=0.2\tGT:DP:PS\t0|1:10:44',sample='S1',meta=meta))
  expect_identical(a[[1]]$event_id,b[[1]]$event_id)
  expect_identical(vcf_robust_source(a[[1]])$projected_info,vcf_robust_source(b[[1]])$projected_info)
  for (row in c(a,b)) expect_identical(Neoantimon2:::phase_source_assertions(row,vcf_robust_source(row))[[1]]$genotype_state,'RIGHT')
  dup <- vcf_robust_parse(vcf_robust_write(rep('chrT\t2\tx\tA\tC\t.\tPASS\t.',2)))
  merged <- Neoantimon2:::event_merge_rows(dup)
  expect_length(merged,1L)
  expect_length(merged[[1]]$event_details_json$source_records,2L)
})
