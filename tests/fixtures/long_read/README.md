# Synthetic T028–T030 fixtures

No fixture is a patient observation. Reference resources contain the frozen
synthetic reconstruction reference plus chrLR/chrLRminus and chrIntergenic.
The added GRCh38-labelled miniature coordinates are test data, not real GRCh38
sequence. Real GENCODE validation is a separate script over authoritative files.

LRREF/LRREF2 deliberately share structure to preserve ambiguous authoritative
start hypotheses. LRCAT1/LRCAT2 contribute different catalog junctions. The
models.gtf cases cover FSM, catalog combination, donor, acceptor, retention,
novel exon, exon skip, terminal changes, lost start, frameshift, single exon,
minus strand, antisense, intergenic and a locus without compatible reference.

sqanti3.gtf/flair.gtf/isoquant.gtf use different provider transcript IDs for the
same exact novel-exon structure. Named metadata includes zero, NA, grouped
counts, SQANTI3 QC/CDS annotations and an unknown future column. Provider labels
are synthetic evidence and do not determine structural novelty.

exact.fa, mismatch.fa and length.fa represent the identical FSM model with an
exact sequence, one changed base, and a truncated sequence respectively. They
are QC evidence only; genome reconstruction remains canonical.

settings-denovo.json explicitly enables ATG enumeration; all other runs default
to no de novo discovery. Separate tests construct AUG/stop/boundary/Kozak/uORF/
alternative-frame examples and invalid provider CDS coordinates. Normal-role
and shared long-read/splice origin fixtures are constructed deterministically
inside the tests to exercise the existing public entry points.
