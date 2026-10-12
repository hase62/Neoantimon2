# Independent T008 gold reference

All sequence below was hand-authored before using the reconstruction code.
Coordinates are 1-based inclusive in the GTF and genome FASTA.

| Transcript | CDS intervals in transcript order | Reference CDS | Protein |
| --- | --- | --- | --- |
| TA.1, plus | chr1:1–8, 21–30 | `ATGGAACCTAAACCCGGG` | `MEPKPG` |
| TB.1, plus | chr1:1–5, 21–30 | `ATGGATAAACCCGGG` | `MDKPG` |
| TC.1, minus | chr2:21–28, 4–13 | `ATGAAACCTGGGAAACCC` | `MKPGKP` |

TA.1 exon 1 contributes `CC` and exon 2 contributes `T` to `CCT` (Pro).
TC.1 high-coordinate exon contributes `CC` and low-coordinate exon contributes
`T` to its split `CCT` (Pro). Both following CDS segments have GTF phase 1;
their first base remains part of the biological CDS.

The complete transcript sequence adds the terminal `TAA` stop after each CDS.
TD.1 is on chrS with a `TGA` codon at position 2, a Sec annotation and a
protein catalog `U` at that position; it must be deferred. TM.1 is mitochondrial and
must also be deferred. TX.1 is a catalog-only Sec record, with no primary-GTF
coordinate authority.

Selected hand-authored variant oracles:

| Change | Tumor CDS | Protein |
| --- | --- | --- |
| chr1:4 G>A | `ATGAAACCTAAACCCGGG` | `MKPKPG` |
| chr2:25 T>C | `ATGGAACCTGGGAAACCC` | `MEPGKP` |
| chr1:4–5 GA>AC | `ATGACACCTAAACCCGGG` | `MTPKPG` |
| chr1:6 A>ACCC | `ATGGAACCCCCTAAACCCGGG` | `MEPPKPG` |
| chr1:3 GGAA>G | `ATGCCTAAACCCGGG` | `MPKPG` |
| chr2:25 T>TCCC | `ATGGGGAAACCTGGGAAACCC` | `MGKPGKP` |
| chr2:22 GTTT>G | `ATGCCTGGGAAACCC` | `MPGKP` |

The chr2 alleles are forward-genomic alleles. No expected value in this file
or the focused test is generated from a T008 function.

## T009 hand-authored translation oracles

The added `TF.1` and `TR.1` transcripts have the same 42-base spliced sequence:

`ATG AAA CCC GGG AAA | CCC GGG AAA CCC TAA CCC GGG TAA GGG`

The bar is the exon junction after CDS base 15. The 10-base genomic intron is
absent from the spliced sequence. The canonical protein is `MKPGKPGKP`, with
stop at amino-acid position 10. `TR.1` is on the minus strand and its source
alleles are forward-genomic. Changing the annotated stop from `TAA` to `CAA`
on either strand yields `MKPGKPGKPQPG`; the new stop is position 13. The
extension is `QPG`, from the spliced 3' transcript sequence. `TO.1` has no
downstream stop in its known extension, giving the bounded prefix
`MKPGKPGKPQPGK` after the same stop loss.

For the central frame-restoration oracle, `TF.1` starts with CDS
`ATGAAACCCGGGAAACCCGGGAAACCC` and protein `MKPGKPGKP`. Inserting `C` after
coding base 6 and deleting coding base 17 on the same side gives CDS
`ATGAAACCCCGGGAAACCGGGAAACCC`, codons
`ATG AAA CCC CGG GAA ACC GGG AAA CCC`, and protein `MKPRETGKP`. Frame is
shifted after the insertion and restored after the deletion. The first fully
restored codon begins at tumor coding base 19, amino acid 7. Amino acids 4–6
change from `GKP` to `RET`; amino acid 3 remains `P` despite changed codon
origin. Deleting coding base 6 then inserting `C` after coding base 17 gives
`ATGAACCCGGGAAACCCCGGGAAACCC` and `MNPGNPGKP`.

`TU.1` has CDS `ATG AAA CCC GGG CTA AAC CCG GGT` and protein `MKPGLNPG`.
Deleting one `A` in codon 2 gives `ATG AAC CCG GGC TAA` before termination,
protein `MNPG`, and stop at position 5. `TV.1` contains a second, downstream
`ATG` after its annotated start; destroying the first start leaves canonical
tumor initiation unresolved.

`TS.1` has a GENCODE Sec annotation at amino acid 8: its reference CDS codons
are `ATG AAA CCC GGG AAA CCC GGG TGA CCC`, protein `MKPGKPGUP`. The exact
annotated eighth `TGA` is `U`; an unannotated nuclear `TGA` remains a stop.
`TM2.1` has mitochondrial `ATA TGA AGA`, producing `MW` followed by stop.
`TM3.1` has `ATA TGA T`; its terminal `T` is a partial stop accepted only by
the validated mitochondrial catalog rule. No nucleotide is appended to the
stored transcript sequence.

`TW.1` shares `TF.1`'s annotated protein and stop at amino acid 10, but has
the hand-authored downstream transcript tail `CCGGGTAAAGGG`. Inserting `C`
after coding base 6 shifts the frame without restoring it. The original stop
is lost in the shifted frame, and the tail's original bases 36–38 form a new
`TAA` stop at mutant amino acid 13. This separates the frame and termination
effects in a real translated fixture.

These literal expected CDS and protein sequences are independent of T009
translation helpers. The T009 tests also use them for compound codon changes,
frame restoration, stop before restoration, and protein difference intervals.
