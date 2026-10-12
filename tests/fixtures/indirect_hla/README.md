# Deterministic indirect-HLA mechanics fixtures

All sequences, allele mappings and provider outputs are synthetic. Allele
spellings exercise the frozen HLA syntax; these are not real HLA proteins or
actual NetMHCIIpan predictions. No executable, model weights or patient data
is included.

The five donor source proteins share a 20-residue sequence: six terminal-inclusive
15-mers per origin, six molecular ligands, thirty distinct occurrences. The five
self proteins share one 15-mer whose primary forward nine-residue core matches
one source core on the explicitly typed DR presenter. All affinities are 500 nM.
The independent oracle is five included core/presenter units and one excluded
self unit, with thirty occurrences and six source pMHCs. DRA is typed explicitly
for presentation and excluded from the five-locus source/self generation space.

Focused tests perturb these records for repetition, shared cores, alternative
registers, threshold boundaries, missingness, ambiguity and provenance/order.
