# numts
R script to detect numts in the COI Leray fragment.

It checks sequences against all genetic codes (or all metazoan codes if it is a metazoan) to detect stop codons. 
For metazoans and seqs of 313 bp it also checks the five aminoacids conserved according to Pentinsaari et al 2016. 
We assume that the codon starts in position 2 in the Leray fragment. Most sequences should be 313 bp. Other lengths are acceptable if codon coherence is kept (i.e., length is 313+-3*n)

