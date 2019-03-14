# lncRNA-ortholog-reconstruction

This repository contains the retrieve-fasta program to retrieve FASTA sequences from a BED12 file; a script to transform GFF3 features into BED12 lines and a script for retrieving full transcript orthologs from a list of ortholog splice sites of lncRNAs.

#buildOrthologs.pl

Usage: buildOrthologs.pl --input MAP_FILE --db REF_BED --map_threshold N --reference SPECIES --blast_address FILE --blast_threshold I --blast_thread T --db_reg_exp EXP --lower_seq_thresh L --upper_seq_thresh U

mandatory parameters:

--input: site MAP_FILE from the SpliceMap Pipeline
--db: REF_BED, or BED12 reference lncRNA database
--map_threshold: splice map threshold for excluding unlikely real splice sites, N can be: undef, -3, 0, 3, etc. (recommended stringent threshold: 3)
--reference: alias of the reference species, example: hg38, ponAbe2, rheMac3, mm10, etc. as in the input MAP_FILE
--blast_address: a file containing two rows divided by a tab: first row with species alias (hg38, mm10, etc.) and second row with full path to index files for BLASTN
--blast_threshold: e-value I threshold for filtering BLASTN hits
--blast_thread: T number of threads to run BLASTN
--db_reg_exp: regular expression EXP used to select lines from the --input MAP_FILE based on the --db REF_BED (e.g. ENST for human Gencode transcripts, ENSMUS for mouse, etc.)
--lower_seq_thresh: lower sequence length threshold L of input sequence to BLAST (anything shorter than L is not submitted to BLAST and marked as invalid reconstruction)
--upper_seq_thresh: upper sequence length threshold U of input sequence to BLAST (anything longer than U is not submitted to BLAST and marked as invalid reconstruction)

#retrieve-fasta

Usage: retrieve-fasta GENOME_FASTA BED_FILE

#GFF32BED12.pl

Usage: perl GFF32BED12.pl --gff GFF3_FILE --type TYPE

#References

buildOrthologs.pl and retrieve-fasta.c are refered in the following publication:

"SSS-test: a novel test for detecting positive selection on RNA secondary structure", Maria Beatriz Walter Costa, Christian Höner zu Siederdissen, Marko Dunjic, Peter F. Stadler and Katja Nowick. BMC Bioinformatics. 2019
  https://doi.org/10.1186/s12859-019-2711-y

To build the ortholog splice sites, please consult the following publication:

"Comparison of splice sites reveals that long noncoding RNAs are evolutionarily well conserved", Anne Nitsche, Dominic Rose, Mario Fasold, Kristin Reiche and Peter F Stadler. RNA. 2015
  https://rnajournal.cshlp.org/content/21/5/801.long
  
