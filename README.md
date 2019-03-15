# LncRNA orthology

This repository contains scripts related to retrieving full-transcript orthologs from a reference list of lncRNAs.

The programs were designed for Unix-based operational systems, so you can use them in any Linux or MacOS computer. For Windows systems you need either specific compilers or a Virtual Machine in addition. To obtain the repository in a Unix machine and use it locally, you can either download the repository directly from the GITHub webpage, or you can clone the repository in the terminal with the following command line:

git clone https://github.com/waltercostamb/lncRNA-ortholog-reconstruction

__***buildOrthologs.pl***__

buildOrthologs.pl retrieves full transcripts from a list of ortholog splice sites of lncRNAs. The main input for the script is a MAP_FILE from the SpliceMap Pipeline (Nitsche et al 2015); the output is in BED12 format.

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

__***retrieve-fasta***__

retrieve-fasta retrieves FASTA sequences from a BED12 file and a reference genome.

Compile the C program to an executable first:  
gcc -o retrieve-fasta retrieve-fasta.c

And then run the executable to retrieve the sequences:  
Usage: retrieve-fasta GENOME_FASTA BED_FILE

__***GFF32BED12.pl***__

GFF32BED12.pl transforms GFF3 features, or types, into BED12 lines. The parameter --type corresponds to the third column of the GFF file that you wish to transform to BED12 (e.g. CDS, mRNA, gene, etc). The script is currently set to ignore introns, so the BED12 will have only 1 block.

Usage: perl GFF32BED12.pl --gff GFF3_FILE --type TYPE

__***References***__

buildOrthologs.pl and retrieve-fasta.c are refered in the following publication:

"SSS-test: a novel test for detecting positive selection on RNA secondary structure", Maria Beatriz Walter Costa, Christian Höner zu Siederdissen, Marko Dunjic, Peter F. Stadler and Katja Nowick. BMC Bioinformatics. 2019
  https://doi.org/10.1186/s12859-019-2711-y

To build the ortholog splice sites, please consult the following publication:

"Comparison of splice sites reveals that long noncoding RNAs are evolutionarily well conserved", Anne Nitsche, Dominic Rose, Mario Fasold, Kristin Reiche and Peter F Stadler. RNA. 2015
  https://rnajournal.cshlp.org/content/21/5/801.long
  
The manuscript on the greedy approach for retrieving full transcripts of lncRNA orthologs is under preparation.
  
__***Contact***__

If you have any questions or find any problems, contact the developer: bia.walter@gmail.com

__***Author***__

Maria Beatriz Walter Costa  
https://github.com/waltercostamb
