// A program to extract portions of molecules in a multi-fasta file, evaluating
// the reverse-complement if the gene are in the minus strand.

// This program will create an index in the first execution on a fasta file.
// Consecutive executions on the same file will read the index.  For a
// fasta file x the index will be named x.rfindx.

// Usage: retrieve-fasta fasta-file bed-file

// Entries in the bed file must be complete, and the number of exons may be 0.

// Guilherme P. Telles and Maria Beatriz Walter Costa, 2014-11-05.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>


void badbed() {
  printf("Bad bed file.\n");
  exit(1);
}


void die() {
  printf("Error: %s.\n",errno ? strerror(errno) : "errno not set");
  exit(errno?errno:1);
}


int is_blank(char* s) {
  if (!s)
    return 0;

  int i = 0;
  while (s[i]) {
    if (s[i] != ' ' && s[i] != '\n' && s[i] != '\t')
      return 0;
    i++;
  }
  return 1;
}


// The fasta file will be indexed.  The index is an array of struct entry:
typedef struct {
  char* tag;
  long header;
  long sequence;
} entry;


// Comparison functions:
int sort_entry(const void *a, const void *b) {
  return strcmp(((entry*)a)->tag, ((entry*)b)->tag);
}

int bs_entry(const void *k, const void *a) {
  return strcmp((char*) k, ((entry*)a)->tag);
}


// The node of a list to hold data on sequence offsets in the genome fasta
// during index construction.  That will avoid reading the fasta file twice.
typedef struct node {
  char* tag;
  long header;
  long sequence;
  struct node* next;
} node;



int main(int argc, char** argv) {

  int i,j,l,k;

  // The index file:
  char* indexf = malloc((strlen(argv[1])+20)*sizeof(char));
  if (!indexf) die();
  sprintf(indexf,"%s.rfindx",argv[1]);

  // The index array and its lenght:
  entry* index = 0;
  int indexl = 0;

  int build_index = 1;
  struct stat fastast, indexst;

  if (stat(argv[1], &fastast) != 0) die();

  if (access(indexf,R_OK) == 0) {
    if (stat(indexf, &indexst) != 0) die();

    if (fastast.st_mtime < indexst.st_mtime)
      build_index = 0;
  }

  if (!build_index) {
    // Loads index into memory:
    FILE* indexh = fopen(indexf,"rb");
    if (!indexh) die();

    fread(&indexl,sizeof(int),1,indexh);

    index = malloc(indexl*sizeof(entry));
    if (!index) die();

    for (i=0; i<indexl; i++) {
      fread(&l,sizeof(int),1,indexh);
      index[i].tag = malloc((l+1)*sizeof(char));
      if (!index[i].tag) die();

      fread(index[i].tag,sizeof(char),l,indexh);
      index[i].tag[l] = 0;
      fread(&(index[i].header),sizeof(long),1,indexh);
      fread(&(index[i].sequence),sizeof(long),1,indexh);
    }

    fclose(indexh);
  }
  else {
    // Builds an index and saves to file:
    node* I = 0;

    // Opens fasta file:
    FILE* fastah = fopen(argv[1],"rt");
    if (!fastah) die();

    // Reads fasta char-by-char and adds sequence offsets to list I:
    char c, pc = '\n';
    while (pc != EOF && (c = fgetc(fastah)) != EOF) {

      if (c == '>' && pc == '\n') {
        c = fgetc(fastah);
        node* seq = malloc(sizeof(node));
        if (!seq) die();

        seq->next = I;
        I = seq;
        seq->header = ftell(fastah)-1;

        char name[64];
        i = 0;
        while (c != ' ' && c != '\t' && c != '\n' && c != '|') {
          name[i] = c;
          c = fgetc(fastah);
          i++;
        }
        name[i] = 0;

        seq->tag = strdup(name);

        while (c != '\n')
          c = fgetc(fastah);

        seq->sequence = ftell(fastah);
        indexl++;
      }

      pc = c;
    }

    fclose(fastah);

    // Builds an array from the list, and releases the list:
    index = malloc(indexl*sizeof(entry));
    if (!index) die();

    node* p = I;
    i = 0;
    while (p) {
      index[i].tag = p->tag;
      index[i].header = p->header;
      index[i].sequence = p->sequence;

      node* d = p;
      p = p->next;
      free(d);

      i++;
    }

    qsort(index,indexl,sizeof(entry),sort_entry);

    // Writes the index file:
    FILE* indexh = fopen(indexf,"wb");
    if (!indexh) die();

    fwrite(&indexl,sizeof(int),1,indexh);

    for (i=0; i<indexl; i++) {
      l = strlen(index[i].tag);
      fwrite(&l,sizeof(int),1,indexh);
      fwrite(index[i].tag,sizeof(char),l,indexh);
      long off = index[i].header;
      fwrite(&off,sizeof(long),1,indexh);
      off = index[i].sequence;
      fwrite(&off,sizeof(long),1,indexh);
    }

    fclose(indexh);
  }


  // Opens fasta:
  FILE* fastah = fopen(argv[1],"rt");
  if (!fastah) die();

  // Opens bed:
  FILE* bed = fopen(argv[2],"rt");
  if (!bed) die();

  // The retrieved molecule:
  char* rna = 0;
  int rnas = 0;
  int rnal = 0;

  // Process bed lines:
  char *line = 0;
  size_t lines = 0;
  size_t linel;

  while ((linel = getline(&line, &lines, bed)) != -1) {

    if (is_blank(line)) continue;

    char chr_name[64];
    char gene_name[64];
    int gene_start, gene_end;
    char strand;

    if (sscanf(line,"%[^\t]\t%d\t%d\t%[^\t]\t%*d\t%c\t%n",
               chr_name,&gene_start,&gene_end,gene_name,&strand,&k) != 5)
      badbed();
    j = k;

    for (i=0; i<3; i++) {
      sscanf(line+j,"%*s%n",&k);
      j += k;
    }

    int exons_number;
    if (sscanf(line+j,"%d %n",&exons_number,&k) != 1)
      badbed();
    j += k;

    int exon_start[exons_number+1];
    int exon_end[exons_number+1];
    int exon_len[exons_number+1];
    int length;

    if (exons_number > 0) {
      for (i=0; i<exons_number; i++) {
        if (sscanf(line+j,"%d%*[,\t]%n",&exon_len[i],&k) != 1)
          badbed();
        j += k;
      }

      for (i=0; i<exons_number; i++) {
        if (sscanf(line+j,"%d%*[,\t\n]%n",&exon_start[i],&k) != 1)
          badbed();
        j += k;
	exon_start[i] += gene_start;
      }

      // Evaluates exons' end:
      for (i=0; i<exons_number; i++)
        exon_end[i] = exon_start[i]+exon_len[i];

      // Evaluates the length of the resulting molecule:
      length = gene_end - gene_start;
      length -= exon_start[0] - gene_start;
      for (i=0; i<exons_number-1; i++)
        length -= exon_start[i+1] - (exon_start[i] + exon_len[i]);
      length -= gene_end - (exon_start[exons_number-1] + exon_len[exons_number-1]);
    }
    else {
      exon_start[0] = gene_start;
      exon_end[0] = gene_end;
      length = gene_end - gene_start + 1;
      exons_number = 1;
    }

    // The rna product and its length:
    if (length+3 > rnas) {
      rnas = length+3;
      rna = malloc(rnas*sizeof(char));
      if (!rna) die();
    }
    rnal = 0;

    // Retrieves chromossome position in the index:
    entry* match = bsearch(chr_name,index,indexl,sizeof(entry),bs_entry);
    if (!match) {
      printf("The sequence %s is not in the genome.\n",chr_name);
    }
    else {
      //for (i=0; i<exons_number; i++)
      //  printf("-%d-%d-%d-\n",exon_start[i],exon_end[i],exon_len[i]);

      fseek(fastah,match->sequence,SEEK_SET);

      // Reads the chromossome and copies characters in exons to rna:
      char c =  fgetc(fastah);
      char pc = '\n';

      int p = 0;  // the current letter in the chromossome.
      int e = 0;  // the current exon.

      while (c != EOF) {
        if (c == '>' && pc == '\n')
          c = EOF;
        else if (c != '\n' && c != ' ' && c != '\t') {
          if (p >= exon_start[e] && p < exon_end[e])
            rna[rnal++] = c;
          else if (p == exon_end[e]) {
            e++;
            if (e < exons_number && p >= exon_start[e] && p < exon_end[e])
              rna[rnal++] = c;
          }

          if (e == exons_number)
            c = EOF;
          else {
            pc = c;
            c = fgetc(fastah);
          }
          p++;
        }
        else {
          pc = c;
          c = fgetc(fastah);
        }
      }

      if (rnal > 0) {
        rna[rnal] = 0;
        printf(">%s:%d-%d:%s:%c\n",chr_name,gene_start,gene_end,gene_name,strand);

        if (strand == '-') {
          char tr[256];

          for (i=0; i<256; i++)
            tr[i] = i;

          tr['A'] = 'T';  tr['a'] = 't';
          tr['C'] = 'G';  tr['c'] = 'g';
          tr['G'] = 'C';  tr['g'] = 'c';
          tr['T'] = 'A';  tr['t'] = 'a';
          tr['R'] = 'Y';  tr['r'] = 'y';
          tr['Y'] = 'R';  tr['y'] = 'r';
          tr['K'] = 'M';  tr['k'] = 'm';
          tr['M'] = 'K';  tr['m'] = 'k';
          tr['B'] = 'V';  tr['b'] = 'v';
          tr['V'] = 'B';  tr['v'] = 'b';
          tr['D'] = 'H';  tr['d'] = 'h';
          tr['H'] = 'D';  tr['h'] = 'd';

          for (i=rnal-1; i>=0; i--)
            printf("%c",tr[(int)rna[i]]);

          printf("\n");
        }
        else
          printf("%s\n",rna);
      }
      else
        printf("The region was not found.\n");
    }
  }
  if (ferror(bed))
    die();


  fclose(bed);
  fclose(fastah);

  free(indexf);
  free(line);
  free(rna);

  for (i=0; i<indexl; i++)
    free(index[i].tag);
  free(index);

  return 0;
}
