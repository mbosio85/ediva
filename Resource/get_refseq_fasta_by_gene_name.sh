#!/bin/bash

BIN=/users/GD/tools/ediva/Resource
GENEFILE=file_containing_genes_per_line
BPS=100 ## amount in bps to extend on both sides of the exons
FASTAOUT=output.fasta ## output fasta file name
TLOC=temp_location ## will be cleared after the script execution

perl $BIN/get_genes_by_db.pl $GENEFILE
perl $BIN/get_ref_fasta_by_gene_name.pl $GENEFILE $BPS $FASTAOUT $TLOC
