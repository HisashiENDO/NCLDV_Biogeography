#!/bin/bash
#===================================================================================================
# Copyright (c) 2020 Hisashi Endo and Yanze Li
#
# Pplacer extracts polB sequences mapped onto the NCLDV references from the environmental sequence.
# Then, parse_jplace_all.py script generate a ID list of these sequences.
# Finally, ectract_treget_line.py script generate abundance table(s) of these genes.
#===================================================================================================

#sepicify reference files
REF_ALN="pplacer_ref211.msa"
RAXML_INFO="RAxML_info.pplacer.txt"
RAXML_RESULT="RAxML_result.pplacer.nwk"

#specify query sequences
QUERY_SEQ="$1"

#Align query sequences with reference
mafft --thread 16 --6merpair --addfragments ${QUERY_SEQ} ${REF_ALN} > ${QUERY_SEQ}.combo.fasta

#Pplacer
pplacer -j 16 --verbosity 0 --keep-at-most 1 -o ${QUERY_SEQ}.combo.jplace -t ${RAXML_RESULT} -s ${RAXML_INFO} ${QUERY_SEQ}.combo.fasta

#Decipher Pplacer output to sequence IDs
python3 parse_jplace.py ${QUERY_SEQ}.combo.jplace > ${QUERY_SEQ}.tit

#Extract target lines from an abundance profiles
python3 ectract_treget_line.py ${QUERY_SEQ}.tit

# Excute with: sh run_pplacer.sh [query]
