#!/bin/bash

ENV ?= ${PATH}=/software/R-3.6.1-el7-x86_64/bin/R:${PATH}
DATA_TOPDIR ?= /project2/gca/Heart_Atlas
ATAC_DIR ?= ${DATA_TOPDIR}/ATAC_seq/
RNA_DIR ?= ${DATA_TOPDIR}/Nuc_seq/

help:
	cat Makefile

process_rna:
	env ${ENV} Rscript scripts/process_RNA.R ${RNA_DIR} 

process_atac:
	env ${ENV} Rscript scripts/process_ATAC.R ${ATAC_DIR} 