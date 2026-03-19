# Simple_BGE_Sample_QC

A Hail-based pipeline for quality control and genetic ancestry prediction.

## Features
- Import and filter VCF files with comprehensive QC metrics
- Sex imputation from X chromosome heterozygosity
- Region-specific QC (exome vs. non-exome)
- Genetic ancestry inference using 1000 Genomes reference populations
- PCA-based projection for ancestry assignment
- Relatedness filtering for unbiased reference panel

## Pipeline Steps
1. VCF import and MatrixTable conversion
2. Hard variant filtering (bi-allelic SNPs, MAF, call rate)
3. Genetic sex imputation (F-statistic method)
4. Sample QC metrics (Ti/Tv, het/hom ratios, variant counts)
5. 1000 Genomes reference processing (QC, LD pruning, relatedness removal)
6. Principal Component Analysis with projection
7. Random Forest ancestry classification

## Requirements
- Hail 0.2+
- PySpark
- gnomAD methods library
- 1000 Genomes Phase 3 VCF files
- BED file for exome capture regions

