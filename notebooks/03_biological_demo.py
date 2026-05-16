#!/usr/bin/env python
"""
Biological Demo: MoDES on stimulated vs control PBMC multiome data.

Requires Python >=3.10 and:
  wget -O /tmp/pbmc_10k_multiome.h5 \
    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/\
    pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

This demo creates pseudobulk from real PBMC data with cell-type-level
aggregation to identify regulatory events per cell type.
"""
print("Biological demo: requires Python >=3.10 and real 10x multiome H5 data.")
print("See notebooks/02_pbmc_real_data_test.py for the full pipeline script.")
