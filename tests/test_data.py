"""Tests for MoDEData container."""

import pytest
import numpy as np
import pandas as pd

from modes.data import MoDEData


class TestMoDEData:
    def test_init(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        assert data.n_samples == 20
        assert data.n_genes == 30
        assert data.n_peaks == 60
        assert data.rna.index.equals(data.atac.index)
        assert data.rna.index.equals(data.obs.index)

    def test_validate(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        issues = data.validate()
        # Should be clean for synthetic data
        assert len(issues) == 0

    def test_properties(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        assert len(data.samples) == data.n_samples
        assert len(data.gene_names) == data.n_genes
        assert len(data.peak_names) == data.n_peaks

    def test_filter_features(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        filtered = data.filter_features(min_counts=5, min_rna_samples=5)
        assert filtered.n_genes <= data.n_genes
        assert filtered.n_peaks <= data.n_peaks
        assert filtered.n_samples == data.n_samples

    def test_normalize_cpm(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        norm = data.normalize_library_size(method="cpm")
        assert norm.n_samples == data.n_samples
        assert norm.n_genes == data.n_genes

    def test_normalize_median_ratio(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        norm = data.normalize_library_size(method="median_ratio")
        assert norm.n_samples == data.n_samples

    def test_get_library_sizes(self, synthetic_bulk_data):
        data, _, _ = synthetic_bulk_data
        rna_ls, atac_ls = data.get_library_sizes()
        assert len(rna_ls) == data.n_samples
        assert len(atac_ls) == data.n_samples


class TestMoDEDataConstructors:
    def test_from_matrices_tsv(self, tmp_path):
        """Test loading from TSV files."""
        rna_path = tmp_path / "rna.tsv"
        atac_path = tmp_path / "atac.tsv"
        meta_path = tmp_path / "meta.tsv"

        # Write test files
        samples = [f"sample_{i}" for i in range(6)]
        rna_df = pd.DataFrame(
            np.random.poisson(10, (6, 4)),
            index=samples,
            columns=["gene_0", "gene_1", "gene_2", "gene_3"],
        )
        atac_df = pd.DataFrame(
            np.random.poisson(5, (6, 6)),
            index=samples,
            columns=["peak_0", "peak_1", "peak_2", "peak_3", "peak_4", "peak_5"],
        )
        meta = pd.DataFrame({
            "condition": ["ctrl", "ctrl", "ctrl", "trt", "trt", "trt"],
            "batch": ["A", "A", "B", "A", "B", "B"],
        }, index=samples)
        meta.index.name = "sample"

        rna_df.to_csv(rna_path, sep="\t")
        atac_df.to_csv(atac_path, sep="\t")
        meta.to_csv(meta_path, sep="\t")

        data = MoDEData.from_matrices(
            str(rna_path), str(atac_path), str(meta_path),
            condition_col="condition",
            index_col=0,
        )
        assert data.n_samples == 6
        assert data.n_genes == 4
        assert data.n_peaks == 6

    def test_from_matrices_mismatch_warns(self, tmp_path):
        """Test alignment when samples don't match."""
        rna_path = tmp_path / "rna.tsv"
        atac_path = tmp_path / "atac.tsv"
        meta_path = tmp_path / "meta.tsv"

        samples_rna = ["s1", "s2", "s3", "s4"]
        samples_atac = ["s2", "s3", "s4"]
        samples_meta = ["s1", "s2", "s3", "s4"]

        rna_df = pd.DataFrame(np.ones((4, 2)), index=samples_rna, columns=["g1", "g2"])
        atac_df = pd.DataFrame(np.ones((3, 2)), index=samples_atac, columns=["p1", "p2"])
        meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]}, index=samples_meta)
        meta.index.name = "sample"

        rna_df.to_csv(rna_path, sep="\t")
        atac_df.to_csv(atac_path, sep="\t")
        meta.to_csv(meta_path, sep="\t")

        data = MoDEData.from_matrices(
            str(rna_path), str(atac_path), str(meta_path),
            condition_col="condition",
            index_col=0,
        )
        # Should align to intersection: s2, s3, s4
        assert data.n_samples == 3
