"""Tests for MoDEData container."""

import numpy as np
import pandas as pd
import pytest

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


def test_from_pseudobulk_minimal_anndata():
    """Test from_pseudobulk with a minimal AnnData object."""
    import anndata
    import numpy as np
    import pandas as pd

    n_cells = 30
    rng = np.random.default_rng(42)

    adata = anndata.AnnData(
        X=rng.poisson(10, (n_cells, 5)).astype(float),
        obs=pd.DataFrame({
            "donor": [f"donor_{i % 3}" for i in range(n_cells)],
            "condition": ["ctrl"] * 15 + ["trt"] * 15,
            "cell_type": ["T_cell"] * 10 + ["B_cell"] * 10 + ["NK"] * 10,
        }),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4", "G5"]),
    )
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    # Add ATAC in obsm
    adata.obsm["atac"] = rng.poisson(5, (n_cells, 4)).astype(float)
    adata.uns["atac_var_names"] = ["peak_1", "peak_2", "peak_3", "peak_4"]

    from modes.data import MoDEData
    data = MoDEData.from_pseudobulk(
        adata,
        groupby=["donor", "condition", "cell_type"],
        condition_col="condition",
        donor_col="donor",
        min_cells_per_group=2,
    )
    assert data.n_samples > 0
    assert "context" in data.obs.columns


def test_zero_library_size_rejected(synthetic_bulk_data):
    data, _, _ = synthetic_bulk_data
    data.rna.iloc[0, :] = 0
    with pytest.raises(ValueError, match="zero RNA library size"):
        data.get_library_sizes()


def test_cpm_normalization_rejects_zero_library(synthetic_bulk_data):
    data, _, _ = synthetic_bulk_data
    data.rna.iloc[0, :] = 0
    with pytest.raises(ValueError):
        data.normalize_library_size(method="cpm")


def test_validate_finds_duplicate_genes(synthetic_bulk_data):
    data, _, _ = synthetic_bulk_data
    cols = list(data.rna.columns)
    cols[1] = cols[0]
    data.rna.columns = cols
    issues = data.validate()
    assert any("Duplicate gene" in i for i in issues)


def test_validate_finds_condition_binary(synthetic_bulk_data_small):
    data, _, _ = synthetic_bulk_data_small
    n = data.n_samples
    data.obs["condition"] = ["A", "B", "C"] * (n // 3) + ["A"] * (n % 3)
    issues = data.validate(condition_col="condition")
    assert any("binary" in i.lower() for i in issues)


def test_from_mudata_toy():
    """Test from_mudata with a toy MuData-like double AnnData setup."""
    import anndata
    import numpy as np
    import pandas as pd

    class ToyMuData:
        def __init__(self, adata_rna, adata_atac):
            self.mod = {"rna": adata_rna, "atac": adata_atac}

    rng = np.random.default_rng(42)
    n_cells = 60

    adata_rna = anndata.AnnData(
        X=rng.poisson(10, (n_cells, 5)).astype(float),
        obs=pd.DataFrame({
            "donor": [f"d{i % 3}" for i in range(n_cells)],
            "condition": ["ctrl"] * 30 + ["trt"] * 30,
            "cell_type": ["T"] * 20 + ["B"] * 20 + ["NK"] * 20,
        }),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4", "G5"]),
    )
    adata_atac = anndata.AnnData(
        X=rng.poisson(5, (n_cells, 4)).astype(float),
        obs=adata_rna.obs,
        var=pd.DataFrame(index=["p1", "p2", "p3", "p4"]),
    )

    mdata = ToyMuData(adata_rna, adata_atac)

    from modes.data import MoDEData
    data = MoDEData.from_mudata(
        mdata,
        rna_mod="rna", atac_mod="atac",
        groupby=["donor", "condition", "cell_type"],
        condition_col="condition", donor_col="donor",
        min_cells_per_group=2,
    )
    assert data.n_samples > 0
    assert data.n_genes == 5
    assert data.n_peaks == 4
    assert "context" in data.obs.columns
