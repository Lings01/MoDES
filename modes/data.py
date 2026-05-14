"""MoDEData: unified data container for bulk, pseudobulk, and single-cell input."""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


class MoDEData:
    """
    Unified internal data container for MoDES.

    Holds RNA and ATAC count matrices with sample metadata,
    regardless of input format (bulk TSV, AnnData, pseudobulk).

    Parameters
    ----------
    rna : pd.DataFrame
        RNA count matrix, shape (n_samples, n_genes).
    atac : pd.DataFrame
        ATAC peak count matrix, shape (n_samples, n_peaks).
    obs : pd.DataFrame
        Sample metadata, shape (n_samples, n_covariates).
    """

    def __init__(
        self,
        rna: pd.DataFrame,
        atac: pd.DataFrame,
        obs: pd.DataFrame,
    ):
        if not (rna.index.equals(atac.index) and rna.index.equals(obs.index)):
            raise ValueError(
                "RNA, ATAC, and obs must have matching sample indices"
            )
        self.rna = rna.astype(float)
        self.atac = atac.astype(float)
        self.obs = obs.copy()

    # -- Properties --

    @property
    def n_samples(self) -> int:
        return self.rna.shape[0]

    @property
    def n_genes(self) -> int:
        return self.rna.shape[1]

    @property
    def n_peaks(self) -> int:
        return self.atac.shape[1]

    @property
    def samples(self) -> pd.Index:
        return self.rna.index

    @property
    def gene_names(self) -> pd.Index:
        return self.rna.columns

    @property
    def peak_names(self) -> pd.Index:
        return self.atac.columns

    # -- Classmethod constructors --

    @classmethod
    def from_matrices(
        cls,
        rna_counts: Union[str, pd.DataFrame],
        atac_counts: Union[str, pd.DataFrame],
        metadata: Union[str, pd.DataFrame],
        condition_col: str,
        donor_col: Optional[str] = None,
        batch_col: Optional[str] = None,
        **kwargs,
    ) -> "MoDEData":
        """
        Build MoDEData from TSV/CSV matrices.

        Parameters
        ----------
        rna_counts : str or DataFrame
            Path to TSV or DataFrame, (samples x genes).
        atac_counts : str or DataFrame
            Path to TSV or DataFrame, (samples x peaks).
        metadata : str or DataFrame
            Path to TSV or DataFrame with sample metadata.
            Must contain at least the sample index and condition_col.
        condition_col : str
            Column name for the condition of interest.
        donor_col : str, optional
            Column name for donor/replicate.
        batch_col : str, optional
            Column name for batch.
        **kwargs
            Passed to pd.read_csv for TSV inputs.
        """
        meta_kwargs = kwargs.copy()
        meta_kwargs.setdefault("index_col", 0)
        rna = _load_matrix(rna_counts, **kwargs)
        atac = _load_matrix(atac_counts, **kwargs)
        obs = _load_matrix(metadata, **meta_kwargs)

        # Align samples
        common = rna.index.intersection(atac.index).intersection(obs.index)
        if len(common) == 0:
            raise ValueError("No common samples across RNA, ATAC, and metadata")
        if len(common) < len(rna.index) or len(common) < len(atac.index):
            warnings.warn(
                f"Aligning to {len(common)} common samples "
                f"(RNA={len(rna.index)}, ATAC={len(atac.index)}, obs={len(obs.index)})"
            )

        rna = rna.loc[common]
        atac = atac.loc[common]
        obs = obs.loc[common]

        # Ensure required columns
        for col in [condition_col, donor_col, batch_col]:
            if col and col not in obs.columns:
                raise ValueError(f"Required column '{col}' not found in metadata")

        return cls(rna=rna, atac=atac, obs=obs)

    @classmethod
    def from_anndata(
        cls,
        adata,
        condition_col: str,
        donor_col: Optional[str] = None,
        batch_col: Optional[str] = None,
        atac_layer: Optional[str] = None,
    ) -> "MoDEData":
        """
        Build MoDEData from a single AnnData with both RNA and ATAC.

        Assumes RNA is in adata.X and ATAC in adata.obsm['atac'] or a layer.

        Parameters
        ----------
        adata : AnnData
            Multiome AnnData object.
        condition_col : str
            Column in adata.obs for the condition.
        donor_col : str, optional
        batch_col : str, optional
        atac_layer : str, optional
            If provided, adata.layers[atac_layer] is the ATAC matrix.
            Otherwise checks adata.obsm['atac'].
        """
        try:
            import anndata
        except ImportError:
            raise ImportError("anndata is required for from_anndata()")

        # RNA matrix
        rna_arr = adata.X
        if hasattr(rna_arr, "toarray"):
            rna_arr = rna_arr.toarray()
        rna = pd.DataFrame(
            rna_arr,
            index=adata.obs_names,
            columns=adata.var_names,
        )

        # ATAC matrix
        if atac_layer is not None:
            atac_arr = adata.layers[atac_layer]
        elif "atac" in adata.obsm:
            atac_arr = adata.obsm["atac"]
        else:
            raise ValueError(
                "ATAC matrix not found. Provide atac_layer or store in adata.obsm['atac']"
            )

        if hasattr(atac_arr, "toarray"):
            atac_arr = atac_arr.toarray()

        # Get ATAC feature names
        atac_var_names = None
        if hasattr(adata, "uns") and "atac_var_names" in adata.uns:
            atac_var_names = adata.uns["atac_var_names"]
        else:
            atac_var_names = [f"peak_{i}" for i in range(atac_arr.shape[1])]

        atac = pd.DataFrame(
            atac_arr,
            index=adata.obs_names,
            columns=atac_var_names,
        )

        obs = adata.obs.copy()

        for col in [condition_col, donor_col, batch_col]:
            if col and col not in obs.columns:
                raise ValueError(f"Required column '{col}' not found in adata.obs")

        return cls(rna=rna, atac=atac, obs=obs)

    @classmethod
    def from_pseudobulk(
        cls,
        adata,
        groupby: List[str],
        condition_col: str,
        donor_col: str,
        min_cells_per_group: int = 10,
        aggregation: str = "sum",
        atac_layer: Optional[str] = None,
    ) -> "MoDEData":
        """
        Build MoDEData by aggregating single-cell data to pseudobulk.

        Parameters
        ----------
        adata : AnnData
            Single-cell multiome AnnData.
        groupby : list of str
            Columns in adata.obs to group by (e.g., ['donor', 'condition', 'cell_type']).
        condition_col : str
        donor_col : str
        min_cells_per_group : int
            Minimum cells per group to include.
        aggregation : str
            Aggregation function ('sum' or 'mean').
        atac_layer : str, optional
            Layer for ATAC matrix.

        Returns
        -------
        MoDEData with pseudobulk samples.
        """
        try:
            import anndata
        except ImportError:
            raise ImportError("anndata is required for from_pseudobulk()")

        obs = adata.obs.copy()

        # Validate columns
        for col in groupby:
            if col not in obs.columns:
                raise ValueError(f"Groupby column '{col}' not found in adata.obs")

        # Remove existing index name to avoid conflict
        obs.index.name = None

        # Create group identifiers
        obs["_group"] = obs[groupby].astype(str).agg("_".join, axis=1)
        group_sizes = obs.groupby("_group").size()
        valid_groups = group_sizes[group_sizes >= min_cells_per_group].index
        if len(valid_groups) == 0:
            raise ValueError(
                f"No groups with >= {min_cells_per_group} cells. "
                f"Max group size: {group_sizes.max()}"
            )

        valid_mask = obs["_group"].isin(valid_groups)
        adata_sub = adata[valid_mask.values].copy()
        adata_sub.obs["_group"] = obs.loc[adata_sub.obs_names, "_group"].astype(str).values

        # RNA matrix
        if hasattr(adata_sub.X, "toarray"):
            rna_cell = pd.DataFrame(
                adata_sub.X.toarray(),
                index=adata_sub.obs_names,
                columns=adata_sub.var_names,
            )
        else:
            rna_cell = pd.DataFrame(
                adata_sub.X,
                index=adata_sub.obs_names,
                columns=adata_sub.var_names,
            )

        rna_cell["_group"] = adata_sub.obs["_group"].values
        agg_func = rna_cell.groupby("_group")
        rna_pb = agg_func.sum() if aggregation == "sum" else agg_func.mean()
        rna_pb.drop(columns=["_group"], inplace=True, errors="ignore")

        # ATAC matrix
        if atac_layer is not None:
            atac_arr = adata_sub.layers[atac_layer]
        elif "atac" in adata_sub.obsm:
            atac_arr = adata_sub.obsm["atac"]
        else:
            raise ValueError("ATAC matrix not found")

        if hasattr(atac_arr, "toarray"):
            atac_arr = atac_arr.toarray()

        atac_var_names = None
        if hasattr(adata, "uns") and "atac_var_names" in adata.uns:
            atac_var_names = adata.uns["atac_var_names"]
        else:
            atac_var_names = [f"peak_{i}" for i in range(atac_arr.shape[1])]

        atac_cell = pd.DataFrame(
            atac_arr,
            index=adata_sub.obs_names,
            columns=atac_var_names,
        )
        atac_cell["_group"] = adata_sub.obs["_group"].values
        agg_func_atac = atac_cell.groupby("_group")
        atac_pb = agg_func_atac.sum() if aggregation == "sum" else agg_func_atac.mean()
        atac_pb.drop(columns=["_group"], inplace=True, errors="ignore")

        # Metadata: aggregate mode / first for categorical, mean for numeric
        obs_sub = adata_sub.obs.copy()
        obs_agg = obs_sub.groupby("_group").agg({
            c: ("first" if obs_sub[c].dtype == object or obs_sub[c].dtype.name == "category" else "mean")
            for c in obs_sub.columns
            if c != "_group"
        })

        # Build context label
        context_parts = [c for c in groupby if c not in (donor_col,)]
        obs_agg["context"] = obs_agg[context_parts].astype(str).agg("_".join, axis=1)

        # Build sample names
        obs_agg.index = obs_agg.index.astype(str)

        return cls(rna=rna_pb, atac=atac_pb, obs=obs_agg)

    # -- Methods --

    def validate(self, condition_col: str = None) -> List[str]:
        """Check data integrity. Returns list of issues found."""
        issues = []
        # Sample count
        if self.n_samples < 3:
            issues.append("Fewer than 3 samples")
        if self.n_samples < 8:
            issues.append(
                f"Only {self.n_samples} samples; >=8 recommended for reliable "
                "variance estimation"
            )
        # Index alignment
        if not self.rna.index.equals(self.atac.index):
            issues.append("RNA and ATAC sample indices differ")
        if not self.rna.index.equals(self.obs.index):
            issues.append("RNA and metadata sample indices differ")
        # Duplicate IDs
        if len(self.rna.columns) != len(set(self.rna.columns)):
            issues.append("Duplicate gene names in RNA matrix")
        if len(self.atac.columns) != len(set(self.atac.columns)):
            issues.append("Duplicate peak IDs in ATAC matrix")
        if len(self.rna.index) != len(set(self.rna.index)):
            issues.append("Duplicate sample IDs in RNA matrix")
        # Missing values
        if self.rna.isnull().any().any():
            issues.append("RNA matrix contains NaN values")
        if self.atac.isnull().any().any():
            issues.append("ATAC matrix contains NaN values")
        if self.obs.isnull().any().any():
            issues.append("Metadata contains NaN values")
        # Negative / non-integer
        if (self.rna.values < 0).any():
            issues.append("RNA matrix contains negative values")
        if (self.atac.values < 0).any():
            issues.append("ATAC matrix contains negative values")
        # Library size
        rna_sum = self.rna.sum(axis=1)
        atac_sum = self.atac.sum(axis=1)
        if (rna_sum <= 0).any():
            bad = list(self.rna.index[rna_sum <= 0])
            issues.append(f"Samples with zero RNA library size: {bad}")
        if (atac_sum <= 0).any():
            bad = list(self.atac.index[atac_sum <= 0])
            issues.append(f"Samples with zero ATAC library size: {bad}")
        # Infinite values
        if np.isinf(self.rna.values).any():
            issues.append("RNA matrix contains infinite values")
        if np.isinf(self.atac.values).any():
            issues.append("ATAC matrix contains infinite values")
        # Condition column
        if condition_col is not None:
            if condition_col not in self.obs.columns:
                issues.append(f"Condition column '{condition_col}' not in metadata")
            else:
                cond = self.obs[condition_col]
                n_cats = len(set(cond.dropna()))
                if n_cats != 2:
                    issues.append(
                        f"Condition '{condition_col}' has {n_cats} unique values; "
                        "binary (2) is required"
                    )
        return issues

    def filter_features(
        self,
        min_counts: int = 1,
        min_rna_samples: int = 3,
        min_atac_samples: int = 3,
    ) -> "MoDEData":
        """Filter low-coverage genes and peaks."""
        rna_mask = (self.rna > 0).sum(axis=0) >= min_rna_samples
        rna_mask &= self.rna.sum(axis=0) >= min_counts

        atac_mask = (self.atac > 0).sum(axis=0) >= min_atac_samples
        atac_mask &= self.atac.sum(axis=0) >= min_counts

        return MoDEData(
            rna=self.rna.loc[:, rna_mask],
            atac=self.atac.loc[:, atac_mask],
            obs=self.obs.copy(),
        )

    def normalize_library_size(self, method: str = "median_ratio") -> "MoDEData":
        """
        Compute library-size-normalized counts (returned as new matrices).
        Does NOT modify in place.

        Parameters
        ----------
        method : str
            'median_ratio' (DESeq2-style) or 'cpm' (counts per million).

        Returns
        -------
        MoDEData with normalized count matrices.
        """
        if method == "median_ratio":
            rna_norm = _median_ratio_normalize(self.rna.values)
            atac_norm = _median_ratio_normalize(self.atac.values)
        elif method == "cpm":
            rna_sum = self.rna.sum(axis=1).values
            atac_sum = self.atac.sum(axis=1).values
            if np.any(rna_sum <= 0) or np.any(atac_sum <= 0):
                raise ValueError(
                    "Cannot CPM-normalize samples with zero library size."
                )
            rna_sf = rna_sum / 1e6
            atac_sf = atac_sum / 1e6
            rna_norm = self.rna.values / rna_sf[:, None]
            atac_norm = self.atac.values / atac_sf[:, None]
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        return MoDEData(
            rna=pd.DataFrame(rna_norm, index=self.rna.index, columns=self.rna.columns),
            atac=pd.DataFrame(atac_norm, index=self.atac.index, columns=self.atac.columns),
            obs=self.obs.copy(),
        )

    def get_library_sizes(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return log library sizes for RNA and ATAC (used as offsets)."""
        rna_sum = self.rna.sum(axis=1).values
        atac_sum = self.atac.sum(axis=1).values
        if np.any(rna_sum <= 0):
            bad = list(self.rna.index[rna_sum <= 0])
            raise ValueError(f"Samples with zero RNA library size: {bad}")
        if np.any(atac_sum <= 0):
            bad = list(self.atac.index[atac_sum <= 0])
            raise ValueError(f"Samples with zero ATAC library size: {bad}")
        return np.log(rna_sum), np.log(atac_sum)


# -- Internal helpers --

def _load_matrix(path_or_df: Union[str, pd.DataFrame], **kwargs) -> pd.DataFrame:
    """Load a matrix from file or return DataFrame as-is."""
    if isinstance(path_or_df, pd.DataFrame):
        return path_or_df.copy()
    # Auto-detect separator
    sep = kwargs.pop("sep", None)
    if sep is None:
        sep = "\t" if str(path_or_df).endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path_or_df, sep=sep, **kwargs)


def _median_ratio_normalize(counts: np.ndarray) -> np.ndarray:
    """
    Median-of-ratios normalization (DESeq2-style).

    Returns normalized counts (not log-transformed).
    """
    counts = np.asarray(counts, dtype=float)
    # Geometric mean per feature (excluding zeros)
    log_counts = np.log(counts.clip(min=1))
    geo_means = np.exp(log_counts.mean(axis=0))
    # Size factor per sample: median of ratios to geometric mean
    ratios = counts / geo_means[None, :]
    size_factors = np.median(ratios, axis=1)
    size_factors = np.clip(size_factors, 0.1, 10.0)
    return counts / size_factors[:, None]
