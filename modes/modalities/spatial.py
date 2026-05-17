"""Native spatial graph engine for MoDES-Spatial (v2.0)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from dataclasses import dataclass


@dataclass
class SpatialEvidence:
    """Per-event spatial evidence metrics."""
    moran_i: float = 0.0
    neighbor_effect: float = 0.0
    edge_score: float = 0.0
    region_autocorr: float = 0.0


class SpatialMoDEData:
    """
    Spatial multi-omics data container with coordinates and neighborhood graph.

    Parameters
    ----------
    rna : pd.DataFrame
        RNA count matrix (spots × genes).
    atac : pd.DataFrame
        ATAC count matrix (spots × peaks).
    obs : pd.DataFrame
        Spot metadata (condition, region, etc.).
    coords : pd.DataFrame
        Spot coordinates (columns: x, y).
    regions : pd.Series, optional
        Anatomical region labels per spot.
    """

    def __init__(
        self,
        rna: pd.DataFrame,
        atac: pd.DataFrame,
        obs: pd.DataFrame,
        coords: pd.DataFrame,
        regions: pd.Series | None = None,
    ):
        self.rna = rna
        self.atac = atac
        self.obs = obs
        self.coords = coords
        self.regions = regions
        self.spatial_graph: csr_matrix | None = None
        self.n_neighbors: int = 6

    @property
    def n_spots(self) -> int:
        return self.rna.shape[0]

    def build_graph(
        self,
        method: str = "knn",
        n_neighbors: int = 6,
        radius: float | None = None,
    ) -> csr_matrix:
        """
        Build spatial neighborhood graph from coordinates.

        Parameters
        ----------
        method : str
            "knn" or "radius".
        n_neighbors : int
            Number of nearest neighbors (knn mode).
        radius : float
            Distance threshold (radius mode).
        """
        from sklearn.neighbors import NearestNeighbors
        coords_arr = self.coords[["x", "y"]].values

        if method == "knn":
            nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, self.n_spots))
            nn.fit(coords_arr)
            distances, indices = nn.kneighbors(coords_arr)
            # Build adjacency (exclude self)
            rows, cols, data = [], [], []
            for i in range(self.n_spots):
                for j_idx in range(1, len(indices[i])):  # skip self
                    j = indices[i][j_idx]
                    if j < self.n_spots:
                        rows.append(i); cols.append(j); data.append(1.0)
                        rows.append(j); cols.append(i); data.append(1.0)
            self.spatial_graph = csr_matrix(
                (data, (rows, cols)), shape=(self.n_spots, self.n_spots)
            )
        elif method == "radius" and radius is not None:
            nn = NearestNeighbors(radius=radius)
            nn.fit(coords_arr)
            adj = nn.radius_neighbors_graph(coords_arr, mode="connectivity")
            self.spatial_graph = adj.astype(float).tocsr()

        self.n_neighbors = n_neighbors
        return self.spatial_graph

    def compute_neighbor_effect(
        self, event_values: np.ndarray
    ) -> np.ndarray:
        """
        Compute per-spot neighbor-averaged event values.

        event_values: (n_spots,) array of per-spot event activity.
        Returns: (n_spots,) neighbor-weighted average.
        """
        if self.spatial_graph is None:
            self.build_graph()
        rowsum = np.asarray(self.spatial_graph.sum(axis=1)).flatten()
        rowsum[rowsum == 0] = 1
        neighbor_avg = self.spatial_graph.dot(event_values) / rowsum
        return neighbor_avg

    def compute_moran_i(self, values: np.ndarray) -> float:
        """Compute Moran's I spatial autocorrelation."""
        n = len(values)
        if self.spatial_graph is None:
            self.build_graph()
        W = self.spatial_graph
        w_sum = W.sum()
        if w_sum == 0:
            return 0.0
        mean = values.mean()
        dev = values - mean
        num = n * dev.dot(W.dot(dev))
        den = w_sum * dev.dot(dev)
        return float(num / den) if den > 0 else 0.0

    def compute_edge_score(self) -> np.ndarray:
        """Compute edge/boundary artifact score per spot."""
        if self.spatial_graph is None:
            self.build_graph()
        degrees = np.asarray(self.spatial_graph.sum(axis=1)).flatten()
        median_deg = np.median(degrees)
        edge = np.where(degrees < median_deg * 0.5, 1.0, 0.0)
        return edge

    def compute_spatial_evidence(
        self, event_values: np.ndarray
    ) -> SpatialEvidence:
        """Compute full spatial evidence for an event."""
        return SpatialEvidence(
            moran_i=self.compute_moran_i(event_values),
            neighbor_effect=float(self.compute_neighbor_effect(event_values).mean()),
            edge_score=float(self.compute_edge_score().mean()),
            region_autocorr=0.0,  # requires region labels
        )
