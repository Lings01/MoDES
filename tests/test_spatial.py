"""Tests for MoDES-Spatial native graph engine."""

import pytest
import numpy as np
import pandas as pd

from modes.modalities.spatial import SpatialMoDEData, SpatialEvidence


class TestSpatialMoDEData:
    @pytest.fixture
    def spatial_data(self):
        rng = np.random.default_rng(42)
        n = 30
        rna = pd.DataFrame(
            rng.poisson(100, (n, 5)),
            index=[f"s{i}" for i in range(n)],
            columns=[f"G{i}" for i in range(5)],
        )
        atac = pd.DataFrame(
            rng.poisson(50, (n, 4)),
            index=[f"s{i}" for i in range(n)],
            columns=[f"chr1:{1000+i*500}-{1500+i*500}" for i in range(4)],
        )
        obs = pd.DataFrame(
            {"condition": ["ctrl"] * 15 + ["trt"] * 15},
            index=[f"s{i}" for i in range(n)],
        )
        coords = pd.DataFrame({
            "x": rng.uniform(0, 100, n),
            "y": rng.uniform(0, 100, n),
        }, index=[f"s{i}" for i in range(n)])
        return SpatialMoDEData(rna=rna, atac=atac, obs=obs, coords=coords)

    def test_init(self, spatial_data):
        assert spatial_data.n_spots == 30
        assert spatial_data.rna.shape == (30, 5)
        assert spatial_data.atac.shape == (30, 4)

    def test_build_knn_graph(self, spatial_data):
        g = spatial_data.build_graph(method="knn", n_neighbors=4)
        assert g.shape == (30, 30)
        assert g.nnz > 0
        assert spatial_data.spatial_graph is not None

    def test_build_radius_graph(self, spatial_data):
        g = spatial_data.build_graph(method="radius", radius=30.0)
        assert g.shape == (30, 30)

    def test_compute_neighbor_effect(self, spatial_data):
        spatial_data.build_graph(method="knn", n_neighbors=4)
        values = np.random.randn(30)
        neighbor = spatial_data.compute_neighbor_effect(values)
        assert len(neighbor) == 30

    def test_compute_moran_i(self, spatial_data):
        spatial_data.build_graph(method="knn", n_neighbors=4)
        values = np.random.randn(30)
        mi = spatial_data.compute_moran_i(values)
        assert -1 <= mi <= 1 or np.isnan(mi)

    def test_compute_edge_score(self, spatial_data):
        spatial_data.build_graph(method="knn", n_neighbors=4)
        edge = spatial_data.compute_edge_score()
        assert len(edge) == 30
        assert set(np.unique(edge)).issubset({0.0, 1.0})

    def test_compute_spatial_evidence(self, spatial_data):
        spatial_data.build_graph(method="knn", n_neighbors=4)
        values = np.random.randn(30)
        ev = spatial_data.compute_spatial_evidence(values)
        assert isinstance(ev, SpatialEvidence)
        assert isinstance(ev.moran_i, float)
