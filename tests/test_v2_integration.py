"""v2.0 integration tests: full pipeline with new modalities."""

import pytest
import numpy as np
import pandas as pd


class TestCUTTAGPipeline:
    """End-to-end: RNA + CUT&Tag through MoDES."""

    def test_epigenomic_pipeline_runs(self):
        from modes.data import MoDEData
        from modes import MoDES

        rng = np.random.default_rng(42)
        n = 10
        condition = np.array(["ctrl"] * (n // 2) + ["trt"] * (n // 2))
        obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n)])

        gene_names = [f"g{i}:chr1:{1000000+i*500000}" for i in range(10)]
        rna = pd.DataFrame(
            rng.poisson(200, (n, 10)), index=obs.index, columns=gene_names,
        )

        peak_names = [f"chr1:{990000+i*500000}-{1010000+i*500000}" for i in range(5)]
        atac = pd.DataFrame(
            rng.poisson(100, (n, 5)), index=obs.index, columns=peak_names,
        )

        epi_names = [
            f"chr1:{990000+i*500000}-{1010000+i*500000}|H3K27ac" for i in range(3)
        ]
        epi = pd.DataFrame(
            rng.poisson(60, (n, 3)), index=obs.index, columns=epi_names,
        )
        epi_feat = pd.DataFrame({
            "feature_id": epi_names,
            "chr": ["chr1"] * 3,
            "start": [990000 + i * 500000 for i in range(3)],
            "end": [1010000 + i * 500000 for i in range(3)],
            "assay": ["CUTTAG"] * 3,
            "target": ["H3K27ac"] * 3,
        })

        # Add CUT&Tag effect in trt
        epi.iloc[n // 2:, 0] = rng.poisson(300, n // 2)

        tss_map = {}
        for g in gene_names:
            parts = g.split(":chr")
            name = parts[0]
            chrom = "chr" + parts[1].split(":")[0]
            pos = int(parts[1].split(":")[1])
            tss_map[g] = (name, chrom, pos)

        data = MoDEData.from_epigenomic_matrices(
            rna_counts=rna,
            epigenomic_counts=epi,
            epigenomic_features=epi_feat,
            metadata=obs,
            condition_col="condition",
            target="H3K27ac",
            atac_counts=atac,
        )

        assert data.n_samples == n
        assert "h3k27ac_cuttag" in data.modalities
        assert "h3k27ac_cuttag" in data.modality_specs
        assert data.modality_specs["h3k27ac_cuttag"].target == "H3K27ac"
        assert data.modality_specs["h3k27ac_cuttag"].is_activating()

        # Run MoDES on RNA+ATAC (cuttag modality available but not yet in core pipeline)
        modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
        result = modes.run()
        assert len(result.event_table) > 0
        assert "state" in result.event_table.columns


class TestProteinPipeline:
    """End-to-end: RNA + ATAC + Protein through MoDES."""

    def test_protein_pipeline_runs(self):
        from modes.data import MoDEData
        from modes import MoDES

        rng = np.random.default_rng(42)
        n = 10
        condition = np.array(["ctrl"] * (n // 2) + ["trt"] * (n // 2))
        obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n)])

        gene_names = [f"g{i}:chr1:{1000000+i*500000}" for i in range(10)]
        rna = pd.DataFrame(
            rng.poisson(200, (n, 10)), index=obs.index, columns=gene_names,
        )
        peak_names = [f"chr1:{990000+i*500000}-{1010000+i*500000}" for i in range(5)]
        atac = pd.DataFrame(
            rng.poisson(100, (n, 5)), index=obs.index, columns=peak_names,
        )
        prot_names = ["CD4", "CD8A", "CD19"]
        prot = pd.DataFrame(
            rng.poisson(80, (n, 3)), index=obs.index, columns=prot_names,
        )
        # Spike protein effect
        prot.iloc[n // 2:, 0] = rng.poisson(400, n // 2)

        prot_links = pd.DataFrame({
            "protein_id": prot_names,
            "gene": gene_names[:3],
        })

        tss_map = {}
        for g in gene_names:
            parts = g.split(":chr")
            name = parts[0]
            chrom = "chr" + parts[1].split(":")[0]
            pos = int(parts[1].split(":")[1])
            tss_map[g] = (name, chrom, pos)

        data = MoDEData.from_protein_matrices(
            rna_counts=rna,
            atac_counts=atac,
            protein_counts=prot,
            protein_gene_links=prot_links,
            metadata=obs,
            condition_col="condition",
        )

        assert data.n_samples == n
        assert "protein" in data.modalities
        assert "protein" in data.modality_specs
        assert data.modalities["protein"].shape == (n, 3)

        # Run MoDES on RNA+ATAC (protein modality available in modalities dict)
        modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
        result = modes.run()
        assert len(result.event_table) > 0
