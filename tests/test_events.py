"""Tests for EventCandidateBuilder."""

import pandas as pd

from modes.events import EventCandidateBuilder, _parse_peak_name


class TestEventCandidateBuilder:
    def test_build_basic(self):
        """Test basic event building with promoter and distal peaks."""
        builder = EventCandidateBuilder(promoter_window=2000, distal_window=250000)

        # Gene named with coordinates: geneA:chr1:10000
        gene_names = ["geneA:chr1:10000", "geneB:chr1:50000"]
        peak_names = [
            "chr1:9000-11000",    # promoter of geneA
            "chr1:10000-12000",   # promoter of geneA (overlapping TSS)
            "chr1:20000-22000",   # distal of geneA (10kb from TSS)
            "chr1:50000-52000",   # promoter of geneB
            "chr1:300000-302000", # outside distal range
            "chr2:5000-7000",     # wrong chromosome
        ]

        events = builder.build(gene_names=gene_names, peak_names=peak_names)

        assert len(events) > 0
        # Should have events for both genes
        assert any("geneA" in eid for eid in events["event_id"])
        assert any("geneB" in eid for eid in events["event_id"])

    def test_build_creates_unique_event_ids(self):
        builder = EventCandidateBuilder()
        gene_names = ["geneX:chr1:1000"]
        peak_names = ["chr1:800-1200"]

        events = builder.build(gene_names=gene_names, peak_names=peak_names)

        assert len(events) == 1
        assert "geneX" in events.iloc[0]["event_id"]
        assert "chr1:800-1200" in events.iloc[0]["event_id"]

    def test_build_with_external_links(self):
        """Test merging external pre-computed links."""
        builder = EventCandidateBuilder()
        gene_names = ["geneX:chr1:1000"]
        peak_names = ["chr1:800-1200"]

        external = pd.DataFrame({
            "peak_id": ["chr1:2000-3000"],
            "gene": ["geneX:chr1:1000"],
            "tf_name": ["STAT1"],
            "source": ["scenic"],
        })

        events = builder.build(
            gene_names=gene_names,
            peak_names=peak_names,
            external_links=external,
        )

        # Should include the external link
        assert any("external" in eid for eid in events["event_id"])

    def test_build_with_motif_annotation(self):
        """Test assigning TF motifs to peaks."""
        builder = EventCandidateBuilder()
        gene_names = ["geneX:chr1:1000"]
        peak_names = ["chr1:800-1200"]

        motifs = pd.DataFrame({
            "peak_id": ["chr1:800-1200"],
            "tf_name": ["IRF1"],
        })

        events = builder.build(
            gene_names=gene_names,
            peak_names=peak_names,
            motif_annotation=motifs,
        )

        assert events.iloc[0]["tf_name"] == "IRF1"

    def test_build_empty_no_peaks(self):
        """Graceful handling when no peaks match."""
        builder = EventCandidateBuilder()
        gene_names = ["geneX:chr1:1000"]
        peak_names = ["chr2:800-1200"]  # wrong chromosome

        events = builder.build(gene_names=gene_names, peak_names=peak_names)
        assert len(events) == 0


class TestParsePeakName:
    def test_standard_format(self):
        result = _parse_peak_name("chr1:100-200")
        assert result[0] == "chr1:100-200"
        assert result[1] == "chr1"
        assert result[2] == 100
        assert result[3] == 200

    def test_underscore_format(self):
        result = _parse_peak_name("chr1_100_200")
        assert result[1] == "chr1"
        assert result[2] == 100
        assert result[3] == 200

    def test_fallback(self):
        result = _parse_peak_name("weird_format")
        assert result[1] == "unknown"
