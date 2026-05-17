"""Event candidate construction: link peaks to genes."""

from __future__ import annotations

import numpy as np
import pandas as pd

from modes._types import EventCandidate


class EventCandidateBuilder:
    """
    Build candidate regulatory events linking peaks to target genes.

    Supports promoter peaks, distal peaks, motif annotation,
    and external pre-computed links.

    Parameters
    ----------
    promoter_window : int
        bp around TSS considered as promoter (default 2000).
    distal_window : int
        bp around TSS for distal candidate search (default 250000).
    """

    def __init__(
        self,
        promoter_window: int = 2000,
        distal_window: int = 250000,
    ):
        self.promoter_window = promoter_window
        self.distal_window = distal_window
        self._tss_map: dict[str, tuple] | None = None

    def build(
        self,
        gene_names: list[str],
        peak_names: list[str],
        external_links: pd.DataFrame | None = None,
        motif_annotation: pd.DataFrame | None = None,
        genome_annotation: str | None = None,
        tss_map: dict[str, tuple] | None = None,
    ) -> pd.DataFrame:
        """
        Build candidate events.

        Parameters
        ----------
        gene_names : list of str
            Gene identifiers.
        peak_names : list of str
            Peak identifiers in format 'chr:start-end' or 'chr_start_end'.
        external_links : DataFrame, optional
            Columns: ['peak_id', 'gene', 'tf_name', 'source'].
            Pre-computed peak-to-gene links from SCENIC+/SCARlink/ArchR.
        motif_annotation : DataFrame, optional
            Columns: ['peak_id', 'tf_name'].
        genome_annotation : str, optional
            Path to GTF/GFF for TSS annotation.
        tss_map : dict, optional
            Manual gene -> (name, chr, tss_position) mapping.

        Returns
        -------
        DataFrame of EventCandidate columns.
        """
        rows = []

        import warnings as _warnings

        # Parse peak coordinates
        peak_coords = [_parse_peak_name(p) for p in peak_names]
        peak_df = pd.DataFrame(
            peak_coords, columns=["peak_id", "chr", "start", "end"]
        )
        n_unknown_peaks = (peak_df["chr"] == "unknown").sum()
        if n_unknown_peaks > 0:
            _warnings.warn(
                f"{n_unknown_peaks}/{len(peak_names)} peaks could not be parsed "
                "as genomic intervals and may be excluded from coordinate-based "
                "event generation.",
                UserWarning,
            )

        # Parse gene TSS from names if format gene:chr:pos is used,
        # or from genome annotation, or from manual map
        if tss_map is not None:
            self._tss_map = tss_map
        elif genome_annotation:
            self._tss_map = _parse_gtf_tss(genome_annotation, gene_names)
        else:
            self._tss_map = _parse_gene_tss_from_names(gene_names)

        # Warn about genes with no genomic coordinates
        n_total_genes = len(gene_names)
        n_missing_genes = 0
        for g in gene_names:
            tss_info = self._tss_map.get(g)
            if tss_info is None:
                # Try lookup by parsed short name (e.g., "geneX" from "geneX:chr1:1000")
                for key, val in self._tss_map.items():
                    if val[0] == g or key == g.split(":")[0].split("_")[0]:
                        tss_info = val
                        break
            if tss_info is None:
                n_missing_genes += 1
                continue
            chrom = tss_info[1] if len(tss_info) > 1 else ""
            if chrom in {"", "unknown", None}:
                n_missing_genes += 1
        if n_missing_genes > 0:
            _warnings.warn(
                f"{n_missing_genes}/{n_total_genes} genes have no genomic "
                "coordinates and may be excluded from coordinate-based event "
                "generation. Provide genome_annotation, tss_map, or "
                "external_links for better coverage.",
                UserWarning,
            )

        # 1. Promoter and distal peaks per gene
        # P0 opt: Build per-chromosome interval index for O(G log P) construction
        peak_index = {}
        for chrom, sub in peak_df.groupby("chr"):
            sub = sub.sort_values("start").reset_index(drop=True)
            peak_index[chrom] = {
                "starts": sub["start"].to_numpy(),
                "ends": sub["end"].to_numpy(),
                "df": sub,
            }

        for gene in gene_names:
            tss_info = self._tss_map.get(gene)
            if tss_info is None:
                for key, val in self._tss_map.items():
                    if val[0] == gene or key == gene.split(":")[0].split("_")[0]:
                        tss_info = val
                        break
            if tss_info is None:
                continue
            _, tss_chr, tss_pos = tss_info

            pi = peak_index.get(tss_chr)
            if pi is None:
                continue

            # Binary search for peaks within TSS ± distal_window
            left = tss_pos - self.distal_window
            right = tss_pos + self.distal_window
            lo = np.searchsorted(pi["starts"], left, side="left")
            hi = np.searchsorted(pi["starts"], right, side="right")
            candidate_peaks = pi["df"].iloc[lo:hi]

            for _, peak in candidate_peaks.iterrows():
                distance = _peak_tss_distance(
                    peak["start"], peak["end"], tss_pos
                )
                abs_dist = abs(distance)

                if abs_dist <= self.promoter_window:
                    source = "promoter"
                elif abs_dist <= self.distal_window:
                    source = "distal_250kb"
                else:
                    continue

                from hashlib import sha1
                event_key = f"{peak['peak_id']}|{gene}|{source}"
                event_id = sha1(event_key.encode()).hexdigest()[:12]
                rows.append(
                    EventCandidate(
                        event_id=event_id,
                        gene=gene,
                        peak_id=peak["peak_id"],
                        peak_chr=peak["chr"],
                        peak_start=peak["start"],
                        peak_end=peak["end"],
                        context="",
                        link_source=source,
                        distance_to_tss=distance,
                    )
                )

        candidates = pd.DataFrame(
            [r.__dict__ for r in rows], columns=_candidate_columns()
        )

        # 2. Merge external links
        if external_links is not None:
            candidates = self._merge_external_links(candidates, external_links)

        # 3. Assign motifs
        if motif_annotation is not None:
            candidates = self._assign_motifs(candidates, motif_annotation)

        # P0 opt: dedup by biological key (peak_id + gene + tf_name), not event_id
        if not candidates.empty:
            candidates["_bio_key"] = (
                candidates["peak_id"] + "|" + candidates["gene"] + "|"
                + candidates["tf_name"].fillna("")
            )
            candidates.drop_duplicates(subset=["_bio_key"], inplace=True)
            candidates.drop(columns=["_bio_key"], inplace=True)
            candidates.reset_index(drop=True, inplace=True)

        return candidates

    def _merge_external_links(
        self,
        candidates: pd.DataFrame,
        external: pd.DataFrame,
    ) -> pd.DataFrame:
        """Merge pre-computed peak-gene links, avoiding duplicates."""
        ext = external.copy()
        required = ["peak_id", "gene"]
        for col in required:
            if col not in ext.columns:
                raise ValueError(f"External links missing required column '{col}'")

        # P0 opt: dedup by biological key (peak_id + gene + tf_name)
        existing_keys = set(
            candidates["peak_id"] + "|" + candidates["gene"] + "|"
            + candidates["tf_name"].fillna("")
        )
        ext_records = []
        for _, row in ext.iterrows():
            row_tf = str(row.get("tf_name", "")) if pd.notna(row.get("tf_name")) else ""
            bio_key = f"{row['peak_id']}|{row['gene']}|{row_tf}"
            if bio_key in existing_keys:
                mask = (
                    (candidates["peak_id"] == row["peak_id"])
                    & (candidates["gene"] == row["gene"])
                    & (candidates["tf_name"].fillna("") == row_tf)
                )
                if "tf_name" in ext.columns and pd.notna(row.get("tf_name")):
                    candidates.loc[mask, "tf_name"] = row["tf_name"]
                # Append source
                candidates.loc[mask, "link_source"] = candidates.loc[mask, "link_source"].astype(str) + ";external"
                continue

            from hashlib import sha1
            eid = sha1(f"{row['peak_id']}|{row['gene']}|external".encode()).hexdigest()[:12]
            ext_records.append(
                EventCandidate(
                    event_id=eid,
                    gene=row["gene"],
                    peak_id=row["peak_id"],
                    peak_chr=row.get("peak_chr", ""),
                    peak_start=int(row.get("peak_start", 0)),
                    peak_end=int(row.get("peak_end", 0)),
                    tf_name=row.get("tf_name"),
                    context=row.get("context", ""),
                    link_source=row.get("source", "external"),
                    distance_to_tss=int(row.get("distance_to_tss", 0)),
                )
            )

        if ext_records:
            ext_df = pd.DataFrame(
                [r.__dict__ for r in ext_records], columns=_candidate_columns()
            )
            candidates = pd.concat([candidates, ext_df], ignore_index=True)

        return candidates

    def _assign_motifs(
        self,
        candidates: pd.DataFrame,
        motifs: pd.DataFrame,
    ) -> pd.DataFrame:
        """Assign TF names to peaks with motif hits."""
        if "peak_id" not in motifs.columns or "tf_name" not in motifs.columns:
            raise ValueError("motif_annotation needs 'peak_id' and 'tf_name' columns")

        motif_map = motifs.groupby("peak_id")["tf_name"].apply(
            lambda x: "|".join(x.dropna().unique())
        )
        candidates["tf_name"] = candidates["tf_name"].fillna(
            candidates["peak_id"].map(motif_map)
        )
        return candidates


def _candidate_columns() -> list[str]:
    return [
        "event_id", "gene", "peak_id", "peak_chr", "peak_start",
        "peak_end", "tf_name", "context", "link_source", "distance_to_tss",
    ]


def _parse_peak_name(name: str) -> tuple:
    """Parse 'chr1:100-200' or 'chr1_100_200' format."""
    import re

    # Try chr:start-end
    m = re.match(r"(chr\w+)[:\-_](\d+)[:\-_](\d+)", str(name))
    if m:
        return (str(name), m.group(1), int(m.group(2)), int(m.group(3)))
    # Fallback
    return (str(name), "unknown", 0, 0)


def _parse_gene_tss_from_names(gene_names: list[str]) -> dict[str, tuple]:
    """
    Try to parse gene names that include coordinate info.
    Falls back to a dummy map if names don't contain coordinates.
    """
    import re

    tss_map = {}
    for g in gene_names:
        m = re.match(r"(.+?)[:\-_](chr\w+)[:\-_](\d+)", str(g))
        if m:
            tss_map[m.group(1)] = (m.group(1), m.group(2), int(m.group(3)))
        else:
            # Placeholder: gene acts as its own name, no coordinate info
            tss_map[str(g)] = (str(g), "", 0)
    return tss_map


def _parse_gtf_tss(gtf_path: str, gene_names: list[str]) -> dict[str, tuple]:
    """Parse TSS positions from a GTF file for requested genes."""
    tss_map = {}
    gene_set = set(gene_names)

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene" and parts[2] != "transcript":
                continue

            # Extract gene name from attributes
            attrs = {}
            for a in parts[8].split(";"):
                a = a.strip()
                if not a:
                    continue
                if " " in a:
                    key, val = a.split(" ", 1)
                    attrs[key] = val.strip('"')
                elif "=" in a:
                    key, val = a.split("=", 1)
                    attrs[key] = val.strip('"')

            gname = attrs.get("gene_name") or attrs.get("gene_id") or attrs.get("ID")
            if gname and gname in gene_set:
                chrom = parts[0]
                if parts[6] == "+":
                    tss = int(parts[3])
                else:
                    tss = int(parts[4])
                tss_map[gname] = (gname, chrom, tss)

    return tss_map


def _peak_tss_distance(peak_start: int, peak_end: int, tss_pos: int) -> int:
    """Distance from peak center to TSS. Negative = peak upstream of TSS."""
    peak_center = (peak_start + peak_end) // 2
    return peak_center - tss_pos


def validate_external_links(links: pd.DataFrame, peak_names: set, gene_names: set) -> list[str]:
    """
    Validate an external_links DataFrame.

    Returns a list of issue strings. Empty list means valid.
    """
    issues = []
    required = ["peak_id", "gene"]
    for col in required:
        if col not in links.columns:
            issues.append(f"Missing required column: {col}")
    if issues:
        return issues

    n_total = len(links)
    n_peak_match = links["peak_id"].isin(peak_names).sum()
    n_gene_match = links["gene"].isin(gene_names).sum()
    if n_peak_match < n_total:
        issues.append(
            f"Only {n_peak_match}/{n_total} peak_ids found in ATAC matrix"
        )
    if n_gene_match < n_total:
        issues.append(
            f"Only {n_gene_match}/{n_total} genes found in RNA matrix"
        )

    dupes = links.duplicated(subset=["peak_id", "gene"]).sum()
    if dupes > 0:
        issues.append(f"{dupes} duplicate peak_id-gene pairs found")

    if "score" in links.columns:
        bad_scores = ((links["score"] < 0) | (links["score"] > 1)).sum()
        if bad_scores > 0:
            issues.append(f"{bad_scores} scores outside [0, 1] range")

    return issues
