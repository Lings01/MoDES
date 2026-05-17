"""MoDES orchestrator and result container."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from modes._types import EventResult
from modes.data import MoDEData
from modes.decompose import ConditionalDecomposition
from modes.effects import EffectEstimator
from modes.events import EventCandidateBuilder
from modes.states import EvidenceBuilder, StateClassifier


def _event_result_columns():
    return [
        "event_id", "tf_name", "peak_id", "gene", "context",
        "state", "state_confidence", "quality_score",
        "atac_coef", "atac_se", "atac_pval", "atac_fdr", "atac_direction",
        "rna_coef", "rna_se", "rna_pval", "rna_fdr", "rna_direction",
        "rna_after_atac_coef", "rna_after_atac_se",
        "rna_after_atac_pval", "rna_after_atac_fdr",
        "artifact_risk", "artifact_reason",
        "event_pval", "event_fdr",
    ]


def _read_optional_tsv(output_dir: str, filename: str):
    """Read a TSV file if it exists, else return None."""
    path = os.path.join(output_dir, filename)
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    return None


class MoDES:
    """
    Multi-Omics Discordance/Event State inference.

    Main orchestrator for the MoDES-RA pipeline (RNA + ATAC).

    Parameters
    ----------
    data : MoDEData
        Input data with RNA and ATAC matrices.
    condition_col : str
        Column in data.obs specifying the condition of interest.
    covariate_cols : list of str, optional
        Additional covariates.
    donor_col : str, optional
        Donor/replicate column.
    batch_col : str, optional
        Batch column.
    fdr_threshold : float
        FDR cutoff for significance. Default 0.1.
    genome_annotation : str, optional
        Path to GTF for TSS annotation.
    external_links : pd.DataFrame, optional
        Pre-computed peak-to-gene links.
    motif_annotation : pd.DataFrame, optional
        Peak-to-TF motif mapping.
    """

    def __init__(
        self,
        data: MoDEData,
        condition_col: str,
        covariate_cols: list[str] | None = None,
        donor_col: str | None = None,
        batch_col: str | None = None,
        fdr_threshold: float = 0.1,
        genome_annotation: str | None = None,
        external_links: pd.DataFrame | None = None,
        motif_annotation: pd.DataFrame | None = None,
        tss_map: dict | None = None,
        contrast: tuple | None = None,
        allow_poisson_fallback: bool = True,
        allow_simplified_fallback: bool = False,
        conditional_mode: str = "auto",
        cov_type: str = "nonrobust",
        min_nonzero_samples: int = 3,
        min_total_count: float = 10.0,
        batch_col_quality: str | None = None,
    ):
        self.data = data
        self.condition_col = condition_col
        self.covariate_cols = covariate_cols or []
        self.donor_col = donor_col
        self.batch_col = batch_col
        self.fdr_threshold = fdr_threshold
        self.genome_annotation = genome_annotation
        self.external_links = external_links
        self.motif_annotation = motif_annotation
        self.tss_map = tss_map
        self.contrast = contrast
        self.allow_poisson_fallback = allow_poisson_fallback
        self.allow_simplified_fallback = allow_simplified_fallback
        self.conditional_mode = conditional_mode
        self.cov_type = cov_type
        self.min_nonzero_samples = min_nonzero_samples
        self.min_total_count = min_total_count
        self.batch_col_quality = batch_col_quality

        # Pipeline state — v2.0: generic multi-modality
        self.events: pd.DataFrame | None = None
        self.effects: dict[str, dict] = {}     # modality_name → {feature_id → ModalityEffect}
        self.atac_effects: dict | None = None   # backward compat alias
        self.rna_effects: dict | None = None    # backward compat alias
        self.conditional_effects: pd.DataFrame | None = None
        self.evidence: pd.DataFrame | None = None
        self.states: pd.DataFrame | None = None
        self.results: MoDESResult | None = None
        self.modality_evidence: pd.DataFrame | None = None  # v2.0 long-format

    def run(self) -> MoDESResult:
        """Execute the full MoDES pipeline."""
        self.build_events()
        self.estimate_effects()
        self.decompose()
        self.build_evidence()
        self.classify_states()
        self.results = self._assemble_results()
        self._build_modality_evidence()  # v2.0: long-format multi-modal evidence
        return self.results

    def build_events(
        self,
        external_links: pd.DataFrame | None = None,
        motif_annotation: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """
        Step 1: Build candidate regulatory events.

        Parameters
        ----------
        external_links : DataFrame, optional
            Override instance's external_links.
        motif_annotation : DataFrame, optional
            Override instance's motif_annotation.

        Returns
        -------
        DataFrame of EventCandidate records.
        """
        builder = EventCandidateBuilder(
            promoter_window=2000,
            distal_window=250000,
        )

        ext_links = self.external_links if external_links is None else external_links
        motifs = self.motif_annotation if motif_annotation is None else motif_annotation
        self.events = builder.build(
            gene_names=list(self.data.gene_names),
            peak_names=list(self.data.peak_names),
            external_links=ext_links,
            motif_annotation=motifs,
            genome_annotation=self.genome_annotation,
            tss_map=self.tss_map,
        )
        if self.events is None or len(self.events) == 0:
            raise ValueError(
                "No candidate events were generated. "
                "For normal gene symbols, provide one of: "
                "genome_annotation, tss_map, or external_links."
            )
        return self.events

    def estimate_effects(self) -> tuple[dict, dict]:
        """
        Step 2: Estimate ATAC and RNA condition effects.

        Returns
        -------
        (atac_effects, rna_effects)
        """
        if self.events is None:
            raise RuntimeError("Call build_events() first")

        estimator = EffectEstimator(
            condition_col=self.condition_col,
            covariate_cols=self.covariate_cols,
            donor_col=self.donor_col,
            batch_col=self.batch_col,
            use_empirical_bayes=True,
            contrast=self.contrast,
            allow_poisson_fallback=self.allow_poisson_fallback,
            allow_simplified_fallback=self.allow_simplified_fallback,
            cov_type=self.cov_type,
            min_nonzero_samples=self.min_nonzero_samples,
            min_total_count=self.min_total_count,
        )

        peak_names = list(self.events["peak_id"].unique())
        gene_names = list(self.events["gene"].unique())

        self.atac_effects, self.rna_effects = estimator.estimate_effects(
            self.data, peak_names, gene_names
        )
        # v2.0: populate generic effects dict
        self.effects["atac"] = self.atac_effects
        self.effects["rna"] = self.rna_effects
        return self.atac_effects, self.rna_effects

    def decompose(self) -> pd.DataFrame:
        """
        Step 3: Conditional decomposition (RNA after ATAC).

        Returns
        -------
        DataFrame of conditional effects.
        """
        if self.atac_effects is None or self.rna_effects is None:
            raise RuntimeError("Call estimate_effects() first")

        decomposer = ConditionalDecomposition(
            condition_col=self.condition_col,
            covariate_cols=self.covariate_cols,
            donor_col=self.donor_col,
            batch_col=self.batch_col,
            contrast=self.contrast,
            conditional_mode=self.conditional_mode,
        )

        self.conditional_effects = decomposer.decompose(
            data=self.data,
            events=self.events,
            atac_effects=self.atac_effects,
            rna_effects=self.rna_effects,
        )
        return self.conditional_effects

    def build_evidence(self) -> pd.DataFrame:
        """
        Step 4: Build evidence vectors.

        Returns
        -------
        DataFrame of evidence vectors.
        """
        if self.conditional_effects is None:
            raise RuntimeError("Call decompose() first")

        builder = EvidenceBuilder(
            batch_col=self.batch_col_quality or self.batch_col
        )
        self.evidence = builder.build(
            events=self.events,
            atac_effects=self.atac_effects,
            rna_effects=self.rna_effects,
            conditional_effects=self.conditional_effects,
            data=self.data,
        )
        return self.evidence

    def classify_states(self) -> pd.DataFrame:
        """
        Step 5: Classify events into regulatory states.

        Returns
        -------
        DataFrame with state assignments and confidence scores.
        """
        if self.evidence is None:
            raise RuntimeError("Call build_evidence() first")

        classifier = StateClassifier(
            fdr_threshold=self.fdr_threshold,
            use_empirical_bayes=True,
        )

        self.states = classifier.classify(self.evidence)
        return self.states

    def _assemble_results(self) -> MoDESResult:
        """Combine all pipeline outputs into MoDESResult."""
        if self.events is None or len(self.events) == 0:
            params = {
                "condition_col": self.condition_col,
                "covariate_cols": self.covariate_cols,
                "donor_col": self.donor_col,
                "batch_col": self.batch_col,
                "fdr_threshold": self.fdr_threshold,
                "n_events": 0,
                "n_samples": self.data.n_samples,
                "n_genes": self.data.n_genes,
                "n_peaks": self.data.n_peaks,
            }
            return MoDESResult(
                event_table=pd.DataFrame(columns=_event_result_columns()),
                params=params,
            )

        # --- Pass 0: Pre-build O(1) lookup maps (P0 opt: was O(E²) per-event filter) ---
        cond_map = {}
        for _, cr in self.conditional_effects.iterrows():
            cond_map[cr["event_id"]] = cr
        state_map = {}
        for _, sr in self.states.iterrows():
            state_map[sr["event_id"]] = sr
        ev_map = {}
        for _, er in self.evidence.iterrows():
            ev_map[er["event_id"]] = er

        # --- Pass 1: collect per-event data and compute event_pval ---
        raw_data = []
        from modes.utils import benjamini_hochberg

        for _, event in self.events.iterrows():
            eid = event["event_id"]
            peak = event["peak_id"]
            gene = event["gene"]

            atac = self.atac_effects.get(peak)
            rna = self.rna_effects.get(gene)

            # Conditional effect (O(1) dict lookup)
            cr = cond_map.get(eid)
            if cr is None:
                rna_after_coef = np.nan
                rna_after_se = np.nan
                rna_after_pval = 1.0
                rna_after_fdr = 1.0
            else:
                rna_after_coef = cr.get("rna_after_atac_coef", np.nan)
                rna_after_se = cr.get("rna_after_atac_se", np.nan)
                rna_after_pval = cr.get("rna_after_atac_pval", 1.0)
                rna_after_fdr = cr.get("rna_after_atac_fdr", 1.0)

            # State (O(1) dict lookup)
            sr = state_map.get(eid)
            if sr is None:
                state = "null"
                confidence = 1.0
                artifact_risk = "low"
                artifact_reason = ""
            else:
                state = sr["state"]
                confidence = sr["state_confidence"]
                artifact_risk = sr.get("artifact_risk", "low")
                artifact_reason = sr.get("artifact_reason", "")

            # Quality (O(1) dict lookup)
            er = ev_map.get(eid)
            quality = er["quality_score"] if er is not None else 0.5

            # Event-level p-value based on state
            atac_pval = atac.p_value if atac else 1.0
            rna_pval = rna.p_value if rna else 1.0
            if state == "concordant":
                event_pval = max(atac_pval, rna_pval)
            elif state == "discordant_opposite":
                event_pval = max(atac_pval, rna_pval)
            elif state == "chromatin_primed":
                event_pval = atac_pval
            elif state == "rna_only":
                event_pval = rna_pval
            else:
                event_pval = 1.0

            raw_data.append({
                "event_id": eid,
                "peak_id": peak,
                "gene": gene,
                "tf_name": event.get("tf_name"),
                "context": event.get("context", ""),
                "atac": atac,
                "rna": rna,
                "rna_after_coef": rna_after_coef,
                "rna_after_se": rna_after_se,
                "rna_after_pval": rna_after_pval,
                "rna_after_fdr": rna_after_fdr,
                "state": state,
                "state_confidence": confidence,
                "artifact_risk": artifact_risk,
                "artifact_reason": artifact_reason,
                "quality": quality,
                "event_pval": event_pval,
            })

        # Event-level BH correction
        event_pvals = np.array([d["event_pval"] for d in raw_data])
        event_fdrs = benjamini_hochberg(event_pvals)

        # --- Pass 2: build EventResult records ---
        records = []
        for d, event_fdr in zip(raw_data, event_fdrs):
            atac = d["atac"]
            rna = d["rna"]
            records.append(
                EventResult(
                    event_id=d["event_id"],
                    tf_name=d["tf_name"],
                    peak_id=d["peak_id"],
                    gene=d["gene"],
                    context=d["context"],
                    atac_coef=atac.coef if atac else np.nan,
                    atac_se=atac.se if atac else np.nan,
                    atac_pval=atac.p_value if atac else 1.0,
                    atac_fdr=atac.fdr if atac else 1.0,
                    atac_direction=atac.direction if atac else 0,
                    rna_coef=rna.coef if rna else np.nan,
                    rna_se=rna.se if rna else np.nan,
                    rna_pval=rna.p_value if rna else 1.0,
                    rna_fdr=rna.fdr if rna else 1.0,
                    rna_direction=rna.direction if rna else 0,
                    rna_after_atac_coef=d["rna_after_coef"],
                    rna_after_atac_se=d["rna_after_se"],
                    rna_after_atac_pval=d["rna_after_pval"],
                    rna_after_atac_fdr=d["rna_after_fdr"],
                    state=d["state"],
                    state_confidence=d["state_confidence"],
                    quality_score=d["quality"],
                    artifact_risk=d["artifact_risk"],
                    artifact_reason=d["artifact_reason"],
                    event_pval=d["event_pval"],
                    event_fdr=float(event_fdr),
                )
            )

        import sys

        import numpy as _np
        import pandas as _pd
        import statsmodels as _sm

        from modes import __version__ as _modes_ver

        params = {
            "modes_version": _modes_ver,
            "python_version": sys.version.split()[0],
            "numpy_version": _np.__version__,
            "pandas_version": _pd.__version__,
            "statsmodels_version": _sm.__version__,
            "condition_col": self.condition_col,
            "contrast": str(self.contrast) if self.contrast else "auto",
            "covariate_cols": self.covariate_cols,
            "donor_col": self.donor_col,
            "batch_col": self.batch_col,
            "fdr_threshold": self.fdr_threshold,
            "allow_poisson_fallback": self.allow_poisson_fallback,
            "allow_simplified_fallback": self.allow_simplified_fallback,
            "cov_type": self.cov_type,
            "conditional_mode": self.conditional_mode,
            "min_nonzero_samples": self.min_nonzero_samples,
            "min_total_count": self.min_total_count,
            "n_events": len(records),
            "n_samples": self.data.n_samples,
            "n_genes": self.data.n_genes,
            "n_peaks": self.data.n_peaks,
            "n_external_links": len(self.external_links) if self.external_links is not None else 0,
        }

        model_diag = self._build_model_diagnostics()

        return MoDESResult(
            event_table=pd.DataFrame([r.__dict__ for r in records]),
            state_probabilities=self.states.copy() if self.states is not None else None,
            layer_effects=self._build_layer_effects_df(),
            evidence_vectors=self.evidence.copy() if self.evidence is not None else None,
            model_diagnostics=model_diag,
            modality_evidence=self.modality_evidence.copy() if self.modality_evidence is not None else None,
            params=params,
        )

    def _build_layer_effects_df(self) -> pd.DataFrame:
        """Assemble per-event layer effect estimates."""
        records = []
        for _, event in self.events.iterrows():
            eid = event["event_id"]
            peak = event["peak_id"]
            gene = event["gene"]

            atac = self.atac_effects.get(peak)
            rna = self.rna_effects.get(gene)

            records.append({
                "event_id": eid,
                "peak_id": peak,
                "gene": gene,
                "atac_coef": atac.coef if atac else np.nan,
                "atac_se": atac.se if atac else np.nan,
                "atac_z": atac.z_score if atac else np.nan,
                "atac_fdr": atac.fdr if atac else 1.0,
                "rna_coef": rna.coef if rna else np.nan,
                "rna_se": rna.se if rna else np.nan,
                "rna_z": rna.z_score if rna else np.nan,
                "rna_fdr": rna.fdr if rna else 1.0,
            })
        return pd.DataFrame(records)

    def _build_modality_evidence(self) -> pd.DataFrame:
        """v2.0: Build long-format event × modality evidence table."""
        rows = []
        for _, event in self.events.iterrows():
            eid = event["event_id"]
            gene = event["gene"]
            peak = event["peak_id"]

            # RNA evidence
            rna_eff = self.rna_effects.get(gene)
            if rna_eff:
                rows.append({
                    "event_id": eid, "modality": "rna", "assay": "RNA",
                    "target": None, "feature_id": gene, "role": "transcript_output",
                    "coef": rna_eff.coef, "se": rna_eff.se,
                    "pval": rna_eff.p_value, "fdr": rna_eff.fdr,
                    "direction": rna_eff.direction,
                    "quality_score": float(rna_eff.model_summary.get("quality_score", 0.5)) if isinstance(rna_eff.model_summary, dict) else 0.5,
                    "model_used": rna_eff.model_summary.get("model_used", "unknown") if isinstance(rna_eff.model_summary, dict) else "unknown",
                    "converged": rna_eff.convergence,
                })

            # ATAC evidence
            atac_eff = self.atac_effects.get(peak)
            if atac_eff:
                rows.append({
                    "event_id": eid, "modality": "atac", "assay": "ATAC",
                    "target": None, "feature_id": peak, "role": "chromatin_accessibility",
                    "coef": atac_eff.coef, "se": atac_eff.se,
                    "pval": atac_eff.p_value, "fdr": atac_eff.fdr,
                    "direction": atac_eff.direction,
                    "quality_score": float(atac_eff.model_summary.get("quality_score", 0.5)) if isinstance(atac_eff.model_summary, dict) else 0.5,
                    "model_used": atac_eff.model_summary.get("model_used", "unknown") if isinstance(atac_eff.model_summary, dict) else "unknown",
                    "converged": atac_eff.convergence,
                })

            # Additional modalities from data.modalities
            for mod_name in self.data.modalities:
                if mod_name in ("rna", "atac"):
                    continue
                spec = self.data.modality_specs.get(mod_name)
                eff_dict = self.effects.get(mod_name, {})
                feature = gene if spec and spec.feature_type == "gene" else peak
                mod_eff = eff_dict.get(feature)
                if mod_eff:
                    rows.append({
                        "event_id": eid, "modality": mod_name,
                        "assay": spec.assay if spec else "unknown",
                        "target": spec.target if spec else None,
                        "feature_id": feature,
                        "role": spec.regulatory_role if spec else "unknown",
                        "coef": mod_eff.coef, "se": mod_eff.se,
                        "pval": mod_eff.p_value, "fdr": mod_eff.fdr,
                        "direction": mod_eff.direction,
                        "quality_score": 0.5,
                        "model_used": mod_eff.model_summary.get("model_used", "unknown") if isinstance(mod_eff.model_summary, dict) else "unknown",
                        "converged": mod_eff.convergence,
                    })

        self.modality_evidence = pd.DataFrame(rows) if rows else pd.DataFrame()
        return self.modality_evidence

    def _build_model_diagnostics(self) -> pd.DataFrame:
        """Collect model diagnostics from all effect estimates."""
        rows = []
        for peak_id, e in (self.atac_effects or {}).items():
            s = e.model_summary or {}
            rows.append({
                "feature_id": peak_id,
                "modality": "ATAC",
                "model_used": s.get("model_used", "unknown"),
                "family": s.get("family", "unknown"),
                "alpha": s.get("alpha"),
                "alpha_estimated": s.get("alpha_estimated", False),
                "converged": s.get("converged", e.convergence),
                "dropped_covariates": s.get("dropped_covariates", False),
                "warning": s.get("warning", ""),
            })
        for gene_name, e in (self.rna_effects or {}).items():
            s = e.model_summary or {}
            rows.append({
                "feature_id": gene_name,
                "modality": "RNA",
                "model_used": s.get("model_used", "unknown"),
                "family": s.get("family", "unknown"),
                "alpha": s.get("alpha"),
                "alpha_estimated": s.get("alpha_estimated", False),
                "converged": s.get("converged", e.convergence),
                "dropped_covariates": s.get("dropped_covariates", False),
                "warning": s.get("warning", ""),
            })
        return pd.DataFrame(rows)


class MoDESResult:
    """
    Container for MoDES output.

    Parameters
    ----------
    event_table : DataFrame
        Main output with per-event effect sizes and state classifications.
    state_probabilities : DataFrame
        State confidence probability per event.
    layer_effects : DataFrame
        Per-layer effect size estimates.
    evidence_vectors : DataFrame
        D_e evidence vectors for all events.
    params : dict
        Run parameters.
    """

    def __init__(
        self,
        event_table: pd.DataFrame,
        state_probabilities: pd.DataFrame | None = None,
        layer_effects: pd.DataFrame | None = None,
        evidence_vectors: pd.DataFrame | None = None,
        model_diagnostics: pd.DataFrame | None = None,
        modality_evidence: pd.DataFrame | None = None,
        params: dict | None = None,
    ):
        self.event_table = event_table
        self.state_probabilities = state_probabilities
        self.layer_effects = layer_effects
        self.evidence_vectors = evidence_vectors
        self.model_diagnostics = model_diagnostics
        self.modality_evidence = modality_evidence
        self.params = params or {}

    def summary(self) -> str:
        """Return a text summary of the results."""
        lines = []
        lines.append("=" * 60)
        lines.append("MoDES Results Summary")
        lines.append("=" * 60)
        lines.append(f"Total events:     {len(self.event_table)}")
        lines.append(f"Significant (FDR < 0.1): "
                     f"{(self.event_table['atac_fdr'] < 0.1).sum()} ATAC, "
                     f"{(self.event_table['rna_fdr'] < 0.1).sum()} RNA")
        lines.append("")

        state_counts = self.event_table["state"].value_counts()
        lines.append("State distribution:")
        for state, count in state_counts.items():
            lines.append(f"  {state:20s}: {count:6d} "
                        f"({count / len(self.event_table) * 100:5.1f}%)")

        n_with_tf = self.event_table["tf_name"].notna().sum()
        if n_with_tf > 0:
            lines.append(f"\nEvents with TF annotation: {n_with_tf}")

        return "\n".join(lines)

    def filter(
        self,
        state: str | None = None,
        states: list[str] | None = None,
        min_confidence: float = 0.0,
        fdr_threshold: float | None = None,
        min_atac_fdr: float | None = None,
        min_rna_fdr: float | None = None,
        exclude_high_artifact: bool = False,
        max_artifact_risk: str | None = None,
        max_event_fdr: float | None = None,
        min_quality_score: float | None = None,
        genes: list[str] | None = None,
        peaks: list[str] | None = None,
        context: str | None = None,
    ) -> pd.DataFrame:
        """
        Filter the event table.

        Parameters
        ----------
        state : str, optional
            Keep only events with this state.
        states : list of str, optional
            Keep events with any of these states.
        min_confidence : float
            Minimum state confidence.
        fdr_threshold : float, optional
            Keep events where atac_fdr < threshold OR rna_fdr < threshold.
        exclude_high_artifact : bool
            Exclude events with artifact_risk == "high".
        max_artifact_risk : str, optional
            Maximum allowed artifact risk.
        max_event_fdr : float, optional
            Maximum event-level FDR threshold.
        min_quality_score : float, optional
            Minimum quality score.
        genes : list of str, optional
            Keep only events for these genes.
        peaks : list of str, optional
            Keep only events for these peaks.
        context : str, optional
            Keep only events matching this context.

        Returns
        -------
        Filtered DataFrame.
        """
        df = self.event_table.copy()

        if state is not None:
            df = df[df["state"] == state]

        if states is not None:
            df = df[df["state"].isin(states)]

        if min_confidence > 0:
            df = df[df["state_confidence"] >= min_confidence]

        if fdr_threshold is not None:
            df = df[
                (df["atac_fdr"] < fdr_threshold) |
                (df["rna_fdr"] < fdr_threshold)
            ]

        if min_atac_fdr is not None:
            df = df[df["atac_fdr"] < min_atac_fdr]

        if min_rna_fdr is not None:
            df = df[df["rna_fdr"] < min_rna_fdr]

        if max_event_fdr is not None and "event_fdr" in df.columns:
            df = df[df["event_fdr"] <= max_event_fdr]

        if exclude_high_artifact and "artifact_risk" in df.columns:
            df = df[df["artifact_risk"] != "high"]

        if max_artifact_risk is not None and "artifact_risk" in df.columns:
            risk_order = {"low": 0, "medium": 1, "high": 2}
            df = df[df["artifact_risk"].map(risk_order).fillna(0) <= risk_order.get(max_artifact_risk, 2)]

        if min_quality_score is not None and "quality_score" in df.columns:
            df = df[df["quality_score"] >= min_quality_score]

        if genes is not None:
            df = df[df["gene"].isin(genes)]

        if peaks is not None:
            df = df[df["peak_id"].isin(peaks)]

        if context is not None:
            df = df[df["context"] == context]

        return df.reset_index(drop=True)

    def save(self, output_dir: str) -> None:
        """Save all result tables to a directory. Alias for to_tsv()."""
        self.to_tsv(output_dir)
        # Also save params as JSON
        import json
        with open(os.path.join(output_dir, "run_params.json"), "w") as f:
            json.dump(self.params, f, indent=2, default=str)

    @classmethod
    def load(cls, output_dir: str) -> MoDESResult:
        """Load results from a directory previously written by save()/to_tsv()."""
        et = pd.read_csv(os.path.join(output_dir, "event_table.tsv"), sep="\t")
        # Restore string columns that pandas may have read as float64 (e.g., "null" → NaN)
        for col in ["state", "artifact_risk", "artifact_reason", "tf_name", "context"]:
            if col in et.columns:
                et[col] = et[col].fillna("").astype(str).replace({"nan": "", "None": ""})
                if col in ("state", "artifact_risk"):
                    et[col] = et[col].replace({"": "null", "nan": "null"})
        sp = _read_optional_tsv(output_dir, "event_state_confidence.tsv")
        le = _read_optional_tsv(output_dir, "event_layer_effects.tsv")
        ev = _read_optional_tsv(output_dir, "event_evidence_vectors.tsv")
        md = _read_optional_tsv(output_dir, "model_diagnostics.tsv")
        me = _read_optional_tsv(output_dir, "event_modality_evidence.tsv")
        import json
        params = {}
        rp_json = os.path.join(output_dir, "run_params.json")
        if os.path.exists(rp_json):
            with open(rp_json) as f:
                params = json.load(f)
        elif os.path.exists(os.path.join(output_dir, "run_params.tsv")):
            rp = pd.read_csv(os.path.join(output_dir, "run_params.tsv"), sep="\t")
            params = dict(zip(rp["parameter"], rp["value"]))
        return cls(
            event_table=et, state_probabilities=sp, layer_effects=le,
            evidence_vectors=ev, model_diagnostics=md,
            modality_evidence=me, params=params,
        )

    def to_tsv(self, output_dir: str) -> None:
        """Write TSV output files to output_dir."""
        os.makedirs(output_dir, exist_ok=True)

        self.event_table.to_csv(
            os.path.join(output_dir, "event_table.tsv"),
            sep="\t", index=False,
        )

        if self.state_probabilities is not None:
            self.state_probabilities.to_csv(
                os.path.join(output_dir, "event_state_confidence.tsv"),
                sep="\t", index=False,
            )

        if self.layer_effects is not None:
            self.layer_effects.to_csv(
                os.path.join(output_dir, "event_layer_effects.tsv"),
                sep="\t", index=False,
            )

        if self.evidence_vectors is not None:
            self.evidence_vectors.to_csv(
                os.path.join(output_dir, "event_evidence_vectors.tsv"),
                sep="\t", index=False,
            )

        # v2.0: long-format modality evidence
        if hasattr(self, "modality_evidence") and self.modality_evidence is not None and len(self.modality_evidence) > 0:
            self.modality_evidence.to_csv(
                os.path.join(output_dir, "event_modality_evidence.tsv"),
                sep="\t", index=False,
            )

        # Write model diagnostics
        if hasattr(self, "model_diagnostics") and self.model_diagnostics is not None:
            self.model_diagnostics.to_csv(
                os.path.join(output_dir, "model_diagnostics.tsv"),
                sep="\t", index=False,
            )

        # Write params as JSON-like TSV
        pd.DataFrame(list(self.params.items()), columns=["parameter", "value"]).to_csv(
            os.path.join(output_dir, "run_params.tsv"),
            sep="\t", index=False,
        )

    def to_graphml(self, path: str) -> None:
        """
        Write event network as GraphML.

        Nodes: genes, peaks, TFs.
        Edges: peak->gene (event link), TF->peak (motif support).
        """
        try:
            import networkx as nx
        except ImportError:
            raise ImportError("networkx is required for to_graphml()")

        G = nx.Graph()

        for _, row in self.event_table.iterrows():
            gene_node = f"gene:{row['gene']}"
            peak_node = f"peak:{row['peak_id']}"

            if gene_node not in G:
                G.add_node(gene_node, type="gene")
            if peak_node not in G:
                G.add_node(peak_node, type="peak")

            G.add_edge(
                gene_node, peak_node,
                type="event",
                state=row["state"],
                state_confidence=float(row.get("state_confidence", np.nan)),
                artifact_risk=str(row.get("artifact_risk", "low")),
                artifact_reason=str(row.get("artifact_reason", "")),
                event_pval=float(row.get("event_pval", 1.0)),
                event_fdr=float(row.get("event_fdr", 1.0)),
                atac_coef=row["atac_coef"],
                rna_coef=row["rna_coef"],
            )

            if pd.notna(row.get("tf_name")):
                tf_node = f"tf:{row['tf_name']}"
                if tf_node not in G:
                    G.add_node(tf_node, type="tf")
                G.add_edge(tf_node, peak_node, type="motif")

        nx.write_graphml(G, path)

    def to_report(self, path: str) -> None:
        """Generate HTML report."""
        from modes.report import generate_report
        generate_report(self, path)


def run_by_context(
    modes: MoDES,
    context_col: str = "context",
    min_samples_per_context: int = 4,
) -> pd.DataFrame:
    """
    Run MoDES separately for each context (e.g., cell type).

    Parameters
    ----------
    modes : MoDES
        Configured MoDES instance (without calling .run() yet).
    context_col : str
        Column in data.obs defining contexts (default: "context").
    min_samples_per_context : int
        Minimum samples required per context.

    Returns
    -------
    DataFrame with all contexts' events, with context column.
    """
    all_events = []
    contexts = modes.data.obs[context_col].unique()
    for ctx in contexts:
        ctx_mask = modes.data.obs[context_col] == ctx
        n_ctx = ctx_mask.sum()
        if n_ctx < min_samples_per_context:
            continue
        ctx_data = MoDEData(
            rna=modes.data.rna.loc[ctx_mask],
            atac=modes.data.atac.loc[ctx_mask],
            obs=modes.data.obs.loc[ctx_mask],
        )
        ctx_modes = MoDES(
            data=ctx_data,
            condition_col=modes.condition_col,
            covariate_cols=modes.covariate_cols,
            donor_col=modes.donor_col,
            batch_col=modes.batch_col,
            fdr_threshold=modes.fdr_threshold,
            external_links=modes.external_links,
            tss_map=modes.tss_map,
            contrast=modes.contrast,
        )
        result = ctx_modes.run()
        result.event_table["context"] = ctx
        all_events.append(result.event_table)
    if all_events:
        return pd.concat(all_events, ignore_index=True)
    return pd.DataFrame(columns=_event_result_columns())
