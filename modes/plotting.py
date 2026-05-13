"""Visualization functions for MoDES results."""

from __future__ import annotations

import io
import base64
from typing import Optional

import numpy as np
import pandas as pd

STATE_COLORS = {
    "concordant": "#2ca02c",
    "chromatin_primed": "#ff7f0e",
    "rna_only": "#1f77b4",
    "discordant_opposite": "#d62728",
    "null": "#7f7f7f",
}


def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def plot_state_distribution(results) -> "plt.Figure":
    """Bar chart of event counts per state."""
    import matplotlib.pyplot as plt

    state_counts = results.event_table["state"].value_counts()

    fig, ax = plt.subplots(figsize=(6, 4))
    colors = [STATE_COLORS.get(s, "#7f7f7f") for s in state_counts.index]
    bars = ax.bar(state_counts.index, state_counts.values, color=colors)

    ax.set_xlabel("State")
    ax.set_ylabel("Number of events")
    ax.set_title("Event State Distribution")

    for bar, count in zip(bars, state_counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            str(count), ha="center", va="bottom", fontsize=9,
        )

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    return fig


def plot_volcano_atac(results, fdr_threshold: float = 0.1) -> "plt.Figure":
    """Volcano plot of ATAC effects."""
    import matplotlib.pyplot as plt

    df = results.event_table.copy()
    df["neg_log10_pval"] = -np.log10(df["atac_pval"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 5))

    sig_mask = df["atac_fdr"] < fdr_threshold
    nonsig = df[~sig_mask]

    ax.scatter(
        nonsig["atac_coef"], nonsig["neg_log10_pval"],
        c="#cccccc", s=8, alpha=0.5, label=f"NS (FDR >= {fdr_threshold})",
    )

    for state in ["concordant", "chromatin_primed", "rna_only", "discordant_opposite"]:
        submask = sig_mask & (df["state"] == state)
        if submask.sum() == 0:
            continue
        sub = df[submask]
        ax.scatter(
            sub["atac_coef"], sub["neg_log10_pval"],
            c=STATE_COLORS[state], s=15, alpha=0.8, label=state,
        )

    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5)
    ax.set_xlabel("ATAC effect (log fold change)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("ATAC Volcano Plot")
    ax.legend(fontsize="small", markerscale=0.8)
    fig.tight_layout()
    return fig


def plot_volcano_rna(results, fdr_threshold: float = 0.1) -> "plt.Figure":
    """Volcano plot of RNA effects."""
    import matplotlib.pyplot as plt

    df = results.event_table.copy()
    df["neg_log10_pval"] = -np.log10(df["rna_pval"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 5))

    sig_mask = df["rna_fdr"] < fdr_threshold
    nonsig = df[~sig_mask]

    ax.scatter(
        nonsig["rna_coef"], nonsig["neg_log10_pval"],
        c="#cccccc", s=8, alpha=0.5, label=f"NS (FDR >= {fdr_threshold})",
    )

    for state in ["concordant", "chromatin_primed", "rna_only", "discordant_opposite"]:
        submask = sig_mask & (df["state"] == state)
        if submask.sum() == 0:
            continue
        sub = df[submask]
        ax.scatter(
            sub["rna_coef"], sub["neg_log10_pval"],
            c=STATE_COLORS[state], s=15, alpha=0.8, label=state,
        )

    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5)
    ax.set_xlabel("RNA effect (log fold change)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("RNA Volcano Plot")
    ax.legend(fontsize="small", markerscale=0.8)
    fig.tight_layout()
    return fig


def plot_atac_vs_rna(results, color_by: str = "state") -> "plt.Figure":
    """Scatter plot: ATAC effect vs RNA effect colored by state."""
    import matplotlib.pyplot as plt

    df = results.event_table.copy()

    fig, ax = plt.subplots(figsize=(7, 6))

    for state in ["concordant", "chromatin_primed", "rna_only", "discordant_opposite", "null"]:
        sub = df[df["state"] == state]
        if sub.empty:
            continue
        ax.scatter(
            sub["atac_coef"], sub["rna_coef"],
            c=STATE_COLORS.get(state, "#7f7f7f"),
            s=15, alpha=0.7, label=state,
        )

    # Diagonal: x = y
    lim = max(abs(df["atac_coef"]).max(), abs(df["rna_coef"]).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", alpha=0.3, linewidth=1)
    ax.axhline(0, color="grey", alpha=0.3, linewidth=0.5)
    ax.axvline(0, color="grey", alpha=0.3, linewidth=0.5)

    ax.set_xlabel("ATAC effect")
    ax.set_ylabel("RNA effect")
    ax.set_title("ATAC vs RNA Effect")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.legend(fontsize="small", markerscale=0.8)
    ax.set_aspect("equal")
    fig.tight_layout()
    return fig


def plot_evidence_heatmap(results, n_top: int = 50, figsize=(8, 10)) -> "plt.Figure":
    """Heatmap of evidence vectors for top events."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    if results.evidence_vectors is None:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No evidence vectors available", ha="center", va="center")
        return fig

    df = results.evidence_vectors.copy()

    # Select top events by evidence magnitude
    ev_cols = ["z_atac", "z_rna", "z_rna_given_atac"]
    df["_score"] = df[ev_cols].abs().sum(axis=1)
    df = df.nlargest(n_top, "_score").sort_values("_score", ascending=False)

    heatmap_data = df[ev_cols + ["quality_score"]].set_index(df["event_id"])

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        heatmap_data, ax=ax, cmap="RdBu_r", center=0,
        cbar_kws={"label": "z-score"},
        xticklabels=True, yticklabels=len(heatmap_data) <= 30,
    )
    ax.set_title(f"Top {min(n_top, len(df))} Events Evidence Vectors")
    fig.tight_layout()
    return fig


def plot_quality_distribution(results) -> "plt.Figure":
    """Histogram of quality scores faceted by state."""
    import matplotlib.pyplot as plt

    df = results.event_table.copy()
    states = [s for s in ["concordant", "chromatin_primed", "rna_only", "discordant_opposite", "null"]
              if s in df["state"].values]

    n_states = len(states)
    fig, axes = plt.subplots(1, n_states, figsize=(3 * n_states, 3), sharey=True)
    if n_states == 1:
        axes = [axes]

    for ax, state in zip(axes, states):
        sub = df[df["state"] == state]
        ax.hist(sub["quality_score"], bins=20, color=STATE_COLORS.get(state, "#7f7f7f"), alpha=0.8)
        ax.axvline(0.3, color="red", linestyle="--", alpha=0.5, label="artifact threshold")
        ax.set_title(state)
        ax.set_xlabel("Quality score")
        if ax == axes[0]:
            ax.set_ylabel("Count")

    fig.suptitle("Quality Score Distribution by State")
    fig.tight_layout()
    return fig


def plot_network(results, min_confidence: float = 0.5) -> "plt.Figure":
    """Visualize event network (TF->peak->gene) with state coloring."""
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
    except ImportError:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "networkx not available", ha="center", va="center")
        return fig

    G = nx.Graph()

    top_events = results.event_table[
        results.event_table["confidence"] >= min_confidence
    ]

    # Show top events only
    if len(top_events) > 100:
        top_events = top_events.nlargest(50, "confidence")

    for _, row in top_events.iterrows():
        gene_node = f"g:{row['gene']}"
        peak_node = f"p:{row['peak_id']}"

        G.add_node(gene_node, ntype="gene", size=8)
        G.add_node(peak_node, ntype="peak", size=5)

        color = STATE_COLORS.get(row["state"], "#7f7f7f")
        G.add_edge(gene_node, peak_node, color=color, state=row["state"])

        if pd.notna(row.get("tf_name")) and row["tf_name"]:
            tf_node = f"tf:{row['tf_name']}"
            G.add_node(tf_node, ntype="tf", size=12)
            G.add_edge(tf_node, peak_node, color="#aaaaaa", style="dashed")

    fig, ax = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42, k=1.5)

    # Draw by node type
    for ntype, ncolor, nsize in [("tf", "#e377c2", 80), ("gene", "#17becf", 50), ("peak", "#bcbd22", 30)]:
        nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == ntype]
        if nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=ncolor,
                                   node_size=nsize, ax=ax, alpha=0.8)

    # Draw edges
    edge_colors = [d["color"] for _, _, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, alpha=0.4, ax=ax)

    # Labels (genes and TFs only)
    labels = {n: n.split(":", 1)[1] for n, d in G.nodes(data=True)
              if d.get("ntype") in ("gene", "tf")}
    nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax)

    ax.set_title("MoDES Event Network")
    ax.axis("off")
    fig.tight_layout()
    return fig
