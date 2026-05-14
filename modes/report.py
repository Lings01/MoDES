"""HTML report generation for MoDES results."""

from __future__ import annotations

import io
import base64
import html
from datetime import datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from modes.core import MoDESResult


def _esc(x):
    """HTML-escape a value."""
    return html.escape(str(x))


ALLOWED_STATES = {
    "concordant", "chromatin_primed", "rna_only",
    "discordant_opposite", "null",
}


TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MoDES Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
         max-width: 1100px; margin: 0 auto; padding: 20px; color: #333;
         background: #fafafa; }}
  h1 {{ color: #1a1a2e; border-bottom: 3px solid #1f77b4; padding-bottom: 10px; }}
  h2 {{ color: #16213e; margin-top: 30px; border-bottom: 1px solid #ddd; }}
  table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }}
  th, td {{ padding: 6px 10px; text-align: left; border-bottom: 1px solid #ddd; }}
  th {{ background: #f0f0f0; font-weight: 600; }}
  tr:hover {{ background: #f5f5ff; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
                   gap: 12px; margin: 15px 0; }}
  .summary-card {{ background: white; border: 1px solid #e0e0e0; border-radius: 8px;
                   padding: 15px; text-align: center; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
  .summary-card .value {{ font-size: 28px; font-weight: bold; color: #1f77b4; }}
  .summary-card .label {{ font-size: 12px; color: #666; margin-top: 4px; }}
  .state-bar {{ display: inline-block; padding: 2px 8px; border-radius: 4px;
                font-size: 12px; font-weight: 600; }}
  .concordant {{ background: #2ca02c33; color: #1b6e1b; }}
  .chromatin_primed {{ background: #ff7f0e33; color: #b85b00; }}
  .rna_only {{ background: #1f77b433; color: #0d3b6e; }}
  .discordant_opposite {{ background: #d6272833; color: #8b1a1a; }}
  .null {{ background: #7f7f7f33; color: #555; }}
  .artifact-risk-low {{ color: #2e7d32; }}
  .artifact-risk-medium {{ color: #ef6c00; }}
  .artifact-risk-high {{ color: #c62828; font-weight: bold; }}
  img.plot {{ max-width: 100%; height: auto; margin: 10px 0;
              border: 1px solid #e0e0e0; border-radius: 4px; }}
  .footer {{ margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd;
             font-size: 12px; color: #999; }}
</style>
</head>
<body>
<h1>MoDES Analysis Report</h1>
<p>Generated: {timestamp}</p>

<h2>Summary</h2>
<div class="summary-grid">
  {summary_cards}
</div>

<h2>State Distribution</h2>
{state_distribution_table}
<img class="plot" src="data:image/png;base64,{plot_state_dist}" alt="State Distribution">
{artifact_summary_html}

<h2>Event Table (Top 50)</h2>
{event_table_html}

<h2>Volcano Plots</h2>
<img class="plot" src="data:image/png;base64,{plot_volcano_atac}" alt="ATAC Volcano">
<img class="plot" src="data:image/png;base64,{plot_volcano_rna}" alt="RNA Volcano">

<h2>ATAC vs RNA Effect</h2>
<img class="plot" src="data:image/png;base64,{plot_atac_vs_rna}" alt="ATAC vs RNA">

<h2>Run Parameters</h2>
<table>
  <tr><th>Parameter</th><th>Value</th></tr>
  {params_rows}
</table>

<div class="footer">
  MoDES v0.1.0 &mdash; Multi-Omics Discordance/Event State inference
</div>
</body>
</html>"""


def generate_report(results: "MoDESResult", output_path: str) -> None:
    """
    Generate a self-contained HTML report.

    Parameters
    ----------
    results : MoDESResult
    output_path : str
        Path to write report.html.
    """
    from modes.plotting import (
        plot_state_distribution,
        plot_volcano_atac,
        plot_volcano_rna,
        plot_atac_vs_rna,
    )

    # Generate plots
    fig_state = plot_state_distribution(results)
    fig_va = plot_volcano_atac(results)
    fig_vr = plot_volcano_rna(results)
    fig_avr = plot_atac_vs_rna(results)

    import matplotlib
    matplotlib.use("Agg")

    plot_state_dist = _fig_to_base64(fig_state)
    plot_volcano_atac_img = _fig_to_base64(fig_va)
    plot_volcano_rna_img = _fig_to_base64(fig_vr)
    plot_atac_vs_rna_img = _fig_to_base64(fig_avr)

    import matplotlib.pyplot as plt
    plt.close("all")

    # Summary cards
    state_counts = results.event_table["state"].value_counts()
    summary_cards = ""
    cards_data = [
        ("Total Events", str(len(results.event_table))),
        ("Concordant", str(state_counts.get("concordant", 0))),
        ("Chromatin Primed", str(state_counts.get("chromatin_primed", 0))),
        ("RNA Only", str(state_counts.get("rna_only", 0))),
        ("Null", str(state_counts.get("null", 0))),
        ("Samples", str(results.params.get("n_samples", "N/A"))),
        ("Genes", str(results.params.get("n_genes", "N/A"))),
    ]
    # Artifact risk distribution if available
    if "artifact_risk" in results.event_table.columns:
        risk_counts = results.event_table["artifact_risk"].value_counts()
        for risk_level in ["high", "medium", "low"]:
            if risk_level in risk_counts.index:
                cards_data.append(
                    (f"Artifact Risk: {risk_level}",
                     str(risk_counts.get(risk_level, 0)))
                )
    for label, value in cards_data:
        summary_cards += f'<div class="summary-card"><div class="value">{_esc(value)}</div><div class="label">{_esc(label)}</div></div>\n'

    # State distribution table
    state_rows = ""
    for state, count in state_counts.items():
        cls = state if state in ALLOWED_STATES else "null"
        state_rows += (
            f'<tr><td><span class="state-bar {cls}">{_esc(state)}</span></td>'
            f'<td>{count}</td>'
            f'<td>{count / len(results.event_table) * 100:.1f}%</td></tr>\n'
        )

    # Artifact risk table
    artifact_rows = ""
    if "artifact_risk" in results.event_table.columns:
        risk_counts = results.event_table["artifact_risk"].value_counts()
        for risk in ["low", "medium", "high"]:
            cnt = risk_counts.get(risk, 0)
            artifact_rows += (
                f'<tr><td><span class="artifact-risk-{risk}">{risk}</span></td>'
                f'<td>{cnt}</td>'
                f'<td>{cnt / len(results.event_table) * 100:.1f}%</td></tr>\n'
            )

    state_table_html = f"""
    <table>
      <tr><th>State</th><th>Count</th><th>Fraction</th></tr>
      {state_rows}
    </table>"""

    # Artifact risk summary (if available)
    artifact_summary_html = ""
    if artifact_rows:
        artifact_summary_html = f"""
        <h2>Artifact Risk Distribution</h2>
        <table>
          <tr><th>Risk Level</th><th>Count</th><th>Fraction</th></tr>
          {artifact_rows}
        </table>"""

    # Event table (top 50 by confidence)
    top_events = results.event_table.nlargest(50, "state_confidence")
    display_cols = [
        "event_id", "gene", "peak_id", "state", "state_confidence",
        "artifact_risk",
        "atac_coef", "rna_coef", "atac_fdr", "rna_fdr",
    ]
    avail_cols = [c for c in display_cols if c in top_events.columns]

    event_rows = ""
    for _, row in top_events[avail_cols].iterrows():
        state = row.get("state", "null")
        cls = state if state in ALLOWED_STATES else "null"
        event_rows += "<tr>"
        for col in avail_cols:
            val = row[col]
            if col == "state":
                event_rows += f'<td><span class="state-bar {cls}">{_esc(val)}</span></td>'
            elif isinstance(val, float):
                event_rows += f"<td>{val:.3f}</td>"
            else:
                event_rows += f"<td>{_esc(val)}</td>"
        event_rows += "</tr>\n"

    event_table_html = f"""
    <table>
      <tr>{"".join(f"<th>{c}</th>" for c in avail_cols)}</tr>
      {event_rows}
    </table>"""

    # Params
    params_rows = ""
    for k, v in results.params.items():
        params_rows += f"<tr><td>{_esc(k)}</td><td>{_esc(v)}</td></tr>\n"

    # Render
    html = TEMPLATE.format(
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        summary_cards=summary_cards,
        state_distribution_table=state_table_html,
        artifact_summary_html=artifact_summary_html,
        plot_state_dist=plot_state_dist,
        plot_volcano_atac=plot_volcano_atac_img,
        plot_volcano_rna=plot_volcano_rna_img,
        plot_atac_vs_rna=plot_atac_vs_rna_img,
        event_table_html=event_table_html,
        params_rows=params_rows,
    )

    with open(output_path, "w") as f:
        f.write(html)


def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")
