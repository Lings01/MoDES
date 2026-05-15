"""MoDES command-line interface."""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="modes",
        description="MoDES: Multi-Omics Discordance/Event State inference",
    )
    subparsers = parser.add_subparsers(dest="command", help="subcommand")

    # modes run
    run_parser = subparsers.add_parser("run", help="Run full MoDES pipeline")
    run_parser.add_argument("--rna", required=True, help="RNA count matrix (TSV/CSV)")
    run_parser.add_argument("--atac", required=True, help="ATAC peak count matrix (TSV/CSV)")
    run_parser.add_argument("--metadata", required=True, help="Sample metadata (TSV/CSV)")
    run_parser.add_argument("--condition", required=True, help="Condition column name")
    run_parser.add_argument("--external-links", default=None, help="Peak-gene links (TSV)")
    run_parser.add_argument("--genome-annotation", default=None, help="GTF/GFF file")
    run_parser.add_argument("--donor", default=None, help="Donor column name")
    run_parser.add_argument("--batch", default=None, help="Batch column name")
    run_parser.add_argument("--covariates", default=None, help="Comma-separated covariate columns")
    run_parser.add_argument("--fdr-threshold", type=float, default=0.1, help="FDR threshold")
    run_parser.add_argument("--out", default="output", help="Output directory")
    run_parser.add_argument("--report", action="store_true", help="Generate HTML report")
    run_parser.add_argument("--network", action="store_true", help="Generate GraphML network")

    # modes validate-input
    val_parser = subparsers.add_parser("validate-input", help="Validate input files")
    val_parser.add_argument("--rna", required=True, help="RNA count matrix (TSV/CSV)")
    val_parser.add_argument("--atac", required=True, help="ATAC peak count matrix (TSV/CSV)")
    val_parser.add_argument("--metadata", required=True, help="Sample metadata (TSV/CSV)")
    val_parser.add_argument("--condition", required=True, help="Condition column name")
    val_parser.add_argument("--external-links", default=None, help="Peak-gene links (TSV)")
    val_parser.add_argument("--out", default=None, help="Output report path")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "run":
        _cmd_run(args)
    elif args.command == "validate-input":
        _cmd_validate(args)
    else:
        parser.print_help()
        sys.exit(1)


def _cmd_run(args):
    import pandas as pd
    from modes import MoDES, MoDEData

    covariates = args.covariates.split(",") if args.covariates else None

    data = MoDEData.from_matrices(
        rna_counts=args.rna,
        atac_counts=args.atac,
        metadata=args.metadata,
        condition_col=args.condition,
        donor_col=args.donor,
        batch_col=args.batch,
        index_col=0,
    )

    links = None
    if args.external_links:
        links = pd.read_csv(args.external_links, sep="\t")

    modes = MoDES(
        data=data,
        condition_col=args.condition,
        covariate_cols=covariates,
        donor_col=args.donor,
        batch_col=args.batch,
        fdr_threshold=args.fdr_threshold,
        genome_annotation=args.genome_annotation,
        external_links=links,
    )

    result = modes.run()
    print(result.summary())
    result.to_tsv(args.out)

    if args.report:
        result.to_report(os.path.join(args.out, "report.html"))
    if args.network:
        result.to_graphml(os.path.join(args.out, "event_network.graphml"))

    print(f"\nOutput written to: {args.out}")


def _cmd_validate(args):
    import pandas as pd
    from modes.data import MoDEData

    issues = []
    warnings_list = []

    try:
        data = MoDEData.from_matrices(
            rna_counts=args.rna,
            atac_counts=args.atac,
            metadata=args.metadata,
            condition_col=args.condition,
            index_col=0,
        )
        issues.extend(data.validate())

        cond = data.obs[args.condition]
        categories = sorted(set(cond))
        if len(categories) != 2:
            issues.append(
                f"Condition '{args.condition}' has {len(categories)} categories: "
                f"{categories}. MoDES requires exactly 2."
            )

        if args.external_links:
            links = pd.read_csv(args.external_links, sep="\t")
            for col in ["peak_id", "gene"]:
                if col not in links.columns:
                    issues.append(f"External links missing required column: {col}")
            n_matched = links["peak_id"].isin(data.peak_names).sum()
            n_total = len(links)
            if n_matched < n_total:
                warnings_list.append(
                    f"Only {n_matched}/{n_total} external links match ATAC peaks"
                )
            n_gene_matched = links["gene"].isin(data.gene_names).sum()
            if n_gene_matched < n_total:
                warnings_list.append(
                    f"Only {n_gene_matched}/{n_total} external links match RNA genes"
                )

    except Exception as e:
        issues.append(f"Failed to load data: {e}")

    report = {
        "ok": len(issues) == 0,
        "errors": issues,
        "warnings": warnings_list,
    }
    # Add summary if data was loaded
    try:
        report["n_samples"] = data.n_samples
        report["n_genes"] = data.n_genes
        report["n_peaks"] = data.n_peaks
    except Exception:
        pass
    try:
        report["n_links"] = n_total
        report["n_links_matched"] = n_matched
    except Exception:
        pass

    if args.out:
        ext = os.path.splitext(args.out)[1]
        if ext == ".json":
            import json
            with open(args.out, "w") as f:
                json.dump(report, f, indent=2)
        else:
            with open(args.out, "w") as f:
                f.write("MoDES Input Validation Report\n")
                f.write("=" * 40 + "\n")
                f.write(f"OK: {report['ok']}\n")
                for i in issues:
                    f.write(f"[FAIL] {i}\n")
                for w in warnings_list:
                    f.write(f"[WARN] {w}\n")
                if not issues and not warnings_list:
                    f.write("All checks passed.\n")

    if issues:
        print("ISSUES:")
        for i in issues:
            print(f"  [FAIL] {i}")
    if warnings_list:
        print("WARNINGS:")
        for w in warnings_list:
            print(f"  [WARN] {w}")
    if not issues and not warnings_list:
        print("All checks passed.")


if __name__ == "__main__":
    main()
