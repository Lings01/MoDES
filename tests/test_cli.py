"""Tests for MoDES CLI."""

import os
import subprocess
import sys


def test_cli_help():
    """CLI should print help without crashing."""
    result = subprocess.run(
        [sys.executable, "-m", "modes.cli", "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "MoDES" in result.stdout


def test_cli_run_help():
    """CLI run subcommand should show help."""
    result = subprocess.run(
        [sys.executable, "-m", "modes.cli", "run", "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--rna" in result.stdout


def test_cli_validate_help():
    """CLI validate-input subcommand should show help."""
    result = subprocess.run(
        [sys.executable, "-m", "modes.cli", "validate-input", "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--rna" in result.stdout


def test_cli_no_command_shows_help():
    """CLI with no subcommand should print help and exit 1."""
    result = subprocess.run(
        [sys.executable, "-m", "modes.cli"],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "usage" in result.stdout or "MoDES" in result.stdout


def test_cli_validate_on_example_data():
    """CLI validate-input on the minimal example data."""
    base = os.path.join(os.path.dirname(__file__), "..", "examples", "minimal_bulk")
    result = subprocess.run(
        [
            sys.executable, "-m", "modes.cli", "validate-input",
            "--rna", os.path.join(base, "rna_counts.tsv"),
            "--atac", os.path.join(base, "atac_counts.tsv"),
            "--metadata", os.path.join(base, "metadata.tsv"),
            "--condition", "condition",
            "--external-links", os.path.join(base, "peak_gene_links.tsv"),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0


def test_cli_entry_point():
    """Verify the main() function is importable."""
    from modes.cli import main
    assert main is not None
