from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from .pipeline import (
    CommandRunner,
    PipelineSettings,
    STAGE_ORDER,
    TBASPipeline,
    check_dependencies,
    required_tools,
)


def _parse_stages(value: str) -> list[str]:
    raw = value.strip()
    if raw.lower() == "all":
        return list(STAGE_ORDER)
    stages = [stage.strip() for stage in raw.split(",") if stage.strip()]
    if not stages:
        raise argparse.ArgumentTypeError(
            "Stage list is empty. Use comma-separated names or 'all'."
        )
    return stages


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the TBAS software pipeline locally (no sbatch generation). "
            "Pipeline starts from demultiplexing calls*.bam."
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Path to sample manifest CSV/TSV.",
    )
    parser.add_argument(
        "--stages",
        default="all",
        type=_parse_stages,
        help=(
            "Comma-separated stages to run or 'all'. "
            f"Default order: {', '.join(STAGE_ORDER)}"
        ),
    )
    parser.add_argument(
        "--output-folder",
        type=Path,
        default=Path("GREGoR_adaptive_sampling"),
        help="Root folder containing per-sample output subfolders.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("GRCh38-2.1.0_genome_mainchrs.fa"),
        help="Reference FASTA path.",
    )
    parser.add_argument(
        "--kit-name",
        default="SQK-NBD114-96",
        help="ONT kit name used in read-group tags.",
    )
    parser.add_argument(
        "--clair3-model",
        type=Path,
        default=Path("clair3/bin/models/r1041_e82_400bps_sup_v520"),
        help="Path to Clair3 model directory.",
    )
    parser.add_argument(
        "--adotto-pheno",
        type=Path,
        default=Path("adaptive_regions/patho_tr/adotto_pheno.bed"),
        help="Pathogenic TR BED used by medaka_patho stage.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only, without executing.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print commands before execution.",
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help=(
            "Check that all external tools required by the selected stages are "
            "on PATH, report the result, and exit without running the pipeline."
        ),
    )
    parser.add_argument(
        "--skip-dep-check",
        action="store_true",
        help="Skip the preflight tool check that runs before a real execution.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    stages = args.stages
    missing = check_dependencies(stages)

    if args.check_deps:
        print(f"Required tools for selected stages: {', '.join(required_tools(stages))}")
        if missing:
            print(f"MISSING from PATH: {', '.join(missing)}")
            return 1
        print("All required tools were found on PATH.")
        return 0

    # Fail fast before a long run if a real execution is missing a tool.
    if missing and not args.dry_run and not args.skip_dep_check:
        print(
            "Error: the following required tools are not on PATH: "
            f"{', '.join(missing)}.\n"
            "Install them (see environment.yml) or pass --skip-dep-check to "
            "override.",
            file=sys.stderr,
        )
        return 1

    settings = PipelineSettings(
        manifest=args.manifest,
        output_folder=args.output_folder,
        reference=args.reference,
        kit_name=args.kit_name,
        clair3_model=args.clair3_model,
        adotto_pheno=args.adotto_pheno,
    )
    runner = CommandRunner(dry_run=args.dry_run, print_commands=not args.quiet)
    pipeline = TBASPipeline.from_manifest(settings, runner=runner)
    pipeline.run(stages)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
