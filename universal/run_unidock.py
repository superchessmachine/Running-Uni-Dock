#!/usr/bin/env python3
"""
Universal Uni-Dock runner.

The script turns a single config file into a fully reproducible Uni-Dock
command. It understands both classic batch docking (one receptor with many
ligands) and the screening-test workflow that ships with Uni-Dock.  Paths in
the config are resolved relative to the config file location, so everything can
be driven from a single folder.
"""

from __future__ import annotations

import argparse
import inspect
import json
import shlex
import subprocess
import sys
import tarfile
from glob import glob
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Sequence


class ConfigError(RuntimeError):
    """Raised when the provided config file is incomplete or invalid."""


def _shlex_join(parts: Sequence[str]) -> str:
    """Return a shell-safe command string."""
    if hasattr(shlex, "join"):
        return shlex.join(parts)  # type: ignore[attr-defined]
    return " ".join(shlex.quote(part) for part in parts)


def _as_path(base_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _vector3(value, label: str) -> List[float]:
    """Parse a 3D vector from either a dict or a list."""
    if value is None:
        raise ConfigError(f"Missing {label} specification in box section.")
    entries: Sequence[float | int | str]
    if isinstance(value, dict):
        entries = [value.get(axis) for axis in ("x", "y", "z")]
    elif isinstance(value, (list, tuple)):
        entries = value
    else:
        raise ConfigError(
            f"{label} must be provided as a dict with x/y/z or a length-3 list."
        )
    if len(entries) != 3:
        raise ConfigError(f"{label} must contain exactly three values.")
    vector: List[float] = []
    for axis_value in entries:
        if axis_value is None:
            raise ConfigError(f"Axis missing for {label}.")
        try:
            vector.append(float(axis_value))
        except ValueError as exc:
            raise ConfigError(f"Could not convert {label} value '{axis_value}' to float") from exc
    return vector


def _unique(sequence: Iterable[Path]) -> List[Path]:
    """Preserve order while removing duplicates."""
    seen: set[Path] = set()
    result: List[Path] = []
    for item in sequence:
        normalized = item.resolve()
        if normalized in seen:
            continue
        seen.add(normalized)
        result.append(normalized)
    return result


@dataclass
class LigandCollection:
    paths: List[Path]
    index_file: Path


class UniDockRunner:
    """Build and run Uni-Dock commands from a config file."""

    def __init__(
        self,
        config_path: Path,
        skip_extract: bool = False,
        verbose: bool = False,
    ) -> None:
        self.config_path = config_path.resolve()
        self.config_dir = self.config_path.parent
        try:
            self.cfg = json.loads(self.config_path.read_text())
        except json.JSONDecodeError as exc:
            raise ConfigError(f"Failed to parse {self.config_path}: {exc}") from exc
        self._normalise_aliases()
        self.skip_extract = skip_extract
        self.verbose = verbose
        self.mode = self.cfg.get("mode", "batch").lower()
        if self.mode not in {"batch", "screening"}:
            raise ConfigError("mode must be 'batch' or 'screening'.")
        module_dir = Path(__file__).resolve().parent
        default_parent = module_dir / "runs"
        output_cfg = self.cfg.get("output", {})
        explicit_dir = _as_path(self.config_dir, output_cfg.get("dir"))
        parent_dir = _as_path(self.config_dir, output_cfg.get("parent_dir")) or default_parent
        run_name = self.cfg.get("run_name") or datetime.now().strftime("%Y%m%d-%H%M%S")
        self.output_dir = _ensure_directory(explicit_dir or (parent_dir / run_name))
        self.unidock_binary = self.cfg.get("unidock_binary", "unidock")
        self.ligand_collection: LigandCollection | None = None

    def _normalise_aliases(self) -> None:
        """Support legacy nt/ns keys by mapping them to modern flags."""
        alias_pairs = (("nt", "exhaustiveness"), ("ns", "max_step"))
        for alias, canonical in alias_pairs:
            if alias in self.cfg and canonical not in self.cfg:
                self.cfg[canonical] = self.cfg[alias]
            search_cfg = self.cfg.get("search")
            if isinstance(search_cfg, dict) and alias in search_cfg and canonical not in search_cfg:
                search_cfg[canonical] = search_cfg[alias]

    # ----------------------------- public API ---------------------------------

    def run(self, dry_run: bool = False) -> None:
        """Prepare inputs and either print or execute the Uni-Dock command."""
        self._log(f"Preparing Uni-Dock run ({self.mode} mode).")
        self.ligand_collection = self._prepare_ligands()
        cmd = self._build_command()
        command_txt = self.output_dir / "command.txt"
        command_txt.write_text(_shlex_join(cmd) + "\n")
        self._log(f"Wrote command file to {command_txt}")
        if dry_run:
            print("Dry run selected; not executing Uni-Dock.")
            print(_shlex_join(cmd))
            return
        self._log("Launching Uni-Dock ...")
        subprocess.run(cmd, check=True)

    # --------------------------- preparation logic ---------------------------

    def _prepare_ligands(self) -> LigandCollection:
        ligands: List[Path]
        if self.mode == "screening":
            ligands = self._prepare_screening_ligands()
        else:
            ligands = self._prepare_batch_ligands()
        if not ligands:
            raise ConfigError("No ligand files were discovered; check your configuration.")
        ligands = _unique(ligands)
        idx_path = self.output_dir / "ligands.index"
        idx_path.write_text("\n".join(str(path) for path in ligands) + "\n")
        self._log(f"Prepared {len(ligands)} ligand(s); index at {idx_path}.")
        return LigandCollection(paths=ligands, index_file=idx_path)

    def _prepare_batch_ligands(self) -> List[Path]:
        lig_cfg = self.cfg.get("ligands")
        if not lig_cfg:
            raise ConfigError("Batch mode requires a 'ligands' section in the config.")
        candidates: List[Path] = []
        # Explicit files listed.
        for entry in lig_cfg.get("paths", []):
            path = _as_path(self.config_dir, entry)
            if path is None:
                continue
            if path.is_file():
                candidates.append(path)
            else:
                raise ConfigError(f"Ligand file {path} does not exist.")
        # Directory scan.
        directory = lig_cfg.get("directory")
        if directory:
            pattern = lig_cfg.get("pattern", "*.pdbqt")
            dir_path = _as_path(self.config_dir, directory)
            if dir_path is None or not dir_path.is_dir():
                raise ConfigError(f"Ligand directory {directory} does not exist.")
            candidates.extend(sorted(dir_path.glob(pattern)))
        # Glob pattern(s)
        glob_patterns = lig_cfg.get("glob") or []
        if isinstance(glob_patterns, str):
            glob_patterns = [glob_patterns]
        for pattern in glob_patterns:
            pattern_path = self.config_dir / pattern
            matches = sorted(
                Path(match).resolve() for match in glob(str(pattern_path), recursive=True)
            )
            candidates.extend(matches)
        # ligand_index file to reuse
        index_file = lig_cfg.get("index_file")
        if index_file:
            idx_path = _as_path(self.config_dir, index_file)
            if idx_path is None or not idx_path.is_file():
                raise ConfigError(f"Ligand index file {index_file} not found.")
            for token in idx_path.read_text().split():
                candidates.append(Path(token))
        return [path for path in candidates if path.is_file()]

    def _prepare_screening_ligands(self) -> List[Path]:
        screening_cfg = self.cfg.get("screening")
        if not screening_cfg:
            raise ConfigError("Screening mode requires a 'screening' section.")
        target = screening_cfg.get("target")
        if not target:
            raise ConfigError("Screening configs must provide a target name.")
        dataset_dir = _as_path(self.config_dir, screening_cfg.get("dataset_dir")) or (
            self.config_dir / "indata"
        )
        _ensure_directory(dataset_dir)
        if not self.skip_extract:
            tarball = _as_path(self.config_dir, screening_cfg.get("tarball")) or (
                dataset_dir / f"{target}.tar.bz2"
            )
            if tarball.exists():
                self._extract_tarball(tarball, dataset_dir)
            else:
                self._log(f"No tarball found at {tarball}; assuming ligands already extracted.")
        patterns = screening_cfg.get("ligand_patterns")
        if not patterns:
            patterns = [
                "{target}_unique/*.pdbqt",
                "ligands_unique/*.pdbqt",
                "{target}_unique_charged/*.pdbqt",
            ]
        ligands: List[Path] = []
        for pattern in patterns:
            formatted = pattern.format(target=target)
            ligands.extend(sorted(dataset_dir.glob(formatted)))
        receptor = self.cfg.get("receptor")
        if not receptor:
            receptor_path = dataset_dir / f"{target}.pdbqt"
            if not receptor_path.exists():
                raise ConfigError(
                    f"Could not find receptor file at {receptor_path}. "
                    "Set 'receptor' explicitly in the config."
                )
            self.cfg["receptor"] = str(receptor_path)
        return ligands

    def _extract_tarball(self, tarball: Path, destination: Path) -> None:
        self._log(f"Extracting {tarball} into {destination} ...")
        try:
            with tarfile.open(tarball) as archive:
                extract_sig = inspect.signature(archive.extractall)
                kwargs = {}
                if "filter" in extract_sig.parameters:
                    kwargs["filter"] = "data"
                archive.extractall(destination, **kwargs)
        except tarfile.TarError as exc:
            raise ConfigError(f"Failed to extract {tarball}: {exc}") from exc

    # --------------------------- command generation ---------------------------

    def _build_command(self) -> List[str]:
        receptor = self.cfg.get("receptor")
        maps = self.cfg.get("maps")
        if receptor:
            receptor_path = _as_path(self.config_dir, receptor)
            if receptor_path is None or not receptor_path.is_file():
                raise ConfigError(f"Receptor file {receptor} not found.")
        elif maps:
            receptor_path = _as_path(self.config_dir, maps)
            if receptor_path is None:
                raise ConfigError("Could not resolve maps path.")
        else:
            raise ConfigError("Either 'receptor' or 'maps' must be provided.")

        box_cfg = self.cfg.get("box")
        if not box_cfg:
            raise ConfigError("Missing 'box' section describing docking box coordinates.")
        center = _vector3(box_cfg.get("center"), "box.center")
        size = _vector3(box_cfg.get("size"), "box.size")

        cmd: List[str] = [self.unidock_binary]
        if receptor:
            cmd.extend(["--receptor", str(receptor_path)])
        else:
            cmd.extend(["--maps", str(receptor_path)])

        assert self.ligand_collection is not None
        cmd.extend(["--ligand_index", str(self.ligand_collection.index_file)])
        cmd.extend(["--center_x", f"{center[0]:.4f}"])
        cmd.extend(["--center_y", f"{center[1]:.4f}"])
        cmd.extend(["--center_z", f"{center[2]:.4f}"])
        cmd.extend(["--size_x", f"{size[0]:.4f}"])
        cmd.extend(["--size_y", f"{size[1]:.4f}"])
        cmd.extend(["--size_z", f"{size[2]:.4f}"])

        cmd.extend(["--dir", str(self.output_dir)])

        flex = self.cfg.get("flex")
        if flex:
            flex_path = _as_path(self.config_dir, flex)
            if flex_path is None or not flex_path.exists():
                raise ConfigError(f"Flex file {flex} does not exist.")
            cmd.extend(["--flex", str(flex_path)])

        options_map = {
            "scoring": "--scoring",
            "search_mode": "--search_mode",
            "num_modes": "--num_modes",
            "exhaustiveness": "--exhaustiveness",
            "max_step": "--max_step",
            "max_evals": "--max_evals",
            "refine_step": "--refine_step",
            "max_gpu_memory": "--max_gpu_memory",
            "seed": "--seed",
            "cpu": "--cpu",
            "verbosity": "--verbosity",
            "min_rmsd": "--min_rmsd",
            "energy_range": "--energy_range",
            "spacing": "--spacing",
            "paired_batch_size": "--paired_batch_size",
        }
        for key, flag in options_map.items():
            value = self._get_option_value(key)
            if value is not None:
                cmd.extend([flag, str(value)])

        bool_flags = {
            "autobox": "--autobox",
        }
        for key, flag in bool_flags.items():
            value = self._get_option_value(key)
            if value:
                cmd.append(flag)

        extra_args = self.cfg.get("extra_args") or []
        if isinstance(extra_args, str):
            extra_args = extra_args.split()
        cmd.extend(str(arg) for arg in extra_args)

        # Retain the exact config used for the run.
        resolved_config = self.output_dir / "resolved_config.json"
        resolved_config.write_text(json.dumps(self.cfg, indent=2) + "\n")
        return cmd

    # ------------------------------ utilities --------------------------------

    def _log(self, message: str) -> None:
        if self.verbose:
            print(f"[universal-runner] {message}")

    def _get_option_value(self, key: str):
        if key in self.cfg:
            return self.cfg[key]
        search_cfg = self.cfg.get("search")
        if isinstance(search_cfg, dict) and key in search_cfg:
            return search_cfg[key]
        return None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Config-driven Uni-Dock runner.")
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="Path to the JSON config that describes the run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare everything but skip launching Uni-Dock.",
    )
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip tarball extraction for screening configs.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show progress messages.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Config file {args.config} does not exist.", file=sys.stderr)
        return 2
    try:
        runner = UniDockRunner(
            config_path=config_path,
            skip_extract=args.skip_extract,
            verbose=args.verbose,
        )
        runner.run(dry_run=args.dry_run)
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 3
    except subprocess.CalledProcessError as exc:
        print(f"Uni-Dock exited with code {exc.returncode}", file=sys.stderr)
        return exc.returncode or 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
