# Running Uni-Dock Scripts

This repo now ships a *universal runner* that hides the boilerplate required to
launch Uni-Dock in both classic batch mode (one receptor, many ligands) and the
screening workflow that DP-Tech shared in the upstream `screening_test`
directory.  Everything is driven from a single JSON config so you can prepare a
run once and relaunch it with one command.

## Folder layout

- `universal/run_unidock.py` – CLI entry point. Reads configs, discovers ligand
  `.pdbqt` files, writes the ligand index that Uni-Dock expects, and finally
  launches `unidock`.
- `universal/configs/*.json` – ready-to-edit examples. `batch_example.json` uses
  the small receptor/ligand set that lives in `paired_batch/`, and
  `screening_example.json` reproduces the official screening test input in
  `screening_test/`.
- `universal/runs/` – created automatically. Each run gets its own subfolder
  that stores the ligand index, resolved config, and the exact Uni-Dock command.

The legacy helper data from the upstream repo is untouched under
`paired_batch/` and `screening_test/`, so you can keep using the official
examples or swap in your own receptor/ligand files.

## Quick start

```bash
# Inspect the batch example without launching Uni-Dock
python universal/run_unidock.py \
  --config universal/configs/batch_example.json \
  --dry-run --verbose

# Run the screening test (extracts def.tar.bz2 the first time) 
python universal/run_unidock.py \
  --config universal/configs/screening_example.json \
  --verbose
```

Command-line options:

- `--dry-run` – prepare ligand indices/output folders but skip the actual
  `unidock` execution. Great for reviewing the generated command line.
- `--skip-extract` – in screening mode, reuse the already-extracted ligand set
  instead of unpacking the `*.tar.*` archive.
- `-v/--verbose` – print progress messages.

## Config reference

Every config is regular JSON. Paths are resolved relative to the config file,
which keeps the setup portable.

- `mode`: `"batch"` (default) or `"screening"`.
- `unidock_binary`: defaults to `unidock`. Point it at `build/unidock` if you
  are running from a build tree.
- `run_name`: optional name that becomes part of the `universal/runs/<name>`
  folder. A timestamp is used when omitted.
- `receptor` or `maps`: path to the receptor PDBQT or AD4 map prefix.
- `box.center` and `box.size`: either arrays `[x, y, z]` or `{ "x": ..., ... }`
  objects. These drive `--center_*` and `--size_*`.
- `ligands`: **batch mode only**. Provide one or more of:
  - `directory` + optional `pattern` (defaults to `*.pdbqt`)
  - `glob`: glob string or list of patterns
  - `paths`: explicit file list
  - `index_file`: reuse a pre-built ligand index file
  The runner deduplicates everything and writes a fresh `ligands.index` file for
  Uni-Dock.
- `screening`: **screening mode only**. Fields:
  - `target`: matches the upstream `def`, `mmp13`, etc.
  - `dataset_dir`: folder containing `<target>.tar.bz2` and `<target>.pdbqt`.
  - `tarball`: optional override when the archive lives elsewhere.
  - `ligand_patterns`: override the default `["{target}_unique/*.pdbqt", ...]`.
- `search`: optional sub-dict for clean configs. Values bubble up to the main
  command, e.g. `{"search_mode": "balance", "exhaustiveness": 384}`.
- Top-level scalar keys mirror Uni-Dock flags – `scoring`, `num_modes`,
  `exhaustiveness`, `max_step`, `refine_step`, `seed`, `max_gpu_memory`,
  `verbosity`, etc. Legacy `nt`/`ns` keys are automatically mapped to
  `exhaustiveness`/`max_step`.
- `output`: optional `{ "dir": "...", "parent_dir": "..." }` section that
  controls where the run folder is created. Omit it to get
  `universal/runs/<timestamp>`.
- `extra_args`: string or list of additional CLI tokens appended verbatim.
- `autobox`: boolean that adds `--autobox` to the command.

Each run folder captures:

1. `ligands.index` – the file passed to `--ligand_index`.
2. `command.txt` – the exact shell command.
3. `resolved_config.json` – the config after relative paths are resolved and
   screening defaults applied.

## Next steps

Copy one of the example configs in `universal/configs/`, point it to your
receptor and ligand directories, tweak the box/search settings, and run the
script. Because everything funnels through the same tool you can keep adding
new configs (e.g., different pockets, parameter sweeps) and launch them all with
the same `python universal/run_unidock.py --config ...` command.
