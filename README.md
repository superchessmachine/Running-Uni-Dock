# Uni-Dock Example Inputs

This repo keeps the two public Uni-Dock workflows in a stripped form: every receptor/ligand file is a dummy placeholder so you can drop in your own assets before launching a run.

## Paired batch example (`example/paired_batch/`)

- `example/paired_batch/paired_batch_config.json` lists three dummy pairs (`example_pair_a/b/c`) so you can see how to scale the schema.
- `example/paired_batch/ligand_config.json` sets a neutral docking box; edit it once you know your binding site.
- Replace the placeholder PDBQT files (`receptor_example_*.pdbqt`, `ligand_example_*.pdbqt`) with your prepared receptors/ligands, then extend or rename the JSON entries as needed.
- Launch from your Uni-Dock build with:

```bash
build/unidock --paired_batch_size <batch_size> \
  --ligand_index example/paired_batch/paired_batch_config.json \
  --center_x <x> --center_y <y> --center_z <z> \
  --size_x <sx> --size_y <sy> --size_z <sz> \
  --dir <output_dir> --exhaustiveness <value> --max_step <value>
```

## Screening example (`example/screening_test/`)

- `run_dock.py` is the stock helper script from DP Tech.
- `config.json` names the placeholder target `example_target`; every file in `indata/` shares that prefix (`example_target.pdbqt`, `example_target.tar.bz2`, and the ligands under `example_target_unique/`).
- Replace the receptor, tarball, and ligand files with your own, update `config.json`, and re-run the script each time you change targets or docking settings.
- Execute directly inside `example/screening_test/`:

```bash
cd example/screening_test
python run_dock.py
```

The script extracts the tarball, writes `<target>_ligands.txt`, runs Uni-Dock, and collects scores in `result/`.
