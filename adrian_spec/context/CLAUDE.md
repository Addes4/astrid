# CLAUDE.md

ASTRID annotates cell types in scRNA-seq / Xenium data, designed to hold up on tumour
samples. Four chained stages: clustering → annotation → validation → chromosomal damage.

## Environment

- **Python lives in `.venv310` (3.10).** Use `.venv310/bin/python`, never the system
  interpreter — the pinned stack (numpy 1.23.4, scipy 1.8.1, scanpy 1.9.3) does not
  install on newer Pythons.
- **`Rscript` must be on PATH.** The annotation stage shells out to R and needs
  SingleR, SingleCellExperiment, tidyverse, Matrix.
- No `requirements.txt` at the repo root yet; the pinned list is in the README.
- There are no tests. See "Verifying a change" below.

## Running it

```bash
.venv310/bin/python ASTRID_v0.01.py --all \
  --input_file  <sample>_preASTRID.h5ad \
  --input_prefix <sample> \
  --output_file <sample>_outASTRID.h5ad \
  --output_clustering_results <sample>_ASTRID_Result.csv \
  --out_dir <dir> \
  --author_type "<obs column with author labels>"
```

- **Always pass `--out_dir`.** Without it `ASTRID_v0.01.py:647` concatenates
  `os.getcwd()` and the prefix with no separator and writes to a sibling of the repo.
- The script `os.chdir`s to its own directory at import (`:15`), so `data/…` paths
  resolve from anywhere — but relative paths passed on the command line are then
  resolved against the repo, not your shell's cwd. Prefer absolute paths.
- Stages can run separately (`--clustering`, `--annotation`, `--validation`,
  `--damage`) but they chain through files derived from `--output_file`, so a later
  stage needs the earlier stage's outputs present under matching names.

## Mental model

**ASTRID never types individual cells.** It clusters recursively, aggregates each
cluster to pseudobulk, and labels the *cluster*; the label is then broadcast back to
cells. One row in the result CSV is one cluster. Weight anything you compute by
`CellCount`.

Two independent methods assign the type and are allowed to disagree: SingleR against a
reference atlas, and hand-written boolean marker formulas in
`data/ExpectedCellTypesMarkers.csv` (`*` = AND, `+` = OR, `!` = NOT). The disagreement
is a feature — `OtherCellTypesPassed` is the most informative column in the output.

Python and R communicate through CSV files on disk, not memory. The pseudobulk matrix
and the SingleR output are the interface, not debug artefacts.

## Config that silently changes results

| File | Effect |
|---|---|
| `data/excludeSingleR.csv` | Drops reference samples before SingleR. Currently removes the whole `Neutrophils`, `Eosinophils`, `pre-adipocytes` and `Chondrocytes` labels (28 labels → 24). |
| `data/ExpectedCellTypesMarkers.csv` | The marker formulas. Several are permissive pure-OR expressions. |
| `data/CellTypeAliases.json` | Lineage groupings used for the CD4 / PTPRC special cases. |
| `data/IntOGen-DriverGenes_CSCC.tsv` | Driver genes for the damage stage — hardcoded to skin cancer (`ASTRID_v0.01.py:576`). |

## Traps

Things that look like results but are not:

- **`cancer_percentage` is not a prediction.** It is the author's own label matched
  against the string `"Cancer"` (`:399`). No author cancer labels → the column is all
  zero by construction.
- **`CNVCorrelation` and `newCNVScore*` are relative.** The correlation is against the
  mean profile of the top 5 % most damaged cells *in the same sample* (`:537–543`).
  Compare clusters within one run; never across runs, never against a fixed threshold.
- **`FormulaPassed` can be the string `"not applicable"`** (`:460`), so the column is
  object-typed. Compare with `== True`.
- **The `unassigned` pseudo-cluster is skipped by validation** (`:446`), so its marker
  columns are `NaN` by design. Check how big it is every run.
- **`clustering_level_11` is a cap, not convergence** (`clustering_depth = 12`, `:110`).

Known bugs, unfixed at time of writing:

- `ASTRID_Tools.py:186` builds the k-NN adjacency with `np.repeat(indices[:, 0], k)`
  as the row index instead of `np.repeat(np.arange(n), k)`. `indices[:, 0]` is an
  arbitrary neighbour, not the cell itself. Isolated nodes become singleton clusters,
  which are marked `unassigned`.
- `RunSingleR.R:12` reads the pseudobulk CSV with R's default `check.names = TRUE`,
  mangling gene symbols (`HLA-DRA` → `HLA.DRA`). ~4000 genes, including all of HLA
  class I and II, then fail to intersect with the reference and are silently dropped.
- `ASTRID_Tools.py:293` calls R with `os.system` and ignores the exit code, then reads
  the output CSV. A failed R run either raises far from the cause or silently reuses a
  stale CSV from a previous run.
- `RunSingleR.R` falls back to the default reference without warning if `--reference`
  points at a path that does not exist.

## Scaling limits

`normalize_data` densifies with `.toarray()` (`ASTRID_Tools.py:227`) and the clustering
recipe builds a full n×n distance matrix (`:278`). Roughly 1.8 GB at 15k cells and
9.8 GB at 35k, and it applies at level 1 where the sample is one cluster.

## Conventions

- Large outputs (`*_outASTRID.h5ad` can be several GB) do not belong in git.
- `ASTRID_Vignette.ipynb` is tracked **with its stored outputs**, which are the
  reference numbers for a known-good p1 run. Running it dirties the file — work on a
  copy outside version control instead.

## Verifying a change

There is no test suite. The practical check is to re-run a known sample and compare:

- composition weighted by `CellCount`
- the NMI lines printed by each stage
- the formula pass rate and the size of the `unassigned` cluster

The vignette's stored outputs give a reference point for p1: NMI 0.753 (level 1) →
0.955 (final), SingleR 0.680, 208 clusters, 101/208 formulas passed.

## Repo status

This is a fork (`Addes4/astrid`, forked from `jhausserlab/astrid`) used for an
ongoing research project, not the upstream repo. Implications:

- Known bugs listed above are upstream's, not introduced here — don't "fix" them
  in passing unless the task is specifically to fix that bug. If a fix seems
  warranted, flag it rather than changing behavior silently; it may need to be
  reported upstream rather than patched locally.
- Changes here are meant to stay local experimentation (e.g. the RMSD-vs-corr
  work in `RunSingleR_RMSD.R` / `compare_corr_vs_rmsd.py`) unless told otherwise —
  there's no upstream PR workflow in place yet.
- Divergence from upstream is not automatically tracked; if a change touches a
  file also present in `jhausserlab/astrid`, mention that it diverges from
  upstream so it isn't mistaken for shared behavior.