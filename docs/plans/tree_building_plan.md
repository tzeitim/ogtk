# Consolidate Tree Generation into cassiopeia_petracer Extension

**Status**: Implemented (Feb 2026, branch `feat/cassiopeape_ext_tree`)

## Context

Tree building previously lived in standalone scripts (`build_tree_sc.py`, `build_tree_single.py`) with hardcoded paths, manually submitted via LSF. The cassiopeia_petracer extension stopped at allele tables — everything after was ad-hoc. This plan consolidated tree generation as a proper extension step following the dorado pattern (prepare + submit + detect).

This is a prerequisite for the smolscells extension plan (`docs/smolscells_extension_plan.md`), which needs structured tree outputs to aggregate into a LineageForest.

---

## Usage

### CLI (standalone / LSF jobs)

```bash
# Full run (SC mode, all cells)
python -m ogtk.ltr.fracture.extensions.cassiopeia_petracer build-tree \
  --input /path/to/alleles_pl_collapsed.parquet \
  --outdir /path/to/trees \
  --mode sc --solver nj

# Quick test run (subsample)
python -m ogtk.ltr.fracture.extensions.cassiopeia_petracer build-tree \
  --input /path/to/alleles_pl_collapsed.parquet \
  --outdir /tmp/tree_test \
  --mode sc --solver nj --skip-branch-lengths \
  --test-mode --top-x-cells 200 --sample-n-cells 50 \
  --top-x-intbcs 10 --sample-n-intbcs 3

# Single-molecule mode (one intBC/sbc group)
python -m ogtk.ltr.fracture.extensions.cassiopeia_petracer build-tree \
  --input /path/to/alleles_pl.parquet \
  --outdir /path/to/trees \
  --mode sm --intbc TCTGAA --sbc SBC1 --solver nj
```

#### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--mode` | required | `sc` (one tree, all intBCs as characters) or `sm` (one tree per intBC/sbc) |
| `--solver` | `nj` | `nj`, `vanilla`, `mcgs` |
| `--spanning-deletions` | `unedited` | `unedited` (state 0), `missing` (null), `use` (keep) |
| `--skip-branch-lengths` | false | Skip convexml branch length estimation |
| `--min-cells` | 10 | Skip if fewer cells/molecules |
| `--test-mode` | false | Enable subsampling |
| `--top-x-cells` / `--sample-n-cells` | None | Two-step cell subsampling: rank by size, take top X, sample N |
| `--top-x-intbcs` / `--sample-n-intbcs` | None | Same for intBCs |
| `--allele-rep-thresh` | 1.0 | Cassiopeia allele representation threshold |
| `--min-branch` | 0.01 | Minimum branch length for convexml |
| `--sample-name` | `sample` | Label used in plot filenames |
| `--n-rcols` | 15 | Number of r columns (r1..rN) |

### Extension pipeline (YAML config)

Add `build_trees` to extension steps:

```yaml
extensions:
  - cassiopeia_petracer

extension_config:
  cassiopeia_petracer:
    # ... existing allele params ...
    solver: nj
    min_cells_for_tree: 10
    spanning_deletions: unedited
    skip_branch_lengths: false
    min_branch_length: 0.01
    tree_mode: auto                # auto: sc if collapse_to_cells else sm
    allele_rep_thresh: 1.0

    # Test mode (set tree_test_mode: false for production)
    tree_test_mode: true
    top_x_cells: 200
    sample_n_cells: 50
    top_x_intbcs: 10
    sample_n_intbcs: 3

    # LSF (set tree_use_lsf: false for inline execution)
    tree_use_lsf: false
    # tree_lsf_queue: gsla-cpu
    # tree_lsf_cores: 1
    # tree_lsf_mem: 4G
    # tree_conda_env: cas11

extension_steps:
  cassiopeia_petracer:
    - parse_contigs
    - classify_cassettes
    - plug_cassiopeia
    - build_trees              # tree generation step
```

With `tree_mode: auto`, SC mode is selected when `collapse_to_cells` is true, SM mode otherwise.

---

## Output directory structure

```
{workdir}/trees/span_{spanning_deletions}/
  # SC mode:
  ├── newick.txt
  ├── character_matrix.parquet
  ├── tdata.h5td
  ├── tdata_md5.txt
  ├── sc_{n_cells}_{sample}_{solver}_span-{mode}_circ.png
  ├── sc_{n_cells}_{sample}_{solver}_span-{mode}_linear.png
  └── done

  # SM mode (one subdir per group):
  ├── {intbc}_{sbc}/
  │   ├── newick.txt
  │   ├── character_matrix.parquet
  │   ├── tdata.h5td
  │   ├── tdata_md5.txt
  │   ├── sm_{n_mols}_{intbc}_{sbc}_{solver}_circ.png
  │   ├── sm_{n_mols}_{intbc}_{sbc}_{solver}_linear.png
  │   └── done
  └── ...
```

Re-running skips groups that already have a `done` file (cache). Delete the `done` file or use `--force` (CLI) to recompute.

---

## Implementation details

**File**: `ogtk/ltr/fracture/extensions/cassiopeia_petracer.py`

### Enum

`CassiopeiaStep.BUILD_TREES = 'build_trees'`

### Config fields on `CassiopeiaConfig`

```python
# Tree generation
solver: str = 'nj'                          # nj | vanilla | mcgs
min_cells_for_tree: int = 10
skip_branch_lengths: bool = False
min_branch_length: float = 0.01
tree_mode: str = 'auto'                     # auto | sc | sm
allele_rep_thresh: float = 1.0

# Test mode — two-step subsampling: rank by size, take top X, then sample N
tree_test_mode: bool = False
top_x_cells: Optional[int] = None
sample_n_cells: Optional[int] = None
top_x_intbcs: Optional[int] = None
sample_n_intbcs: Optional[int] = None

# LSF submission (follows dorado pattern)
tree_use_lsf: bool = True
tree_lsf_queue: str = 'gsla-cpu'
tree_lsf_cores: int = 1
tree_lsf_mem: str = '4G'
tree_wait_for_completion: bool = False
tree_conda_env: Optional[str] = None
```

### Methods added to `CassiopeiaLineageExtension`

| Method | Description |
|---|---|
| `_build_trees()` | Orchestrator: loads alleles, applies test subsampling, dispatches to SC/SM or LSF |
| `_build_tree_sc(allele_table, rcols, outdir)` | SC mode: one tree per sample, all intBCs as characters. Handles spanning deletions, convexml branch lengths, wide allele pivot, TreeData, plots |
| `_build_tree_sm(allele_table, rcols, outdir, intbc, sbc)` | SM mode: one tree per (intBC, sbc) group |
| `_get_solver(solver_name)` | Returns configured Cassiopeia solver (NJ uses `weighted_hamming_distance`) |
| `_estimate_branch_lengths(cas_tree, config)` | `reconstruct_ancestral_characters()` then `convexml` branch length estimation |
| `_make_tree_cmd(job, allele_path, rcols)` | Builds the CLI command string for LSF submission |
| `_submit_tree_jobs_lsf(jobs, allele_path, rcols)` | Submits bsub jobs via stdin (follows tree_qc.py pattern) |

### Module-level functions (used by both extension methods and CLI)

| Function | Description |
|---|---|
| `pivot_alleles_wide(allele_table, rcols)` | Pivots from long (cell x intBC) to wide (cell) with `intBC_rN` column naming. Handles multiple UMIs via mode. |
| `add_depth(tdata, tree_key)` | Edge-counting depth via `pycea.pp.add_depth` plus normalized depth |
| `add_weighted_depth(tdata, tree_key)` | Branch-length-weighted depth via `nx.single_source_dijkstra_path_length` |
| `write_done(outdir, reason)` | Writes `done` marker (or `done_{reason}` for early exits) |
| `_build_tree_standalone(args)` | Entry point for CLI / LSF jobs, mirrors the tested `build_tree_sc.py` / `build_tree_single.py` logic |

### Wiring in `process()`

`BUILD_TREES` runs after `SEGMENTED_ALLELE` (the last allele-table step):

```python
if self.should_run_step(CassiopeiaStep.BUILD_TREES.value):
    result = self._build_trees()
    final_results.update(result.results)
    final_metrics.update(result.metrics)
```

### Dead code removed

`_pycea_explore()` — superseded by `_build_tree_sc` / `_build_tree_sm`.

---

## Integration with existing scripts

After consolidation, `build_tree_sc.py` and `build_tree_single.py` (in `~/projects/lt/scripts/`) are optional. The primary path is through the extension step or the CLI entry point. The standalone scripts can be retired once the extension is stable.
