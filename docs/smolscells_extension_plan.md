# smolscells Extension Plan

## Problem

The fracture pipeline produces solved trees through the `build_trees` step of the `cassiopeia_petracer` extension (see `docs/tree_building_plan.md`):
- **Single-cell (SC)**: one tree per sample using all intBCs as characters
- **Single-molecule (SM)**: one tree per (intBC, sbc) group

Tree outputs live under `{workdir}/trees/span_{mode}/` with `newick.txt`, `character_matrix.parquet`, `tdata.h5td`, and a `done` marker per tree.

There's no structured way to aggregate these fractions into a single object for comparison. The smolscells `LineageForest` is the natural container (`'sc'` + `'sm_*'` CassiopeiaTree slots) but lacks import methods for cassiopeia pipeline outputs.

## Goal

1. Add import methods to `LineageForest` (smolscells side) for loading solved trees
2. Add tree merging methods (combine SM trees from multiple sbcs)
3. Create a fracture extension (`smolscells`) that orchestrates the aggregation
4. Track pipeline parameters in `lf.uns` for comparing runs with different settings

---

## Part 1: smolscells — New methods on `LineageForest`

**File**: `~/src/smolscells/smolscells/lineage_forest.py`

### 1a. `from_tree_outputs()` — classmethod

Scans directories for solved tree outputs (from `build_tree_sc.py` / `build_tree_single.py`) and loads them into a LineageForest.

```python
@classmethod
def from_tree_outputs(
    cls,
    sc_dir: Path | str | None = None,
    sm_dirs: dict[str, Path | str] | None = None,
) -> "LineageForest":
    """Create LineageForest from pre-solved tree outputs.

    Parameters
    ----------
    sc_dir
        Directory containing SC tree outputs:
        - newick.txt
        - character_matrix.parquet
        Typically: {outdir}/{sample}/span_{mode}/
    sm_dirs
        Dict mapping intBC identifiers to directories, each containing:
        - newick.txt
        - character_matrix.parquet
        Typically: {outdir}/{tube}/{cell_line}/{tube}_{intbc}_{sbc}/
    """
```

**Implementation notes**:
- Read `newick.txt` → pass to `CassiopeiaTree(tree=newick_str)`
- Read `character_matrix.parquet` → pandas DataFrame → pass to `CassiopeiaTree(character_matrix=df)`
- SC goes to `lf.add_sc_tree(tree)`
- Each SM goes to `lf.add_sm_tree_with_id(intbc_key, tree)`
- Check for `done` file — skip dirs without it (failed/incomplete runs)
- Store source paths in `lf.uns['source_paths']`

**Expected directory layout for SC** (from `cassiopeia_petracer` `build_trees` step, SC mode):
```
{workdir}/trees/span_{spanning_deletions}/
├── newick.txt
├── character_matrix.parquet
├── tdata.h5td
├── tdata_md5.txt
└── done
```

**Expected directory layout for SM** (from `cassiopeia_petracer` `build_trees` step, SM mode):
```
{workdir}/trees/span_{spanning_deletions}/{intbc}_{sbc}/
├── newick.txt
├── character_matrix.parquet
├── tdata.h5td
├── tdata_md5.txt
└── done
```

### 1b. `merge_sm_trees()` — instance method

Concatenates character matrices from multiple SM trees into a single tree slot. Use case: multiple sbcs from the same biological sample should become a single tree.

```python
def merge_sm_trees(
    self,
    source_keys: list[str],
    target_key: str,
    remove_sources: bool = False,
) -> CassiopeiaTree:
    """Merge multiple SM trees by concatenating their character matrices.

    The resulting CassiopeiaTree is UNSOLVED (no topology, just character matrix).
    The user must solve it after merging.

    Parameters
    ----------
    source_keys
        Tree keys to merge (e.g., ['sm_0', 'sm_1', 'sm_2'])
    target_key
        Key for the merged tree (e.g., 'sm_merged' or 'sc')
    remove_sources
        If True, delete source trees after merge
    """
```

**Implementation notes**:
- Get `character_matrix` from each source tree
- `pd.concat([cm1, cm2, ...], axis=0)` — rows are cells/molecules, columns are characters
- Handle column mismatch: if SM trees have different character columns (different intBCs), use outer join and fill missing with `-1` (missing state)
- Create new `CassiopeiaTree(character_matrix=merged_cm)`
- Store merge info in `lf.uns['merged_trees'][target_key] = {'source_keys': [...], 'n_cells_per_source': {...}}`
- If `remove_sources`, delete source trees from `self.trees` and invalidate their tdata cache

### 1c. `store_allele_table()` — instance method

Store a raw allele table in `uns` for reference. Useful for re-running cassiopeia conversion with different parameters.

```python
def store_allele_table(
    self,
    allele_table,  # pl.DataFrame or pd.DataFrame
    key: str = 'allele_table',
):
    """Store allele table in uns for reference."""
    import pandas as pd
    if hasattr(allele_table, 'to_pandas'):
        allele_table = allele_table.to_pandas()
    self.uns[key] = allele_table
```

---

## Part 2: ogtk fracture extension — `SmolscellsExtension`

### 2a. New file

**File**: `~/src/ogtk/ogtk/ltr/fracture/extensions/smolscells_ext.py`

```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Type, Set
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults


@dataclass
class SmolscellsConfig(ExtensionConfig):
    # Where to find solved tree outputs
    sc_tree_dir: str | None = None        # Path to SC tree dir (newick.txt + character_matrix.parquet)
    sm_tree_dir: str | None = None        # Parent dir with per-(intBC,sbc) subdirs
    output_name: str = 'lineage_forest'   # Output filename stem


class SmolscellsExtension(PostProcessorExtension):

    @property
    def name(self) -> str:
        return "smolscells"

    def get_config_class(self) -> Type[ExtensionConfig]:
        return SmolscellsConfig

    def process(self, contigs_path: Path) -> StepResults:
        from smolscells import LineageForest

        workdir = contigs_path.parent
        config = self.config
        metrics = {}

        # 1. DISCOVER tree directories
        sc_dir = Path(config.sc_tree_dir) if config.sc_tree_dir else None
        sm_dirs = self._discover_sm_dirs(config.sm_tree_dir) if config.sm_tree_dir else None

        # 2. BUILD LineageForest
        lf = LineageForest.from_tree_outputs(sc_dir=sc_dir, sm_dirs=sm_dirs)

        # 3. RECORD pipeline parameters from cassiopeia extension
        cas_config = self.xp.extension_config.get('cassiopeia_petracer', {})
        lf.uns['pipeline_params'] = {
            'spanning_deletions': cas_config.get('spanning_deletions', 'use'),
            'min_umis_per_cell': cas_config.get('min_umis_per_cell', 4),
            'min_umi_agreement': cas_config.get('min_umi_agreement', 0.5),
            'allele_rep_thresh': cas_config.get('allele_rep_thresh', 1.0),
            'collapse_to_cells': cas_config.get('collapse_to_cells', True),
            'top_n_cells': cas_config.get('top_n_cells', None),
        }

        metrics['n_trees'] = lf.n_trees
        metrics['tree_keys'] = lf.tree_keys

        # 4. SAVE
        output_path = workdir / f'{config.output_name}.pkl'
        lf.save(output_path)

        return StepResults(
            results={'lineage_forest_path': str(output_path)},
            metrics=metrics,
        )

    def _discover_sm_dirs(self, parent_dir: str) -> dict[str, Path]:
        """Scan parent_dir for subdirectories containing solved trees."""
        sm_dirs = {}
        parent = Path(parent_dir)
        for d in sorted(parent.rglob('done')):
            tree_dir = d.parent
            newick = tree_dir / 'newick.txt'
            if newick.exists():
                # Use directory name as key
                sm_dirs[tree_dir.name] = tree_dir
        return sm_dirs
```

### 2b. Register extension

**File**: `~/src/ogtk/ogtk/ltr/fracture/extensions/__init__.py`

Add:
```python
from . import smolscells_ext      # Auto-register smolscells extension
```

---

## Part 3: Pipeline parameter tracking

The extension records cassiopeia parameters in `lf.uns['pipeline_params']` so users can compare across runs:

```python
lf1 = LineageForest.load('run_span_unedited/lineage_forest.pkl')
lf2 = LineageForest.load('run_span_missing/lineage_forest.pkl')

print(lf1.uns['pipeline_params']['spanning_deletions'])  # 'unedited'
print(lf2.uns['pipeline_params']['spanning_deletions'])  # 'missing'
```

The tracked parameters:
- `spanning_deletions`: unedited / missing / use
- `min_umis_per_cell`: UMI count threshold for cell collapse
- `min_umi_agreement`: consensus threshold
- `allele_rep_thresh`: cassiopeia character matrix threshold
- `collapse_to_cells`: whether UMI collapse was applied
- `top_n_cells`: cell count cap

---

## Files summary

| File | Action |
|------|--------|
| `~/src/smolscells/smolscells/lineage_forest.py` | Add `from_tree_outputs()`, `merge_sm_trees()`, `store_allele_table()` |
| `~/src/ogtk/ogtk/ltr/fracture/extensions/smolscells_ext.py` | **Create** |
| `~/src/ogtk/ogtk/ltr/fracture/extensions/__init__.py` | Add import line |

## Verification

1. Manually create a LineageForest from an existing tree output dir (e.g. the 4T1 SC tree):
   ```python
   from smolscells import LineageForest
   lf = LineageForest.from_tree_outputs(sc_dir='path/to/span_unedited/')
   print(lf)  # Should show 1 tree ('sc')
   lf.csc     # Should return CassiopeiaTree
   ```
2. Test `merge_sm_trees()` with two mock SM trees
3. Run the extension via fracture pipeline config on a sample workdir
