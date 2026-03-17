# Plan: Implement Neovim-style `deep_extend` for Config Merging

## Problem

The `consolidate_conf` method in `db.py` uses shallow merge. When an experiment config defines `pp: {dorado: {model: "hac"}}`, the entire template's `pp` dict is replaced, losing keys like `bin_path`, `device`, `lsf_queue`, etc.

This forces experiments to redundantly specify `pp.dorado.template` to reload defaults.

## Solution

Implement `deep_extend(behavior, base, *overrides)` similar to Neovim's `vim.tbl_deep_extend()`.

## Changes

### 1. Add `deep_extend` function to `db.py` (after imports, ~line 12)

```python
from copy import deepcopy
from typing import Literal

def deep_extend(
    behavior: Literal["force", "keep", "error"],
    base: dict,
    *overrides: dict
) -> dict:
    """
    Recursively merge dicts, similar to vim.tbl_deep_extend().

    - "force": override wins on conflict
    - "keep": base wins on conflict
    - "error": raise on conflict
    - Lists are REPLACED, not appended
    - Only recurses when BOTH values are dicts
    """
    if not isinstance(base, dict):
        raise TypeError(f"base must be dict, got {type(base).__name__}")

    result = deepcopy(base)

    for override in overrides:
        if override is None:
            continue
        if not isinstance(override, dict):
            raise TypeError(f"override must be dict, got {type(override).__name__}")
        result = _deep_extend_pair(behavior, result, override)

    return result


def _deep_extend_pair(behavior: str, base: dict, override: dict) -> dict:
    """Internal: merge two dicts recursively. Modifies base in-place."""
    for key, override_val in override.items():
        if key not in base:
            base[key] = deepcopy(override_val)
        elif isinstance(base[key], dict) and isinstance(override_val, dict):
            _deep_extend_pair(behavior, base[key], override_val)
        else:
            if behavior == "force":
                base[key] = deepcopy(override_val)
            elif behavior == "keep":
                pass
            elif behavior == "error":
                raise ValueError(f"Conflict for '{key}': {base[key]!r} vs {override_val!r}")
    return base
```

### 2. Modify `consolidate_conf` method (~lines 739-744)

Replace the shallow assignment loop with deep merge for specific keys:

```python
# Keys that should be deep-merged rather than replaced
DEEP_MERGE_KEYS = {'pp', 'fracture', 'extension_config'}

for k, v in xp_template.items():
    if k not in vars(self) and k not in self.special_vars:
        setattr(self, k, v)
    elif k in DEEP_MERGE_KEYS:
        # Deep merge: template as base, experiment as override
        template_val = v
        experiment_val = getattr(self, k, {})
        if isinstance(template_val, dict) and isinstance(experiment_val, dict):
            merged = deep_extend("force", template_val, experiment_val)
            setattr(self, k, merged)
            logger.debug(f'deep merged {k}')
        else:
            logger.debug(f'kept {k} from experiment (type mismatch)')
    else:
        logger.debug(f'kept {k} from experiment conf')
```

### 3. Simplify `run_dorado` (~lines 136-153)

Remove redundant template re-loading since `consolidate_conf` now handles deep merge:

```python
# Before (remove this block):
template_path = dorado_conf['template'].replace('${prefix}', xp.prefix)
template_data = yaml.load(open(template_path), Loader=yaml.FullLoader)
if 'pp' in template_data and 'dorado' in template_data['pp']:
    dorado_template = template_data['pp']['dorado'].copy()
# ... merge loop ...

# After (simplified):
dorado_template = pp['dorado'].copy()  # Already deep-merged by consolidate_conf

# Keep **args runtime override:
if args:
    for k, v in args.items():
        dorado_template[k] = v
```

## Files to Modify

| File | Lines | Change |
|------|-------|--------|
| `/home/projects/nyosef/pedro/src/ogtk/ogtk/utils/db.py` | ~12 | Add `deep_extend` function |
| `/home/projects/nyosef/pedro/src/ogtk/ogtk/utils/db.py` | 739-744 | Use deep merge in `consolidate_conf` |
| `/home/projects/nyosef/pedro/src/ogtk/ogtk/utils/db.py` | 136-153 | Simplify `run_dorado` template loading |

## Backwards Compatibility

- Existing configs work unchanged
- If experiment fully specifies all keys, result is identical to shallow merge
- Direct key access (`xp.fracture['start_k']`) still works

## Verification

1. Run existing experiment config:
   ```bash
   cd ~/projects/lt
   python -c "
   from ogtk.utils.db import Xp
   xp = Xp('conf/20260126_pet_groups_11_12_ont.yml')
   print('lsf_queue:', xp.pp['dorado'].get('lsf_queue'))
   print('bin_path:', xp.pp['dorado'].get('bin_path'))
   "
   ```
   Should show `gsla_high_gpu` and the bin_path from template.

2. Test deep_extend directly:
   ```python
   from ogtk.utils.db import deep_extend
   base = {'pp': {'dorado': {'bin': '/path', 'model': 'sup'}}}
   override = {'pp': {'dorado': {'model': 'hac'}}}
   result = deep_extend("force", base, override)
   assert result['pp']['dorado'] == {'bin': '/path', 'model': 'hac'}
   ```

3. Run the dorado step to ensure it still submits to correct queue.
