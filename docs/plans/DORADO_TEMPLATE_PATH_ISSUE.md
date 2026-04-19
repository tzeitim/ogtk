# Dorado Template Path Issue - Summary

## Issue Description

Pipeline fails with the following error when running the dorado step:

```
FileNotFoundError: [Errno 2] No such file or directory: "f'{self.prefix}/src/ogtk/templates/dorado_template.yaml'"
```

The error occurs at `ogtk/utils/db.py:138` in the `run_dorado` function.

## Root Cause

The issue occurs because of a mismatch between the template path format in the configuration files and what the code expects:

1. **XP Template Files** (e.g., `ontCarlin_xp_template.yml`) use f-string syntax:
   ```yaml
   pp:
     dorado:
       template: "f'{self.prefix}/src/ogtk/templates/dorado_template.yaml'"
   ```

2. **Code Expectation** (`db.py:137`) expects `${prefix}` placeholder format:
   ```python
   template_path = dorado_conf['template'].replace('${prefix}', xp.prefix)
   ```

3. **The Problem**: The `pp` (preprocessing) section is NOT a "special pattern" (like `xp_`, `pro_`, `sample_`), so values in it are assigned directly without evaluating f-strings. This means the literal string `"f'{self.prefix}/..."` is passed to the code, which then tries to replace `${prefix}` but finds nothing.

## Why It Worked Before

The Nov 4, 2025 successful run used config `/home/projects/nyosef/pedro/projects/lt/conf/20251030_pet_pilot_ont.yml`, which contained an **override** in the experiment config itself:

```yaml
pp:
  dorado:
    template: "${prefix}/projects/lt/conf/ont_pet_10mer_xp_template.yml"
```

This override used the correct `${prefix}` format, bypassing the problematic f-string in the xp_template.

## Why It's Failing Now

The current config `/home/projects/nyosef/pedro/projects/lt/conf/20251201_carlin_4-5-13-44_ont.yml` is **missing** the `pp:` section override, so it's using the unevaluated f-string from `ontCarlin_xp_template.yml`.

## Solutions

### Solution 1: Add Override in Experiment Config (Recommended)

Add this to your experiment config file (e.g., `20251201_carlin_4-5-13-44_ont.yml`):

```yaml
pp:
  dorado:
    template: "${prefix}/src/ogtk/templates/dorado_template.yaml"
```

### Solution 2: Fix XP Template Files

Change all xp_template files to use `${prefix}` format instead of f-string format:

**Before:**
```yaml
template: "f'{self.prefix}/src/ogtk/templates/dorado_template.yaml'"
```

**After:**
```yaml
template: "${prefix}/src/ogtk/templates/dorado_template.yaml"
```

Files to update:
- `/home/projects/nyosef/pedro/projects/lt/conf/ontCarlin_xp_template.yml:25`
- `/home/projects/nyosef/pedro/projects/lt/conf/ont_pet_10mer_xp_template.yml:42`

### Solution 3: Update Code to Handle F-strings

Modify `ogtk/utils/db.py:137` to handle both formats:

```python
# Load template and merge configurations
template_path_str = dorado_conf['template']
# Handle both ${prefix} placeholder format and f-string format
if 'self.' in template_path_str:
    # F-string format: evaluate and replace self. with xp.
    template_path = eval(template_path_str.replace('self.', 'xp.'))
else:
    # Placeholder format: replace ${prefix}
    template_path = template_path_str.replace('${prefix}', xp.prefix)
template_data = yaml.load(open(template_path), Loader=yaml.FullLoader)
```

This mirrors the approach already used for `temp_symlink_dir` at `db.py:424-426`.

## Recommendation

**Use Solution 1** (add override in experiment config) as it's the quickest fix and follows the pattern already established in working configs. For long-term consistency, also consider Solution 2 to fix the xp_template files.

## Related Files

- **Code**: `ogtk/utils/db.py:137`
- **Working config**: `projects/lt/conf/20251030_pet_pilot_ont.yml`
- **Failing config**: `projects/lt/conf/20251201_carlin_4-5-13-44_ont.yml`
- **XP templates**:
  - `projects/lt/conf/ontCarlin_xp_template.yml`
  - `projects/lt/conf/ont_pet_10mer_xp_template.yml`

## Timeline

- **Nov 4, 2025**: Dorado step worked successfully with `20251030_pet_pilot_ont.yml`
- **Dec 2, 2025**: Dorado step fails with `20251201_carlin_4-5-13-44_ont.yml`
- **Key difference**: Missing `pp:` override in the failing config
