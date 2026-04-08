# Plan: Add SLURM support to fracture pipeline distributed compute

## Context

The fracture pipeline (`ogtk.ltr.fracture.pipeline`) currently only supports LSF for distributed compute. LSF is not actually used inside the `pipeline/` module itself — the distributed-compute entry points live in `ogtk/utils/db.py` (dorado basecalling) and `ogtk/ltr/fracture/post/tree_qc.py` (tree QC batch jobs), both called from pipeline steps. All scheduler calls (`bsub`, `bjobs`, job-id regex, status codes) are hardcoded with no abstraction layer. We want to run the same pipeline on SLURM clusters without forking code paths, and keep commits atomic on the `slurm` branch (no co-author trailers).

## Approach

Introduce a minimal scheduler abstraction in a new module `ogtk/utils/scheduler.py`, refactor the existing LSF code paths to use it, then add a SLURM implementation. Backend is selected via a new `scheduler` config key (values: `lsf`, `slurm`), with `use_lsf: true` preserved as a deprecated alias that maps to `scheduler: lsf` so existing YAMLs keep working.

### Scheduler interface

`ogtk/utils/scheduler.py` (new):

- `@dataclass JobSpec` — fields: `command: str | list[str]`, `job_name: str`, `stdout: str`, `stderr: str`, `queue: str | None`, `gpu: str | None`, `memory: str | None`, `cores: int | None`, `extra: list[str]`, `stdin_script: str | None`.
- `class Scheduler(ABC)`:
  - `submit(spec: JobSpec) -> str` — returns job id
  - `poll(job_ids: list[str]) -> dict[str, str]` — returns normalized statuses: `PENDING|RUNNING|DONE|FAILED|UNKNOWN`
  - `wait(job_ids, poll_interval, max_wait_hours) -> int` — shared concrete implementation in base class built on `poll()`; returns 0/1. Lifts the existing logic from `_monitor_lsf_jobs` (db.py:604-659).
  - `monitor_hint(job_ids) -> str` — e.g. `bjobs ...` / `squeue -j ...`
- `class LSFScheduler(Scheduler)`:
  - `submit`: builds `bsub -q … -gpu … -R rusage[mem=…] -o … -e … -J …` from `JobSpec`; parses `Job <12345>`.
  - `poll`: `bjobs -noheader -o "jobid stat" <ids>`; maps `PEND→PENDING`, `RUN→RUNNING`, `DONE→DONE`, `EXIT→FAILED`, missing→`DONE` (matches current behavior at db.py:635-644).
- `class SLURMScheduler(Scheduler)`:
  - `submit`: builds `sbatch --partition=<queue> --gres=gpu:<n> --mem=<mem> --cpus-per-task=<cores> -J <name> -o <out> -e <err> [--time=…] <script>`; parses `Submitted batch job (\d+)`.
    - GPU translation: if `JobSpec.gpu` is set, parse LSF-style `num=N:gmem=XG` → `--gres=gpu:N` (+ warn that `gmem` is ignored) OR accept SLURM-native `slurm_gres` override from the template. Prefer an explicit `slurm_gres` key when present.
    - When a command is passed (not a script path), submit via `sbatch --wrap="<cmd>"`.
  - `poll`: `squeue -h -j <ids> -o "%i %T"` + fallback `sacct -j <id> -n -o State` for jobs that already left the queue; map `PENDING/RUNNING/COMPLETED/FAILED/CANCELLED/TIMEOUT/NODE_FAIL` to normalized states.
- `def make_scheduler(config: dict) -> Scheduler`:
  - Read `config.get('scheduler')`; if missing, fall back to `'lsf' if config.get('use_lsf') else None`.
  - Emit a `DeprecationWarning` via logger when `use_lsf` is used without `scheduler`.

### Refactor existing LSF code to use the interface

**`ogtk/utils/db.py`**
- `_submit_dorado_lsf_jobs` (line 450) → rename to `_submit_dorado_jobs(xp, commands, dorado_template, scheduler)`. Build a `JobSpec` per command using existing fields (`lsf_queue`, `lsf_gpu`, `lsf_mem`) plus new optional SLURM-specific keys (`slurm_partition`, `slurm_gres`, `slurm_mem`, `slurm_time`, `slurm_cpus`). `JobSpec.queue/gpu/memory` are filled from whichever set matches the active backend, so the scheduler-specific knobs stay in the YAML.
- `_submit_iterative_dorado_lsf` (line 526) → `_submit_iterative_dorado(...)` with same refactor; still writes the script via `_generate_iterative_dorado_script` and submits its path.
- `_monitor_lsf_jobs` (line 604) → delete; callers use `scheduler.wait(...)`.
- Call sites at lines 236 and 250: replace `if dorado_template.get('use_lsf', False):` with `scheduler = make_scheduler(dorado_template); if scheduler is not None:`.
- Log hints at lines 518 and 595 become `scheduler.monitor_hint(...)`.
- Config keys list at line 156 (`boolean_keys`) stays; add `'scheduler'` to the string-keys handling alongside existing template keys (verify where string keys are whitelisted).

**`ogtk/ltr/fracture/post/tree_qc.py`**
- `batch_tree_qc` (line ~354): add `scheduler: str | None = None` parameter; keep `use_lsf`/`lsf_queue`/`lsf_memory` for back-compat; add `slurm_partition`, `slurm_memory`, `slurm_cpus`. Build a dict and call `make_scheduler`. The bsub block at lines 462-471 and job-id parse at 476 become `scheduler.submit(JobSpec(...))`.
- `_wait_for_lsf_jobs` (line 534) → delete; use `scheduler.wait(...)`.

### Config / template changes

- `templates/dorado_template.yaml`: add commented-out example block:
  ```yaml
  scheduler: lsf   # or: slurm
  # SLURM-only keys (ignored by LSF):
  # slurm_partition: gpu
  # slurm_gres: "gpu:1"
  # slurm_mem: 64G
  # slurm_cpus: 4
  # slurm_time: "24:00:00"
  ```
  Leave `use_lsf`/`lsf_*` keys in place and documented as the LSF-specific names.
- No change needed to `examples/multi_flowcell_example.yaml` unless we want a SLURM example — out of scope for this PR.

### Files touched

- `ogtk/utils/scheduler.py` (new)
- `ogtk/utils/db.py` (refactor lines ~156, 236, 250, 450-660)
- `ogtk/ltr/fracture/post/tree_qc.py` (refactor lines ~354-492, 534-560)
- `templates/dorado_template.yaml` (doc/keys)
- `tests/` — add a unit test for `scheduler.py` that monkeypatches `subprocess.run` and verifies (a) LSF `bsub` command construction and job-id parsing, (b) SLURM `sbatch` command construction and `Submitted batch job` parsing, (c) `make_scheduler` back-compat for `use_lsf: true`.

### Atomic commits (branch: `slurm`, no co-author trailer)

1. **Add `ogtk/utils/scheduler.py` with `Scheduler` base + `LSFScheduler`** — pure new code, no behavior change. Includes unit tests for LSF path.
2. **Refactor `db.py` dorado submission to use `LSFScheduler`** — behavior-preserving; `use_lsf: true` still works unchanged. Delete `_monitor_lsf_jobs`.
3. **Refactor `tree_qc.py` batch_tree_qc to use `LSFScheduler`** — behavior-preserving. Delete `_wait_for_lsf_jobs`.
4. **Add `SLURMScheduler` + `make_scheduler` backend selection** — new SLURM class, `scheduler:` config key, `use_lsf` deprecation warning, tests for SLURM path.
5. **Docs: update `dorado_template.yaml` with `scheduler` key and SLURM example.**

Each commit should build and leave tests green on its own.

## Verification

- **Unit tests** (`pytest tests/test_scheduler.py`): mock `subprocess.run`, assert the exact argv for LSF and SLURM submission from a shared `JobSpec`; assert job-id parsing for both; assert `make_scheduler({'use_lsf': True})` returns `LSFScheduler` with a deprecation warning; assert `make_scheduler({'scheduler': 'slurm'})` returns `SLURMScheduler`.
- **LSF regression**: on the LSF cluster, run the existing dorado pipeline end-to-end against a small sample with the unchanged `use_lsf: true` config and confirm identical behavior (job submission, polling, exit path).
- **SLURM smoke test**: on a SLURM cluster (or a login node with `sbatch`/`squeue` available), run the dorado step with `scheduler: slurm` + `slurm_partition`/`slurm_gres` set against a tiny input; confirm `sbatch` is invoked, job id is captured, `squeue` polling transitions to `COMPLETED`, and downstream pipeline steps resume.
- **tree_qc smoke test**: run `batch_tree_qc(scheduler='slurm', ...)` on a couple of tree dirs and confirm result collection works after `wait=True`.
