# GateMeClass Module

## What this module does

Runs the GateMeClass R wrapper for cell-type annotation.

- Wrapper: `gatemeclass_wrapper.R`
- Local runner: `run_gatemeclass.sh`
- Output: `gatemeclass_predicted_labels.tar.gz`

The wrapper reads train/test tar archives, applies marker transformation,
predicts labels per sample, and maps labels back to benchmark ids.

## Run locally

```bash
bash models/gatemeclass/run_gatemeclass.sh
```

Optional exclusion list:

```bash
GATEMECLASS_EXCLUDED_DATASETS=FR-FCM-Z3YR,flowcyt bash models/gatemeclass/run_gatemeclass.sh
```

## Run as part of benchmark

Configured in `benchmark/Clustering_conda.yml` analysis stage; run with:

```bash
just benchmark
```

## What `run_gatemeclass.sh` needs

- Conda installed and env `gatemeclass` available
- Preprocessing outputs at `models/gatemeclass/out/data/data_preprocessing/default`
- `Rscript` and GateMeClass dependencies in the conda env
- Writable output directory `models/gatemeclass/out/data/analysis/default/gatemeclass`
