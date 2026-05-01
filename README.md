# GateMeClass Module

## What this module does

Runs the GateMeClass R wrapper for cell-type annotation.

- Wrapper: `gatemeclass_wrapper.R`
- Local runner: `run_gatemeclass.sh`
- Output: `gatemeclass_predicted_labels.tar.gz`

The wrapper reads train/test tar archives, predicts labels per sample,
and maps labels back to benchmark ids.

## Run locally

```bash
bash models/gatemeclass/run_gatemeclass.sh
```

Optional exclusion list:

```bash
GATEMECLASS_EXCLUDED_DATASETS=FR-FCM-Z3YR,flowcyt bash models/gatemeclass/run_gatemeclass.sh
```

The local runner now matches the benchmark defaults more closely:

- `GATEMECLASS_SAMPLING=0.1`
- `GATEMECLASS_K=20`
- `GATEMECLASS_CORES=1`
- `GATEMECLASS_BLAS_THREADS=1`

Those conservative defaults are deliberate. GateMeClass can spike RSS when it
annotates multiple samples in parallel, so raising worker count should be
treated as a memory tradeoff rather than a free speedup.

Dataset exclusions are resolved from the passed label-key/order metadata first,
then from path fallbacks. Matching is case-insensitive and ignores separators,
so `FlowCyt`, `flowcyt`, and `flow-cyt` all match the same exclusion entry.

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
