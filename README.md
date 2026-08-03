# GateMeClass Module

## What this module does

Runs the GateMeClass R wrapper for cell-type annotation.

- Wrapper: `gatemeclass_wrapper.R`
- Local runner: `run_gatemeclass.sh`
- Output: `data_import_predicted_labels.tar.gz` by default

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

Run against specific benchmark outputs:

```bash
GATEMECLASS_TRAIN_MATRIX=/path/to/data_import.train.matrix.tar.gz \
GATEMECLASS_TRAIN_LABELS=/path/to/data_import.train.labels.tar.gz \
GATEMECLASS_TEST_MATRIX=/path/to/data_import.test.matrices.tar.gz \
GATEMECLASS_METADATA=/path/to/data_import.metadata.json.gz \
GATEMECLASS_OUTPUT_DIR=/tmp/gatemeclass-check \
GATEMECLASS_EXCLUDED_DATASETS=FR-FCM-Z3YR,FlowCyt \
bash models/gatemeclass/run_gatemeclass.sh
```

The local runner now matches the benchmark defaults more closely:

- `GATEMECLASS_SAMPLING=0.1`
- `GATEMECLASS_K=20`
- `GATEMECLASS_KNN_BACKEND=caret`
- `GATEMECLASS_KNN_QUERY_CHUNK_SIZE=50000`
- `GATEMECLASS_CORES=1`
- `GATEMECLASS_BLAS_THREADS=1`
- `GATEMECLASS_NAME=data_import`

Those conservative defaults are deliberate. GateMeClass can spike RSS when it
annotates multiple samples in parallel, so raising worker count should be
treated as a memory tradeoff rather than a free speedup.

The wrapper supports three KNN backends:

- `caret` uses the compatibility implementation with matrix-based fitting and
  direct final-model prediction.
- `class` removes caret's model wrapper but retains brute-force neighbor search.
- `kmknn` uses the exact indexed `BiocNeighbors` KMKNN search and bounded query
  chunks. It is normally faster, but caret's fuzzy boundary-distance ties and
  randomized class-vote ties can produce different predictions on tied data.

Backend selection is explicit so accepted benchmark runs cannot silently change
KNN semantics. The wrapper only requests labels from the vendored annotation
function; callers of `GateMeClass_annotate()` retain full cell signatures by
default.

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

The input paths, output directory, run name, and model parameters can all be
overridden with `GATEMECLASS_*` environment variables shown above.
