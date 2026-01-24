#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install conda to run GateMeClass." >&2
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "gatemeclass"; then
  echo "ERROR: gatemeclass conda env not found. Create it with:" >&2
  echo "  conda env create -f ${script_dir}/../benchmark/envs/gateme_env.yml -n gatemeclass" >&2
  exit 1
fi

conda run -n gatemeclass Rscript "${script_dir}/gatemeclass_wrapper.R" \
  --name "gatemeclass" \
  --output_dir "${script_dir}/out/data/analysis/default/gatemeclass" \
  --data.train_matrix "${script_dir}/out/data/data_preprocessing/default/data_import.train.matrix.tar.gz" \
  --data.train_labels "${script_dir}/out/data/data_preprocessing/default/data_import.train.labels.tar.gz" \
  --data.test_matrix "${script_dir}/out/data/data_preprocessing/default/data_import.test.matrices.tar.gz" \
  --GMM_parameterization "V" \
  --sampling "1.0" \
  --k "20"
