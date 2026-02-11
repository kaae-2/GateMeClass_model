#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
model_name="gatemeclass"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install conda to run GateMeClass." >&2
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "${model_name}"; then
  echo "ERROR: ${model_name} conda env not found. Create it with:" >&2
  echo "  conda env create -f ${script_dir}/../benchmark/envs/gateme_env.yml -n ${model_name}" >&2
  exit 1
fi

train_matrix="${script_dir}/out/data/data_preprocessing/default/data_import.train.matrix.tar.gz"
train_labels="${script_dir}/out/data/data_preprocessing/default/data_import.train.labels.tar.gz"
test_matrix="${script_dir}/out/data/data_preprocessing/default/data_import.test.matrices.tar.gz"

for required_file in "$train_matrix" "$train_labels" "$test_matrix"; do
  if [ ! -f "$required_file" ]; then
    echo "ERROR: missing input file: ${required_file}" >&2
    exit 1
  fi
done

output_dir="${script_dir}/out/data/analysis/default/${model_name}"
excluded_datasets="${GATEMECLASS_EXCLUDED_DATASETS:-}"

cmd=(
  conda run --no-capture-output -n "${model_name}" Rscript "${script_dir}/gatemeclass_wrapper.R"
  --name "${model_name}"
  --output_dir "${output_dir}"
  --data.train_matrix "${train_matrix}"
  --data.train_labels "${train_labels}"
  --data.test_matrix "${test_matrix}"
  --GMM_parameterization "V"
  --sampling "1.0"
  --k "20"
)

if [ -n "${excluded_datasets}" ]; then
  cmd+=(--excluded-datasets "${excluded_datasets}")
fi

"${cmd[@]}"
