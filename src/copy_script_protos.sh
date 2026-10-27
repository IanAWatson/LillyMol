#!/usr/bin/env bash

# Install standalone-script protos without requiring LillyMol Python bindings.
# Run `bazel build //lillymol_protos:generated_protos` first.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd -P)
destdir=${1:-${here}/../lib}
package_dir="${destdir}/lillymol_protos"
python_protos=(xgboost_model_pb2.py random_forest_model_pb2.py class_label_translation_pb2.py linear_scaling_pb2.py)
ruby_protos=(xgboost_model_pb.rb random_forest_model_pb.rb)
sources=("${here}/bazel-bin/lillymol_protos/__init__.py" "${here}/lillymol_protos/check_runtime.py")
for binding in "${python_protos[@]}"; do
  sources+=("${here}/bazel-bin/lillymol_protos/${binding}")
done
for binding in "${ruby_protos[@]}"; do
  sources+=("${here}/bazel-bin/xgboost/${binding}")
done

# Check the complete set before changing the installation.
for source_file in "${sources[@]}"; do
  if [[ ! -s "${source_file}" ]]; then
    echo "Missing ${source_file}; build //lillymol_protos:generated_protos first" >&2
    exit 1
  fi
done
mkdir -p "${package_dir}"
for source_file in "${sources[@]}"; do
  destination="${package_dir}/$(basename "${source_file}")"
  if ! cmp -s "${source_file}" "${destination}"; then
    echo "Copying lillymol_protos/$(basename "${source_file}")"
    install -m 644 "${source_file}" "${destination}"
  fi
done
