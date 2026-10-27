# Shared environment for the XGBoost and random-forest wrappers. Source this file.
# These scripts use the caller's model-building Python environment, without
# loading LillyMol's native extension or enforcing its Python ABI.
model_bin_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P) || return 1
if [[ ! -v LILLYMOL_HOME ]]; then
  export LILLYMOL_HOME=$(cd -- "${model_bin_dir}/../.." && pwd -P)
fi
export PATH="${model_bin_dir}:${PATH}"
export PYTHONPATH="${LILLYMOL_HOME}/lib:${model_bin_dir}${PYTHONPATH:+:${PYTHONPATH}}"

check_model_protos() {
  if [[ ! -s "${LILLYMOL_HOME}/lib/lillymol_protos/check_runtime.py" ]]; then
    echo "Model protos are not installed. Build //lillymol_protos:generated_protos and run src/copy_script_protos.sh." >&2
    return 1
  fi
  "${PYTHON:-python}" -m lillymol_protos.check_runtime
}
