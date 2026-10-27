#!/usr/bin/env bash

# Stage Bazel-generated bindings beside the executable-test runner.
# Run `bazel build //test:generated_protos` first.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)

# Validate both inputs before changing either installed binding.
for binding in lillymol_tests_pb2.py lillymol_tests_pb.rb; do
  if [[ ! -s "${here}/bazel-bin/test/${binding}" ]]; then
    echo "Missing bazel-bin/test/${binding}; build //test:generated_protos first" >&2
    exit 1
  fi
done

for binding in lillymol_tests_pb2.py lillymol_tests_pb.rb; do
  source_file="${here}/bazel-bin/test/${binding}"
  destination="${here}/test/${binding}"
  # Timestamps can be misleading after copying files between machines.
  if ! cmp -s "${source_file}" "${destination}"; then
    echo "Copying test/${binding}"
    install -m 644 "${source_file}" "${destination}"
  fi
done
