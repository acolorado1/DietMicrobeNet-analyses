#!/usr/bin/env bash
# Run DM[H]N on every IBD pair directory, recording runtime and peak memory.
# Run from the DietMicrobeNet run-workflow directory.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # DietMicrobeNet-analyses
DMN_DIR="/Users/burkhang/Code_Projs/DietMicrobeNet"                 # DietMicrobeNet repository
IBD_DIR="${REPO_DIR}/data/IBD"

mkdir -p "${REPO_DIR}/data/benchmarking"

python "${REPO_DIR}/scripts/preprocessing/run_workflow_benchmark.py" \
    --directories "${IBD_DIR}"/participants/*/hg_mgx \
                  "${IBD_DIR}"/participants/*/hg_mtx \
                  "${IBD_DIR}"/participants/*/htx_mgx \
                  "${IBD_DIR}"/participants/*/htx_mtx \
    --script "${DMN_DIR}/run_workflow.py" \
    --output "${REPO_DIR}/data/benchmarking/benchmark_IBD.csv" \
    --extra-args \
        --dietmicrobehost \
        --e-weights \
        --n-weights \
        --include-orgs \
        --abundance-col "Abundance_RPKs" \
        --all-food
