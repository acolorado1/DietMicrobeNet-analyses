#!/usr/bin/env bash
# Run DM[H]N on every IBD pair directory, recording runtime and peak memory.
# Run from the DietMicrobeNet run-workflow directory.
set -euo pipefail

IBD_DIR="/Users/burkhang/Desktop/ibd_dmn"
DMN_DIR="/Users/burkhang/Code_Projs/DietMicrobeNet"

python "${IBD_DIR}/run_workflow_benchmark.py" \
    --directories "${IBD_DIR}"/participants/*/hg_mgx \
                  "${IBD_DIR}"/participants/*/hg_mtx \
                  "${IBD_DIR}"/participants/*/htx_mgx \
                  "${IBD_DIR}"/participants/*/htx_mtx \
    --script "${DMN_DIR}/run_workflow.py" \
    --output "${IBD_DIR}/benchmark_results.csv" \
    --extra-args \
        --dietmicrobehost \
        --e-weights \
        --n-weights \
        --include-orgs \
        --abundance-col "Abundance_RPKs" \
        --all-food
