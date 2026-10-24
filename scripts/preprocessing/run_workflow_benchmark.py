#!/usr/bin/env python3
"""
Runs run_workflow.py once per directory (per participant), benchmarking
runtime and peak memory usage for each run, and writes results to a CSV.

Usage:
    python run_workflow_benchmark.py \
        --directories "/path/to/ZIM089/" "/path/to/ZIM091/" ... \
        --script run_workflow.py \
        --output benchmark_results.csv \
        --extra-args --dietmicrobehost --e-weights --n-weights --include-orgs \
                      --abundance-col "Abundance_RPKs" --all-food
"""

import argparse
import csv
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import psutil


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark run_workflow.py across multiple directories, one run per directory."
    )
    parser.add_argument(
        "--directories",
        nargs="+",
        required=True,
        help="List of directories to process, one run_workflow.py call per directory.",
    )
    parser.add_argument(
        "--script",
        default="run_workflow.py",
        help="Path to run_workflow.py (default: run_workflow.py, assumed on PATH or in cwd).",
    )
    parser.add_argument(
        "--output",
        default="benchmark_results.csv",
        help="Path to write the benchmark CSV (default: benchmark_results.csv).",
    )
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=0.2,
        help="Seconds between memory polls (default: 0.2).",
    )
    # Everything after -- gets passed straight through to run_workflow.py
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Remaining flags/args to pass through to run_workflow.py "
             "(e.g. --dietmicrobehost --e-weights ...). Put this LAST.",
    )
    return parser.parse_args()


def monitor_memory(proc, poll_interval, peak_mem_holder, stop_event):
    """
    Polls RSS memory of proc and all its children, tracking the peak
    combined memory usage (in MB) until stop_event is set or the process ends.
    """
    peak_bytes = 0
    try:
        ps_proc = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        peak_mem_holder["peak_mb"] = 0.0
        return

    while not stop_event.is_set():
        try:
            if not ps_proc.is_running():
                break
            mem = ps_proc.memory_info().rss
            # Include memory from any child processes spawned by the script
            for child in ps_proc.children(recursive=True):
                try:
                    mem += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            peak_bytes = max(peak_bytes, mem)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            break
        time.sleep(poll_interval)

    peak_mem_holder["peak_mb"] = peak_bytes / (1024 * 1024)


def run_and_benchmark(script, directory, extra_args, poll_interval):
    """
    Runs: python <script> --directories <directory> <extra_args...>
    Returns a dict of benchmark stats for this single run.
    """
    cmd = [
        sys.executable,
        script,
        "--directories",
        directory,
        *extra_args,
    ]

    start_time = time.time()
    start_ts = datetime.now().isoformat(timespec="seconds")

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    peak_mem_holder = {"peak_mb": 0.0}
    stop_event = threading.Event()
    monitor_thread = threading.Thread(
        target=monitor_memory,
        args=(proc, poll_interval, peak_mem_holder, stop_event),
        daemon=True,
    )
    monitor_thread.start()

    stdout, stderr = proc.communicate()
    stop_event.set()
    monitor_thread.join(timeout=poll_interval * 2)

    end_time = time.time()
    runtime_sec = end_time - start_time

    return {
        "directory": directory,
        "start_time": start_ts,
        "runtime_sec": round(runtime_sec, 3),
        "peak_memory_mb": round(peak_mem_holder["peak_mb"], 2),
        "return_code": proc.returncode,
        "success": proc.returncode == 0,
        "stderr_tail": stderr[-500:].replace("\n", " ") if stderr else "",
    }


def main():
    args = parse_args()

    script_path = Path(args.script)
    if not script_path.exists():
        print(f"WARNING: script path '{args.script}' not found relative to cwd; "
              f"will still attempt to run it (may be on PATH).", file=sys.stderr)

    results = []
    total_dirs = len(args.directories)

    for i, directory in enumerate(args.directories, start=1):
        print(f"[{i}/{total_dirs}] Running on: {directory}")
        result = run_and_benchmark(
            script=args.script,
            directory=directory,
            extra_args=args.extra_args,
            poll_interval=args.poll_interval,
        )
        results.append(result)

        status = "OK" if result["success"] else f"FAILED (code {result['return_code']})"
        print(
            f"    -> {status} | runtime: {result['runtime_sec']}s | "
            f"peak memory: {result['peak_memory_mb']} MB"
        )
        if not result["success"] and result["stderr_tail"]:
            print(f"    stderr (tail): {result['stderr_tail']}")

    # Write results to CSV
    fieldnames = [
        "directory",
        "start_time",
        "runtime_sec",
        "peak_memory_mb",
        "return_code",
        "success",
        "stderr_tail",
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nBenchmark complete. Results written to: {args.output}")

    n_success = sum(r["success"] for r in results)
    n_failed = total_dirs - n_success
    total_runtime = sum(r["runtime_sec"] for r in results)
    print(f"Summary: {n_success} succeeded, {n_failed} failed, "
          f"total wall time across runs: {round(total_runtime, 2)}s")


if __name__ == "__main__":
    main()