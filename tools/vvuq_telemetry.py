#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
# Copyright (C) 2026 p-j-r
"""Quiet-runner best-of-3 wall-time, peak RSS, and ASME V&V 20 accuracy |E| telemetry ledger.

Invokes release validation and benchmark test suites across `trials` repetitions (default 3),
records minimum (best-of-N) and median wall-clock duration, exact per-process peak RSS
(`os.wait4` `ru_maxrss` on Linux/macOS; `K32GetProcessMemoryInfo` on Windows), and parses
`[VV20 ...]` / `[VV20-EP ...]` comparison errors |E| from test output without imposing flaky
wall-time hard-fails.
"""

from __future__ import annotations

import argparse
import ctypes
import datetime
import json
import os
import pathlib
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Any


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

DEFAULT_SUITES: list[dict[str, Any]] = [
    {
        "name": "validation_thermal_nafems",
        "package": "nyx-engine",
        "test": "validation_thermal_nafems",
        "args": [],
    },
    {
        "name": "validation_photonics",
        "package": "nyx-engine",
        "test": "validation_photonics",
        "args": [],
    },
    {
        "name": "validation_mems_piezo",
        "package": "nyx-engine",
        "test": "validation_mems_piezo",
        "args": [],
    },
    {
        "name": "validation_cqed_cross_code",
        "package": "nyx-engine",
        "test": "validation_cqed_cross_code",
        "args": [],
    },
    {
        "name": "validation_cqed_experiment",
        "package": "nyx-engine",
        "test": "validation_cqed_experiment",
        "args": [],
    },
    {
        "name": "validation_cqed_tls_scqubits",
        "package": "nyx-engine",
        "test": "validation_cqed_tls_scqubits",
        "args": [],
    },
    {
        "name": "examples_quantitative",
        "package": "nyx-engine",
        "test": "examples",
        "args": [],
    },
    {
        "name": "r22_frozen_fmsh_oracle",
        "package": "nyx-ems",
        "test": "golden",
        "args": ["golden_soi_wire_and_hybrid_frozen_solver_oracles_r22"],
    },
]

RE_VV20 = re.compile(
    r"\[VV20\s+(PASS|FAIL)\s+([^\]]+)\]\s+[^:]+:\s+"
    r"S=([^\s]+)\s+.*?D=([^\s]+)\s+.*?E=([^\s]+)\s+\(([^%]+)%\).*?"
    r"u_num=([^\s,]+),\s+u_inp=([^\s,]+),\s+u_D=([^\s,]+),\s+u_val=([^\s,]+),.*?"
    r"tol=([^\s,]+)"
)

RE_VV20_EP = re.compile(
    r"\[VV20-EP\s+(PASS|FAIL)\s+([^\]]+)\]\s+[^:]+:\s+"
    r"\[S_low,\s*S_up\]=\[([^,]+),\s*([^\]]+)\].*?"
    r"S_mid=([^\s,]+),\s+D=([^\s,]+),\s+delta_ep=([^\s,]+),\s+u_ale=([^\s,]+),\s+"
    r"\|E_mid\|=([^\s]+)\s+\(([^%]+)%\),\s+tol=([^\s,]+)"
)


def _win32_peak_memory(handle: int) -> tuple[int | None, int | None]:
    """Reads (PeakWorkingSetSize, PeakPagefileUsage) in bytes via kernel32.K32GetProcessMemoryInfo."""

    class ProcessMemoryCounters(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_uint32),
            ("PageFaultCount", ctypes.c_uint32),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    counters = ProcessMemoryCounters()
    counters.cb = ctypes.sizeof(ProcessMemoryCounters)
    kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
    ok = kernel32.K32GetProcessMemoryInfo(
        ctypes.c_void_p(handle), ctypes.byref(counters), counters.cb
    )
    if not ok:
        return None, None
    return int(counters.PeakWorkingSetSize), int(counters.PeakPagefileUsage)


def _git_head() -> str | None:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short=12", "HEAD"],
            cwd=REPO_ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return out or None
    except Exception:
        return None


def _suite_cargo_cmd(suite: dict[str, Any], no_run: bool = False) -> list[str]:
    cmd = [
        "cargo",
        "test",
        "--release",
        "-p",
        suite["package"],
        "--test",
        suite["test"],
    ]
    if no_run:
        cmd.append("--no-run")
    else:
        cmd.append("--")
        cmd.extend(suite["args"])
        cmd.append("--nocapture")
    return cmd


def _parse_vv20_output(output: str) -> list[dict[str, Any]]:
    targets: list[dict[str, Any]] = []
    for line in output.splitlines():
        m = RE_VV20.search(line)
        if m:
            status, tid, s_val, d_val, e_val, rel_pct, u_num, u_inp, u_d, u_val, tol = m.groups()
            try:
                e_f = float(e_val)
                targets.append(
                    {
                        "id": tid,
                        "kind": "vv20",
                        "passed": status == "PASS",
                        "s": float(s_val),
                        "d": float(d_val),
                        "e": e_f,
                        "abs_e": abs(e_f),
                        "rel_e": float(rel_pct) / 100.0,
                        "u_num": float(u_num),
                        "u_input": float(u_inp),
                        "u_d": float(u_d),
                        "u_val": float(u_val),
                        "tol": float(tol),
                    }
                )
            except ValueError:
                pass
            continue
        m_ep = RE_VV20_EP.search(line)
        if m_ep:
            status, tid, s_lo, s_hi, s_mid, d_val, d_ep, u_ale, abs_e, rel_pct, tol = m_ep.groups()
            try:
                targets.append(
                    {
                        "id": tid,
                        "kind": "vv20_epistemic",
                        "passed": status == "PASS",
                        "s_low": float(s_lo),
                        "s_up": float(s_hi),
                        "s": float(s_mid),
                        "d": float(d_val),
                        "abs_e": float(abs_e),
                        "rel_e": float(rel_pct) / 100.0,
                        "delta_ep": float(d_ep),
                        "u_aleatoric": float(u_ale),
                        "tol": float(tol),
                    }
                )
            except ValueError:
                pass
    return targets


def _run_once(cmd: list[str], env: dict[str, str]) -> dict[str, Any]:
    with tempfile.TemporaryFile(mode="w+b") as log_f:
        t0 = time.perf_counter()
        proc = subprocess.Popen(
            cmd,
            cwd=REPO_ROOT,
            env=env,
            stdout=log_f,
            stderr=subprocess.STDOUT,
        )
        peak_rss_bytes: int | None = None
        peak_commit_bytes: int | None = None
        if sys.platform == "win32":
            rc = proc.wait()
            peak_rss_bytes, peak_commit_bytes = _win32_peak_memory(int(proc._handle))
        else:
            _, status, rusage = os.wait4(proc.pid, 0)
            rc = os.waitstatus_to_exitcode(status)
            if sys.platform == "darwin":
                peak_rss_bytes = int(rusage.ru_maxrss)
            else:
                peak_rss_bytes = int(rusage.ru_maxrss) * 1024
        wall_s = time.perf_counter() - t0
        log_f.seek(0)
        output = log_f.read().decode("utf-8", errors="replace")
    return {
        "returncode": rc,
        "wall_s": wall_s,
        "peak_rss_bytes": peak_rss_bytes,
        "peak_commit_bytes": peak_commit_bytes,
        "output": output,
    }


def _find_test_executable(package: str, test_name: str, env: dict[str, str]) -> list[str] | None:
    """Finds the compiled test binary via `cargo test --message-format=json` so RSS/wall-time
    measure the solver process itself rather than the `cargo` wrapper."""
    cmd = [
        "cargo",
        "test",
        "--release",
        "--no-run",
        "--message-format=json",
        "-p",
        package,
        "--test",
        test_name,
    ]
    res = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        return None
    for line in res.stdout.splitlines():
        if not line.startswith("{"):
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue
        if (
            msg.get("reason") == "compiler-artifact"
            and msg.get("target", {}).get("name") == test_name
            and msg.get("executable")
        ):
            return [msg["executable"]]
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Best-of-3 quiet-runner wall-time, peak RSS, and ASME V&V 20 |E| telemetry ledger."
    )
    parser.add_argument(
        "--trials",
        type=int,
        default=3,
        help="Number of repetitions per suite for best-of-N timing (default: 3).",
    )
    parser.add_argument(
        "--suite",
        action="append",
        default=[],
        help="Run only named suite(s); repeatable. Defaults to all validation/example/R22 suites.",
    )
    parser.add_argument(
        "--out",
        type=pathlib.Path,
        default=None,
        help="Write JSON telemetry ledger to this file path.",
    )
    parser.add_argument(
        "--append-jsonl",
        type=pathlib.Path,
        default=None,
        help="Append one compact JSON summary line to a historical JSONL ledger.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned suites and emit a dry-run telemetry ledger without compiling/running tests.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip initial `cargo test --release --no-run` pre-build step.",
    )
    parser.add_argument(
        "--cargo-target-dir",
        type=str,
        default=None,
        help="Optional CARGO_TARGET_DIR override for telemetry runs.",
    )
    args = parser.parse_args()

    if args.trials < 1:
        parser.error("--trials must be >= 1")

    selected = DEFAULT_SUITES
    if args.suite:
        wanted = set(args.suite)
        known = {s["name"] for s in DEFAULT_SUITES}
        unknown = sorted(wanted - known)
        if unknown:
            parser.error(f"unknown suite(s): {', '.join(unknown)} (valid: {', '.join(sorted(known))})")
        selected = [s for s in DEFAULT_SUITES if s["name"] in wanted]

    env = os.environ.copy()
    if args.cargo_target_dir:
        env["CARGO_TARGET_DIR"] = args.cargo_target_dir
    env.setdefault("NYX_QUEUE_LEASE_DEAD_SECS", "60")

    host_info = {
        "os": sys.platform,
        "system": platform.system(),
        "release": platform.release(),
        "arch": platform.machine(),
        "cpu_count": os.cpu_count(),
        "python": platform.python_version(),
        "nyx_threads": env.get("NYX_THREADS", "auto"),
        "nyx_factor_threads": env.get("NYX_FACTOR_THREADS", "auto"),
        "nyx_queue_lease_dead_secs": env.get("NYX_QUEUE_LEASE_DEAD_SECS"),
        "git_commit": _git_head(),
    }

    suite_entries: list[dict[str, Any]] = []
    all_targets: list[dict[str, Any]] = []
    any_failed = False

    if args.dry_run:
        for suite in selected:
            cmd = _suite_cargo_cmd(suite, no_run=False)
            print(f"[dry-run] {suite['name']}: {' '.join(cmd)}")
            suite_entries.append(
                {
                    "name": suite["name"],
                    "package": suite["package"],
                    "test": suite["test"],
                    "command": cmd,
                    "trials": args.trials,
                    "status": "dry_run",
                    "best_wall_s": None,
                    "median_wall_s": None,
                    "trial_walls_s": [],
                    "peak_rss_bytes": None,
                    "peak_rss_mb": None,
                    "vv20_targets": [],
                }
            )
    else:
        for suite in selected:
            exec_base = None if args.skip_build else _find_test_executable(suite["package"], suite["test"], env)
            if exec_base is not None:
                run_cmd = exec_base + suite["args"] + ["--nocapture"]
            else:
                run_cmd = _suite_cargo_cmd(suite, no_run=False)

            trial_walls: list[float] = []
            peak_rss: int | None = None
            peak_commit: int | None = None
            suite_targets: list[dict[str, Any]] = []
            rc_all = 0
            last_output = ""

            for t in range(args.trials):
                res = _run_once(run_cmd, env)
                trial_walls.append(round(res["wall_s"], 6))
                if res["peak_rss_bytes"] is not None:
                    peak_rss = max(peak_rss or 0, res["peak_rss_bytes"])
                if res["peak_commit_bytes"] is not None:
                    peak_commit = max(peak_commit or 0, res["peak_commit_bytes"])
                if res["returncode"] != 0 and rc_all == 0:
                    rc_all = res["returncode"]
                last_output = res["output"]
                if t == 0:
                    suite_targets = _parse_vv20_output(last_output)

            if rc_all != 0:
                any_failed = True
                print(last_output, file=sys.stderr)

            best_s = min(trial_walls)
            med_s = statistics.median(trial_walls)
            rss_mb = round(peak_rss / (1024.0 * 1024.0), 2) if peak_rss is not None else None
            status = "passed" if rc_all == 0 else "failed"
            print(
                f"[telemetry] {suite['name']:30s} status={status:6s} "
                f"best={best_s:.3f}s median={med_s:.3f}s peak_rss={rss_mb if rss_mb is not None else 'n/a'} MiB "
                f"vv20_targets={len(suite_targets)}"
            )

            for tgt in suite_targets:
                row = dict(tgt)
                row["suite"] = suite["name"]
                all_targets.append(row)

            suite_entries.append(
                {
                    "name": suite["name"],
                    "package": suite["package"],
                    "test": suite["test"],
                    "command": run_cmd,
                    "trials": args.trials,
                    "status": status,
                    "best_wall_s": best_s,
                    "median_wall_s": round(med_s, 6),
                    "trial_walls_s": trial_walls,
                    "peak_rss_bytes": peak_rss,
                    "peak_rss_mb": rss_mb,
                    "peak_commit_bytes": peak_commit,
                    "vv20_targets": suite_targets,
                }
            )

    max_rel = max((t["rel_e"] for t in all_targets), default=None)
    ledger = {
        "format": "vvuq_telemetry.v1",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "dry_run": args.dry_run,
        "trials": args.trials,
        "host": host_info,
        "suites": suite_entries,
        "accuracy_summary": {
            "target_count": len(all_targets),
            "all_passed": all(t["passed"] for t in all_targets) if all_targets else True,
            "max_rel_error": max_rel,
            "targets": all_targets,
        },
    }

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Wrote telemetry ledger: {args.out}")

    if args.append_jsonl:
        args.append_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.append_jsonl.open("a", encoding="utf-8") as fp:
            fp.write(json.dumps(ledger, sort_keys=True) + "\n")

    return 1 if any_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
