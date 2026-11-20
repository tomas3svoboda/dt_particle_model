"""Run one DTDC command with the owner's workstation reserve (2026-11-18).

Windows CPU affinity is inherited by descendants, including rg, Git and kernel
workers: https://learn.microsoft.com/en-us/windows/win32/procthread/inheritance
All guarded commands share one nonblocking lock and the same CPU subset.
The memory floor is a preflight check and a cooperative model stop request,
not an allocation limit or a promise to terminate Python during a solver step.

Usage: python scripts/run_with_resource_budget.py [--stop-file PATH] -- COMMAND ...
"""

from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / "scratchpad/resource_budget/heavy_command.lock"
FREE_MEMORY_PERCENT = 35.0
MINIMUM_FREE_BYTES = 25 * 1024**3  # Retain the owner's earlier 25 GiB reserve too.
BUSY_OR_RESOURCE_REFUSAL = 75
AFFINITY_ENV = "DTDC_RESOURCE_AFFINITY_MASK"


class MemoryStatus(ctypes.Structure):
    _fields_ = [("length", wintypes.DWORD), ("load", wintypes.DWORD)] + [
        (name, ctypes.c_ulonglong)
        for name in (
            "total_physical",
            "available_physical",
            "total_pagefile",
            "available_pagefile",
            "total_virtual",
            "available_virtual",
            "available_extended",
        )
    ]


def memory_snapshot():
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GlobalMemoryStatusEx.argtypes = [ctypes.POINTER(MemoryStatus)]
    kernel.GlobalMemoryStatusEx.restype = wintypes.BOOL
    state = MemoryStatus()
    state.length = ctypes.sizeof(state)
    if not kernel.GlobalMemoryStatusEx(ctypes.byref(state)):
        raise ctypes.WinError(ctypes.get_last_error())
    return {
        "total_bytes": state.total_physical,
        "available_bytes": state.available_physical,
        "free_percent": 100.0 * state.available_physical / state.total_physical,
    }


def apply_cpu_budget():
    if os.name != "nt":
        raise RuntimeError("this workstation resource guard requires Windows")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.GetProcessAffinityMask.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_size_t),
    ]
    kernel.GetProcessAffinityMask.restype = wintypes.BOOL
    kernel.SetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.c_size_t]
    kernel.SetProcessAffinityMask.restype = wintypes.BOOL
    kernel.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel.SetPriorityClass.restype = wintypes.BOOL
    kernel.GetActiveProcessorGroupCount.restype = wintypes.WORD
    if kernel.GetActiveProcessorGroupCount() != 1:
        raise RuntimeError("multi-group CPU budgeting is unsupported; no command launched")
    process = kernel.GetCurrentProcess()
    current, system = ctypes.c_size_t(), ctypes.c_size_t()
    if not kernel.GetProcessAffinityMask(process, ctypes.byref(current), ctypes.byref(system)):
        raise ctypes.WinError(ctypes.get_last_error())
    inherited = os.environ.get(AFFINITY_ENV)
    if (
        inherited is not None
        and int(inherited) == current.value
        and current.value.bit_count() <= system.value.bit_count() // 2
    ):
        mask = current.value
    else:
        bits = [1 << i for i in range(system.value.bit_length()) if current.value & (1 << i)]
        if system.value.bit_count() < 2:
            raise RuntimeError("cannot reserve CPU cores on a single-core system")
        mask = sum(bits[: max(1, len(bits) // 2)])
    if not kernel.SetProcessAffinityMask(process, mask):
        raise ctypes.WinError(ctypes.get_last_error())
    if not kernel.SetPriorityClass(process, subprocess.BELOW_NORMAL_PRIORITY_CLASS):
        raise ctypes.WinError(ctypes.get_last_error())
    os.environ[AFFINITY_ENV] = str(mask)
    return {
        "affinity_mask": mask,
        "allowed_logical_cpus": mask.bit_count(),
        "system_logical_cpus": system.value.bit_count(),
        "priority": "BelowNormal",
    }


def acquire_lock():
    import msvcrt

    LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    stream = LOCK_PATH.open("a+b")
    if stream.tell() == 0:
        stream.write(b"0")
        stream.flush()
    stream.seek(0)
    try:
        msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
    except OSError:
        stream.close()
        return None
    return stream  # Closing releases the OS lock, including after a station restart.


def memory_is_safe(memory):
    return (
        memory["free_percent"] >= FREE_MEMORY_PERCENT
        and memory["available_bytes"] >= MINIMUM_FREE_BYTES
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    budget, memory = apply_cpu_budget(), memory_snapshot()
    print(
        json.dumps(
            {
                "resource_budget": budget,
                "memory": memory,
                "minimum_free_memory_percent": FREE_MEMORY_PERCENT,
                "minimum_free_bytes": MINIMUM_FREE_BYTES,
            }
        ),
        flush=True,
    )
    if args.report_only:
        return 0
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command is required")
    if not memory_is_safe(memory) or (args.stop_file is not None and args.stop_file.exists()):
        print("Resource/stop-file preflight refused; no command launched", flush=True)
        return BUSY_OR_RESOURCE_REFUSAL
    if any(
        token == "--pool"
        and command[i + 1 : i + 2] != ["1"]
        or token.startswith("--pool=")
        and token != "--pool=1"
        for i, token in enumerate(command)
    ):
        parser.error("this workstation permits --pool 1 only")
    lock = acquire_lock()
    if lock is None:
        print("Another guarded command is active; no command launched", flush=True)
        return BUSY_OR_RESOURCE_REFUSAL
    try:
        env = dict(os.environ)
        env.update(
            PYTHONHASHSEED="1",
            PYTHONPATH=str(ROOT / "src"),
            OPENBLAS_NUM_THREADS="1",
            OMP_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            NUMEXPR_NUM_THREADS="1",
            RAYON_NUM_THREADS="1",
        )
        if Path(command[0]).name.lower() in ("rg", "rg.exe"):
            command[1:1] = ["--threads", "1"]
        child = subprocess.Popen(
            command, cwd=ROOT, env=env, creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS
        )
        print(json.dumps({"event": "launched", "child_pid": child.pid}), flush=True)
        stop_requested = False
        while child.poll() is None:
            if not stop_requested and not memory_is_safe(memory_snapshot()):
                if args.stop_file is not None:
                    args.stop_file.parent.mkdir(parents=True, exist_ok=True)
                    args.stop_file.touch()
                print(
                    json.dumps(
                        {
                            "event": "memory_reserve_low",
                            "stop_file_requested": args.stop_file is not None,
                            "action": "cooperative stop; Python is not killed",
                        }
                    ),
                    flush=True,
                )
                stop_requested = True
            time.sleep(0.5)
        return child.returncode
    finally:
        lock.close()


if __name__ == "__main__":
    raise SystemExit(main())
