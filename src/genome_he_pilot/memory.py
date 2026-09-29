"""Memory measurement helpers for benchmark scripts."""

from __future__ import annotations

from dataclasses import dataclass
import ctypes
from ctypes import wintypes
import os
import tracemalloc


@dataclass(frozen=True)
class ProcessMemory:
    rss_bytes: int
    peak_rss_bytes: int


@dataclass(frozen=True)
class MemoryMetrics:
    rss_before_bytes: int
    rss_after_bytes: int
    rss_delta_bytes: int
    peak_rss_bytes: int
    python_current_bytes: int
    python_peak_bytes: int


class _ProcessMemoryCounters(ctypes.Structure):
    _fields_ = [
        ("cb", ctypes.c_ulong),
        ("PageFaultCount", ctypes.c_ulong),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


def process_memory() -> ProcessMemory:
    """Return current and peak resident memory for the running process."""

    if os.name == "nt":
        counters = _ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)

        get_current_process = ctypes.windll.kernel32.GetCurrentProcess
        get_current_process.restype = wintypes.HANDLE
        get_process_memory_info = ctypes.windll.psapi.GetProcessMemoryInfo
        get_process_memory_info.argtypes = [
            wintypes.HANDLE,
            ctypes.POINTER(_ProcessMemoryCounters),
            wintypes.DWORD,
        ]
        get_process_memory_info.restype = wintypes.BOOL

        handle = get_current_process()
        ok = get_process_memory_info(
            handle,
            ctypes.byref(counters),
            counters.cb,
        )
        if not ok:
            return ProcessMemory(rss_bytes=0, peak_rss_bytes=0)
        return ProcessMemory(
            rss_bytes=int(counters.WorkingSetSize),
            peak_rss_bytes=int(counters.PeakWorkingSetSize),
        )

    return ProcessMemory(rss_bytes=0, peak_rss_bytes=0)


def start_python_memory_trace() -> ProcessMemory:
    """Start Python allocation tracing and return the initial process memory."""

    tracemalloc.start()
    return process_memory()


def stop_process_memory_measurement(before: ProcessMemory) -> MemoryMetrics:
    """Return process memory metrics without Python allocation tracing overhead."""

    after = process_memory()
    return MemoryMetrics(
        rss_before_bytes=before.rss_bytes,
        rss_after_bytes=after.rss_bytes,
        rss_delta_bytes=after.rss_bytes - before.rss_bytes,
        peak_rss_bytes=after.peak_rss_bytes,
        python_current_bytes=0,
        python_peak_bytes=0,
    )


def stop_python_memory_trace(before: ProcessMemory) -> MemoryMetrics:
    """Stop Python allocation tracing and return combined memory metrics."""

    python_current, python_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    after = process_memory()
    return MemoryMetrics(
        rss_before_bytes=before.rss_bytes,
        rss_after_bytes=after.rss_bytes,
        rss_delta_bytes=after.rss_bytes - before.rss_bytes,
        peak_rss_bytes=after.peak_rss_bytes,
        python_current_bytes=python_current,
        python_peak_bytes=python_peak,
    )
