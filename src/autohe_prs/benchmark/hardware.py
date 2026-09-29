"""Hardware/runtime metadata captured with every benchmark."""

from __future__ import annotations

import os
from pathlib import Path
import platform
import struct
import sys
from typing import Any


def hardware_metadata() -> dict[str, Any]:
    memory_bytes: int | None = None
    physical_cores: int | None = None
    try:
        import psutil  # type: ignore[import-not-found]
    except ImportError:
        pass
    else:
        memory_bytes = int(psutil.virtual_memory().total)
        physical_cores = psutil.cpu_count(logical=False)
    if memory_bytes is None:
        try:
            if sys.platform == "win32":
                import ctypes

                class MemoryStatus(ctypes.Structure):
                    _fields_ = [
                        ("length", ctypes.c_ulong),
                        ("memory_load", ctypes.c_ulong),
                        ("total_physical", ctypes.c_ulonglong),
                        ("available_physical", ctypes.c_ulonglong),
                        ("total_page_file", ctypes.c_ulonglong),
                        ("available_page_file", ctypes.c_ulonglong),
                        ("total_virtual", ctypes.c_ulonglong),
                        ("available_virtual", ctypes.c_ulonglong),
                        ("available_extended_virtual", ctypes.c_ulonglong),
                    ]

                status = MemoryStatus()
                status.length = ctypes.sizeof(MemoryStatus)
                if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                    memory_bytes = int(status.total_physical)
            else:
                memory_bytes = int(
                    os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                )
        except (AttributeError, OSError, ValueError):
            memory_bytes = None
    if physical_cores is None and sys.platform.startswith("linux"):
        try:
            blocks = Path("/proc/cpuinfo").read_text(encoding="utf-8").split("\n\n")
            identifiers = set()
            for block in blocks:
                values = dict(
                    line.split(":", 1)
                    for line in block.splitlines()
                    if ":" in line
                )
                identifiers.add(
                    (
                        values.get("physical id", "").strip(),
                        values.get("core id", "").strip(),
                    )
                )
            identifiers.discard(("", ""))
            physical_cores = len(identifiers) or None
        except OSError:
            physical_cores = None
    if physical_cores is None and sys.platform == "win32":
        try:
            import ctypes

            required = ctypes.c_ulong(0)
            kernel32 = ctypes.windll.kernel32
            kernel32.GetLogicalProcessorInformationEx(
                0, None, ctypes.byref(required)
            )
            buffer = ctypes.create_string_buffer(required.value)
            if kernel32.GetLogicalProcessorInformationEx(
                0, buffer, ctypes.byref(required)
            ):
                offset = 0
                count = 0
                while offset < required.value:
                    relationship, size = struct.unpack_from("II", buffer.raw, offset)
                    if size <= 0:
                        break
                    if relationship == 0:  # RelationProcessorCore
                        count += 1
                    offset += size
                physical_cores = count or None
        except (AttributeError, OSError, struct.error, ValueError):
            physical_cores = None
    try:
        import openfhe  # type: ignore[import-not-found]
    except ImportError:
        openfhe_version = None
    else:
        openfhe_version = getattr(openfhe, "__version__", "unknown")
    docker = Path("/.dockerenv").exists()
    if not docker and sys.platform.startswith("linux"):
        try:
            docker = any(
                token in Path("/proc/1/cgroup").read_text(encoding="utf-8").lower()
                for token in ("docker", "containerd", "kubepods")
            )
        except OSError:
            pass
    wsl = "microsoft" in platform.release().lower()
    return {
        "cpu_model": platform.processor() or platform.machine(),
        "logical_core_count": os.cpu_count(),
        "physical_core_count": physical_cores,
        "ram_bytes": memory_bytes,
        "operating_system": platform.platform(),
        "docker_or_wsl": "docker" if docker else ("wsl" if wsl else "native"),
        "openfhe_version": openfhe_version,
        "python_version": sys.version.split()[0],
    }
