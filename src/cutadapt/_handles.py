"""
Wrappers around multiprocessing.reduction.send_handle/recv_handle
that
- fix a bug in CPython’s standard library preventing send_handle from
working on Windows
- and translate between Windows handles and Unix file descriptors as necessary.
"""

import multiprocessing
import os

try:
    import msvcrt
except ImportError:
    msvcrt = None  # type: ignore
try:
    import _winapi
except ImportError:
    _winapi = None  # type: ignore


class _PatchedDupHandle:
    """Used by _patched_send_handle"""

    def __init__(self, handle, access, pid=None, options=0):
        if pid is None:
            pid = os.getpid()
        proc = _winapi.OpenProcess(_winapi.PROCESS_DUP_HANDLE, False, pid)
        try:
            self._handle = _winapi.DuplicateHandle(
                _winapi.GetCurrentProcess(), handle, proc, access, False, options
            )
        finally:
            _winapi.CloseHandle(proc)
        self._access = access
        self._pid = pid

    def detach(self):
        if self._pid == os.getpid():
            return self._handle
        proc = _winapi.OpenProcess(_winapi.PROCESS_DUP_HANDLE, False, self._pid)
        try:
            return _winapi.DuplicateHandle(
                proc,
                self._handle,
                _winapi.GetCurrentProcess(),
                self._access,
                False,
                _winapi.DUPLICATE_CLOSE_SOURCE,
            )
        finally:
            _winapi.CloseHandle(proc)


def _patched_send_handle(conn, handle, destination_pid):
    """
    A patched version of multiprocessing.reduction.send_handle that works around
    bug https://github.com/python/cpython/issues/82369
    Adapted from code posted by Cameron Kennedy (m3rc1fulcameron) in that issue.
    """
    dh = _PatchedDupHandle(handle, 0, destination_pid, _winapi.DUPLICATE_SAME_ACCESS)
    conn.send(dh)


def send_handle(conn, fd, destination_pid):
    if _winapi is None:
        return multiprocessing.reduction.send_handle(conn, fd, destination_pid)
    else:
        handle = msvcrt.get_osfhandle(fd)
        return _patched_send_handle(conn, handle, destination_pid)


def recv_handle(conn):
    handle = multiprocessing.reduction.recv_handle(conn)
    if msvcrt:
        return msvcrt.open_osfhandle(handle, os.O_RDONLY | os.O_BINARY)
    else:
        return handle
