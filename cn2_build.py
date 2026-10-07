"""Portable build and launch helpers for the CN2 engines.

compile_cmd: $CC if set, else clang if present, else gcc; tune for the host CPU with -mcpu=native on
             ARM (Apple Silicon, aarch64) and -march=native elsewhere.  On an Apple Silicon Mac this gives
             exactly `clang -O3 -mcpu=native`, the command the census engines were built with.
keep_awake:  prefix `caffeinate -i` on macOS so a long run is not interrupted by sleep; elsewhere, nothing.
"""
import os
import platform
import shutil


def compile_cmd(src, out, extra=()):
    cc = os.environ.get('CC') or ('clang' if shutil.which('clang') else 'gcc')
    arm = platform.machine().lower() in ('arm64', 'aarch64')
    return [cc, '-O3', '-mcpu=native' if arm else '-march=native', *extra, '-o', out, src, '-lpthread', '-lm']


def keep_awake(cmd):
    return ['caffeinate', '-i', *cmd] if shutil.which('caffeinate') else list(cmd)
