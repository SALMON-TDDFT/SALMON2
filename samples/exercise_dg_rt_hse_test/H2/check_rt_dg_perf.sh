#!/bin/bash
set -euo pipefail
LOG="${1:?run.log path required}"
rg -n "total calculation time," "$LOG"
rg -n "hamiltonian module -> stencil|stencil" "$LOG" | tail -n 1
if rg -n "SIGSEGV|Program received signal" "$LOG"; then
  echo "fatal signal detected" >&2
  exit 1
fi
