#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--enable-d5" ]]; then
  echo "D5 was not implemented in this D1-D4 delivery; no nonlinear solve started." >&2
  exit 2
fi
echo "D5 disabled. No nonlinear solve started."
