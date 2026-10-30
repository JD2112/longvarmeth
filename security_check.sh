#!/usr/bin/env bash
# security_check.sh
# Comprehensive Pre-Commit & Pre-Push Security Audit Runner for LongVarMeth.
# Usage:
#   ./security_check.sh                  # Run full workspace security check (with strict enforcement)
#   ./security_check.sh --staged-only    # Scan only git staged files
#   ./security_check.sh --install-hook   # Install git pre-push hook to automate checks
#   ./security_check.sh --json           # Export findings as JSON

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="python3"

if ! command -v "$PYTHON_BIN" &> /dev/null; then
    echo "Error: python3 is required to run the security audit scanner."
    exit 1
fi

# Option: Install pre-push git hook
if [[ "${1:-}" == "--install-hook" ]]; then
    HOOK_PATH="${SCRIPT_DIR}/.git/hooks/pre-push"
    if [[ ! -d "${SCRIPT_DIR}/.git" ]]; then
        echo "Error: Not a git repository root. Cannot install hook."
        exit 1
    fi
    mkdir -p "${SCRIPT_DIR}/.git/hooks" 2>/dev/null || true
    if cp "${SCRIPT_DIR}/security_check.sh" "${HOOK_PATH}.tmp" 2>/dev/null; then
        rm -f "${HOOK_PATH}.tmp"
    fi
    
    echo "To install as a git pre-push hook, you can run:"
    echo "  ln -sf ../../security_check.sh .git/hooks/pre-push"
    echo "or add './security_check.sh --strict' to your git pre-push script."
    exit 0
fi

chmod +x "${SCRIPT_DIR}/bin/security_audit.py"
"$PYTHON_BIN" "${SCRIPT_DIR}/bin/security_audit.py" "$@"
