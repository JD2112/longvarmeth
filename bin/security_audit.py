#!/usr/bin/env python3
"""
LongVarMeth Pre-Commit & Pre-Push Security Audit Scanner
Author: Jyotirmoy Das / LongVarMeth Core Team
Version: 1.0.0

Comprehensive defensive security scanner for:
- API keys, credentials, tokens, and private cryptographic keys
- Personal Health Information (PHI) & Personally Identifiable Information (PII)
- Hardcoded local system paths and internal infrastructure IP addresses
- Sensitive / unignored binary artifacts and oversized files (>25MB)
- Shell and code injection risks (e.g. unsafe eval)

Usage:
    python3 bin/security_audit.py [--all | --staged-only] [--strict] [--json]
"""

import os
import sys
import re
import json
import subprocess
import argparse
from pathlib import Path

# ANSI Color codes for readable terminal feedback
RED = "\033[91m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
BLUE = "\033[94m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"

# Maximum allowable file size for code/text files in repo (25 MB)
MAX_FILE_SIZE_BYTES = 25 * 1024 * 1024

# Files & directories strictly excluded from the scan
EXCLUDE_DIRS = {
    ".git",
    "work",
    "results",
    "results_*",
    ".nextflow",
    ".pytest_cache",
    "__pycache__",
    "venv",
    ".venv",
    "review",  # Local reviewer notes
    "node_modules",
}

# File extensions to skip (binary formats)
SKIP_EXTENSIONS = {
    ".png", ".jpg", ".jpeg", ".gif", ".ico", ".svg",
    ".gz", ".tar", ".zip", ".bam", ".bai", ".cram", ".crai",
    ".tbi", ".csi", ".sif", ".img", ".pyc"
}

# ==============================================================================
# AUDIT PATTERNS (Update this dictionary whenever new security checks are added)
# ==============================================================================
AUDIT_PATTERNS = [
    # --- 1. Cloud & Platform Secrets ---
    {
        "id": "SEC-001",
        "category": "Secret",
        "name": "AWS Access Key ID",
        "regex": re.compile(r"\b(AKIA[0-9A-Z]{16})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-002",
        "category": "Secret",
        "name": "Google Cloud API Key",
        "regex": re.compile(r"\b(AIza[0-9A-Za-z\-_]{35})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-003",
        "category": "Secret",
        "name": "GitHub Personal Access Token",
        "regex": re.compile(r"\b(gh[pousr]_[A-Za-z0-9_]{36,255}|github_pat_[A-Za-z0-9_]{82})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-004",
        "category": "Secret",
        "name": "GitLab Personal Access Token",
        "regex": re.compile(r"\b(glpat-[0-9a-zA-Z\-_]{20,})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-005",
        "category": "Secret",
        "name": "Slack Token",
        "regex": re.compile(r"\b(xox[baprs]-[0-9a-zA-Z]{10,48})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-006",
        "category": "Secret",
        "name": "Private Cryptographic Key",
        "regex": re.compile(r"-----BEGIN (?:RSA |EC |DSA |OPENSSH |PGP )?PRIVATE KEY-----"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-007",
        "category": "Secret",
        "name": "OpenAI / Anthropic API Key",
        "regex": re.compile(r"\b(sk-[a-zA-Z0-9]{20,})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-008",
        "category": "Secret",
        "name": "HuggingFace Access Token",
        "regex": re.compile(r"\b(hf_[a-zA-Z0-9]{34,})\b"),
        "severity": "CRITICAL"
    },
    {
        "id": "SEC-009",
        "category": "Secret",
        "name": "Hardcoded Password / Auth Secret Assignment",
        "regex": re.compile(r'''(?i)(?:password|passwd|api_key|secret_key|auth_token)\s*[:=]\s*["']([^"'\s]{8,})["']'''),
        "severity": "HIGH",
        "allowlist_values": {"password", "secret", "token", "default", "none", "password123"}
    },

    # --- 2. Environment & Infrastructure Leakage ---
    {
        "id": "INF-001",
        "category": "Environment",
        "name": "Hardcoded User Filesystem Path",
        "regex": re.compile(r'''(?:['"]|^|\s)(/Users/[a-zA-Z0-9_\-]+|/home/[a-zA-Z0-9_\-]+|/data/[A-Z][a-zA-Z0-9_\-]+)(?:/|\b)'''),
        "severity": "HIGH",
        # Allow paths in documentation examples or audit logs if specified
        "file_allowlist": {"ROADMAP_RESOLVE.md", "ROADMAP.md", "bin/security_audit.py"}
    },
    {
        "id": "INF-002",
        "category": "Environment",
        "name": "Internal / RFC1918 Private IP Address",
        "regex": re.compile(r"\b(?:10\.\d{1,3}\.\d{1,3}\.\d{1,3}|172\.(?:1[6-9]|2\d|3[01])\.\d{1,3}\.\d{1,3}|192\.168\.\d{1,3}\.\d{1,3})\b"),
        "severity": "MEDIUM",
        "file_allowlist": {"bin/security_audit.py"}
    },

    # --- 3. PHI / PII & Clinical Identifiers ---
    {
        "id": "PHI-001",
        "category": "Clinical PHI",
        "name": "US Social Security Number Pattern",
        "regex": re.compile(r"\b\d{3}-\d{2}-\d{4}\b"),
        "severity": "HIGH",
        "file_allowlist": {"bin/security_audit.py"}
    },
    {
        "id": "PHI-002",
        "category": "Clinical PHI",
        "name": "UK NHS Number Pattern",
        "regex": re.compile(r"\b\d{3}\s?\d{3}\s?\d{4}\b"),
        "severity": "HIGH",
        "file_allowlist": {"bin/security_audit.py"}
    },

    # --- 4. Unsafe Code Constructs ---
    {
        "id": "COD-001",
        "category": "Code Safety",
        "name": "Unsafe Python eval() on Dynamic Input",
        "regex": re.compile(r"\beval\s*\([^\)]+\)"),
        "severity": "MEDIUM",
        "file_allowlist": {"bin/security_audit.py"}
    }
]

# Sensitive file names that should never be committed
PROHIBITED_FILE_PATTERNS = [
    re.compile(r"^\.env.*$"),
    re.compile(r"^.*\.pem$"),
    re.compile(r"^.*\.key$"),
    re.compile(r"^.*\.p12$"),
    re.compile(r"^.*\.pfx$"),
    re.compile(r"^id_rsa.*$"),
    re.compile(r"^id_ed25519.*$"),
    re.compile(r"^credentials$"),
    re.compile(r"^.*\.secret$")
]


def get_git_root() -> Path:
    """Resolve repository root using git."""
    try:
        out = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True)
        return Path(out.strip())
    except Exception:
        return Path(__file__).resolve().parent.parent


def get_staged_files(repo_root: Path) -> list:
    """Retrieve list of files staged in git."""
    try:
        out = subprocess.check_output(["git", "diff", "--cached", "--name-only"], text=True)
        lines = [repo_root / line.strip() for line in out.splitlines() if line.strip()]
        return [p for p in lines if p.exists() and p.is_file()]
    except Exception:
        return []


def get_all_scannable_files(repo_root: Path) -> list:
    """Discover all non-ignored, non-binary repository files using git ls-files."""
    try:
        out = subprocess.check_output(
            ["git", "ls-files", "-co", "--exclude-standard"],
            cwd=str(repo_root),
            text=True
        )
        files = []
        for line in out.splitlines():
            line = line.strip()
            if not line:
                continue
            path = repo_root / line
            if path.exists() and path.is_file():
                files.append(path)
        return files
    except Exception:
        # Fallback to filesystem walk if git is unavailable
        scannable = []
        for root, dirs, files in os.walk(repo_root):
            dirs[:] = [
                d for d in dirs
                if d not in EXCLUDE_DIRS and not any(d.startswith(prefix) for prefix in ["results_", "work_"])
            ]
            for file in files:
                path = Path(root) / file
                scannable.append(path)
        return scannable


def audit_repository(repo_root: Path, staged_only: bool = False) -> dict:
    """Executes full security audit on repository files."""
    files_to_scan = get_staged_files(repo_root) if staged_only else get_all_scannable_files(repo_root)
    findings = []
    prohibited_files_found = []
    oversized_files_found = []

    for file_path in files_to_scan:
        rel_path = file_path.relative_to(repo_root).as_posix()

        # 1. Check Prohibited Filenames
        for p_pat in PROHIBITED_FILE_PATTERNS:
            if p_pat.match(file_path.name):
                # Check if ignored by git
                try:
                    subprocess.check_call(["git", "check-ignore", "-q", str(file_path)])
                    is_ignored = True
                except subprocess.CalledProcessError:
                    is_ignored = False

                if not is_ignored:
                    prohibited_files_found.append({
                        "file": rel_path,
                        "rule": "Sensitive / Credential Filename",
                        "severity": "CRITICAL"
                    })

        # 2. Check File Size
        if file_path.exists():
            try:
                size = file_path.stat().st_size
                if size > MAX_FILE_SIZE_BYTES and not rel_path.startswith(".git"):
                    oversized_files_found.append({
                        "file": rel_path,
                        "size_mb": round(size / (1024 * 1024), 2),
                        "severity": "HIGH"
                    })
            except Exception:
                pass

        # 3. Skip binaries from content regex scanning
        if file_path.suffix.lower() in SKIP_EXTENSIONS or file_path.name.endswith(".tbi"):
            continue

        if not file_path.exists() or not file_path.is_file():
            continue

        # 4. Scan Text Content Line by Line
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                for line_idx, line in enumerate(f, 1):
                    for rule in AUDIT_PATTERNS:
                        # Check file allowlist
                        allowlisted = rule.get("file_allowlist", set())
                        if any(rel_path.endswith(allowed) or rel_path == allowed for allowed in allowlisted):
                            continue

                        match = rule["regex"].search(line)
                        if match:
                            matched_str = match.group(0).strip()
                            
                            # Check rule-specific value allowlist
                            val_allowlist = rule.get("allowlist_values", set())
                            if matched_str.lower() in val_allowlist:
                                continue

                            # Mask sensitive match for display
                            if rule["category"] == "Secret" and len(matched_str) > 8:
                                masked = matched_str[:4] + "*" * (len(matched_str) - 8) + matched_str[-4:]
                            else:
                                masked = matched_str

                            findings.append({
                                "rule_id": rule["id"],
                                "category": rule["category"],
                                "name": rule["name"],
                                "severity": rule["severity"],
                                "file": rel_path,
                                "line": line_idx,
                                "match": masked,
                                "context": line.strip()[:120]
                            })
        except Exception:
            pass

    return {
        "scanned_files_count": len(files_to_scan),
        "findings": findings,
        "prohibited_files": prohibited_files_found,
        "oversized_files": oversized_files_found
    }


def main():
    parser = argparse.ArgumentParser(description="LongVarMeth Security & Pre-Push Audit Scanner")
    parser.add_argument("--staged-only", action="store_true", help="Scan only git staged files")
    parser.add_argument("--strict", action="store_true", help="Fail with exit code 1 on WARNINGS or HIGH severity")
    parser.add_argument("--json", action="store_true", help="Output results as JSON")
    args = parser.parse_args()

    repo_root = get_git_root()

    if not args.json:
        print(f"\n{BOLD}{CYAN}{'='*70}{RESET}")
        print(f"{BOLD}{CYAN} LongVarMeth Pre-Push Security & Compliance Audit Scanner {RESET}")
        print(f"{BOLD}{CYAN}{'='*70}{RESET}")
        print(f"Repository Root : {repo_root}")
        print(f"Scan Mode       : {'Staged files only' if args.staged_only else 'Full workspace audit'}\n")

    results = audit_repository(repo_root, staged_only=args.staged_only)

    if args.json:
        print(json.dumps(results, indent=2))
        sys.exit(1 if results["findings"] or results["prohibited_files"] else 0)

    # Print Findings
    findings = results["findings"]
    prohibited = results["prohibited_files"]
    oversized = results["oversized_files"]

    critical_count = sum(1 for f in findings if f["severity"] == "CRITICAL") + len(prohibited)
    high_count = sum(1 for f in findings if f["severity"] == "HIGH") + len(oversized)
    medium_count = sum(1 for f in findings if f["severity"] == "MEDIUM")

    # 1. Prohibited File Violations
    if prohibited:
        print(f"{BOLD}{RED}[!] PROHIBITED / SENSITIVE FILES DETECTED:{RESET}")
        for p in prohibited:
            print(f"  {RED}✖ [{p['severity']}] {p['file']} -> {p['rule']}{RESET}")
        print()

    # 2. Oversized Files
    if oversized:
        print(f"{BOLD}{YELLOW}[!] OVERSIZED FILES DETECTED (>25MB):{RESET}")
        for o in oversized:
            print(f"  {YELLOW}⚠ [{o['severity']}] {o['file']} ({o['size_mb']} MB){RESET}")
        print()

    # 3. Content Findings
    if findings:
        print(f"{BOLD}{RED}[!] SECURITY ISSUES IDENTIFIED IN SOURCE CODE:{RESET}")
        for f in findings:
            color = RED if f["severity"] in ["CRITICAL", "HIGH"] else YELLOW
            print(f"  {color}✖ [{f['severity']}] {f['rule_id']}: {f['name']}{RESET}")
            print(f"    File: {BOLD}{f['file']}:{f['line']}{RESET}")
            print(f"    Match: {f['match']}")
            print(f"    Code : {f['context']}\n")

    # Summary
    print(f"{BOLD}{CYAN}{'-'*70}{RESET}")
    print(f"Scanned Files : {results['scanned_files_count']}")
    print(f"Critical      : {RED if critical_count else GREEN}{critical_count}{RESET}")
    print(f"High          : {YELLOW if high_count else GREEN}{high_count}{RESET}")
    print(f"Medium        : {BLUE if medium_count else GREEN}{medium_count}{RESET}")
    print(f"{BOLD}{CYAN}{'-'*70}{RESET}")

    if critical_count > 0 or (args.strict and high_count > 0):
        print(f"\n{BOLD}{RED}❌ AUDIT FAILED: Security violations must be resolved before pushing to GitHub.{RESET}\n")
        sys.exit(1)
    else:
        print(f"\n{BOLD}{GREEN}✅ AUDIT PASSED: No credentials, secrets, or security violations found. Safe to push.{RESET}\n")
        sys.exit(0)


if __name__ == "__main__":
    main()
