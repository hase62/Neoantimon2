#!/usr/bin/env python3
"""Check the paper's frozen inputs and stage one replay under a runs directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

PAGE = Path(__file__).resolve().parent
PAYLOAD = PAGE / "payload"
SOURCE_NAME = "restored20api-502-files-v001.tar.gz"
SOURCE = PAYLOAD / "shared" / SOURCE_NAME
SOURCE_SHA256 = "af58b053bc57fcb1f46fb8378ee16173621579e4aec156cdb54767d26792486a"
SCOPES = ("mr11", "study6", "r09")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify():
    lock = json.loads((PAGE / "payload-lock.json").read_text())
    files = list(PAYLOAD.rglob("*"))
    check(not any(p.is_symlink() for p in files), "Symlink in frozen payload")
    actual = {p.relative_to(PAYLOAD).as_posix() for p in files if p.is_file()}
    check(actual == lock["files"].keys(), "Extra or missing payload files")
    for name, expected in lock["files"].items():
        path = PAYLOAD / name
        check(path.stat().st_size == expected["bytes"] and sha256(path) == expected["sha256"],
              f"Size or SHA-256 mismatch: {name}")

    check(sha256(SOURCE) == SOURCE_SHA256, "Fixed Core SHA-256 mismatch")
    expected_source = json.loads(
        (PAYLOAD / "mr11/source/RESTORED_SOURCE502_HASHES.json").read_text())
    observed_source = {}
    namespace, description = None, None
    with tarfile.open(SOURCE) as bundle:
        for item in bundle:
            components = Path(item.name).parts
            check(components and components[0] == "restored20api"
                  and ".." not in components and (item.isdir() or item.isfile()),
                  f"Unsafe Core member: {item.name}")
            if item.isdir():
                continue
            name = "/".join(components[1:])
            check(name and name not in observed_source, f"Duplicate Core member: {name}")
            with bundle.extractfile(item) as stream:
                data = stream.read()
            observed_source[name] = hashlib.sha256(data).hexdigest()
            if name == "NAMESPACE":
                namespace = data
            if name == "DESCRIPTION":
                description = data
    check(observed_source == expected_source and len(observed_source) == 502,
          "Fixed Core file inventory mismatch")
    check(namespace and sum(line.startswith(b"export(") for line in namespace.splitlines()) == 20,
          "Fixed Core must export exactly 20 APIs")
    check(description and b"Version: 0.3.0.9000" in description.splitlines(),
          "Fixed Core version mismatch")
    return len(actual)


def inventory(directory):
    return [{"path": path.relative_to(directory).as_posix(),
             "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in sorted(directory.rglob("*")) if path.is_file()]


def write_checksums(path, files):
    path.write_text("".join(f"{item['sha256']}  {item['path']}\n" for item in files),
                    encoding="utf-8", newline="\n")


def command(items, cwd=None, env=None):
    subprocess.run([str(item) for item in items], cwd=cwd, env=env, check=True)


def stage(scope, run_dir, pvac_python):
    # Existing drivers consume a fresh, exact-scope directory. Shared Core bytes
    # are materialized in that directory under the driver's original filename.
    candidate = run_dir / "candidate"
    shutil.copytree(PAYLOAD / scope, candidate)
    inner = {"mr11": candidate / "source",
             "study6": candidate / "study6/sources",
             "r09": candidate / "r09/R09_CANDIDATE/sources"}[scope]
    inner.mkdir(parents=True, exist_ok=True)
    source_name = SOURCE_NAME if scope == "mr11" else "patched.tar.gz"
    shutil.copyfile(SOURCE, inner / source_name)
    check(sha256(inner / source_name) == SOURCE_SHA256, "Staged Core SHA-256 mismatch")

    manifest_name = {"mr11": "MANIFEST.json", "study6": "NEW_CANDIDATE_MANIFEST.json",
                     "r09": "NEW_MANIFEST.json"}[scope]
    if scope != "mr11":
        subtree = inner.parent
        write_checksums(subtree / "SHA256SUMS.txt", inventory(subtree))
    members = inventory(candidate)
    (candidate / manifest_name).write_text(json.dumps({"files": members}, indent=2) + "\n",
                                           encoding="utf-8", newline="\n")
    write_checksums(candidate / "SHA256SUMS", inventory(candidate))
    if scope == "study6":
        command([sys.executable, "-I", candidate / "verify_restoration.py"], cwd=candidate)
    if scope == "r09":
        command([sys.executable, "-I", candidate / "verify_candidate.py"], cwd=candidate)
        # Keep a venv interpreter's symlink. Resolving it could select a different Python.
        command([sys.executable, "-B", candidate / "r09/prepare_existing_pvac.py",
                 "--python", pvac_python], cwd=candidate)
        env = dict(os.environ, PYTHONPATH=str(candidate / "r09/pvac-isolated"),
                   PYTHONDONTWRITEBYTECODE="1")
        command([pvac_python, "-B", "-c",
                 "import pvactools.lib.aggregate_all_epitopes"], cwd=candidate, env=env)
    return candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scope", choices=SCOPES)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--verify-only", action="store_true", help="Verify all payloads without writing")
    mode.add_argument("--runs-dir", type=Path, help="New scope-specific working directory under runs")
    parser.add_argument("--pvac-python", type=Path, help="Explicit pVACtools 7.1.5 interpreter for R09")
    args = parser.parse_args()

    count = verify()
    if args.verify_only:
        print(f"PASS: {count} payload files; fixed Core SHA-256, 502 files, 20 exports")
        return
    check(sys.version_info >= (3, 12), "Replay preparation requires Python 3.12+")
    run_dir = args.runs_dir.resolve()
    check(not run_dir.exists() and not run_dir.is_relative_to(PAGE),
          "Use a new runs directory outside the committed paper tree")
    check(run_dir.parent.is_dir(), "Create the parent runs directory first")
    check(shutil.disk_usage(run_dir.parent).free >= 25 * 1024**3, "25 GiB disk reserve required")
    interpreter = None
    if args.scope == "r09":
        check(args.pvac_python is not None, "R09 requires --pvac-python")
        interpreter = args.pvac_python.expanduser().absolute()
        check(interpreter.is_file() and os.access(interpreter, os.X_OK),
              "The specified pVACtools interpreter is not executable")
    else:
        check(args.pvac_python is None, "--pvac-python is only used for R09")
    check(shutil.which("R") and shutil.which("Rscript"), "R and Rscript are required")
    command(["Rscript", "--vanilla", "-e",
             "stopifnot(getRversion() >= '4.1'); "
             "for (p in c('jsonlite','digest','stringi','yaml','xml2')) "
             "stopifnot(requireNamespace(p, quietly=TRUE))"])
    run_dir.mkdir()
    candidate = stage(args.scope, run_dir, interpreter)
    print(f"Prepared {args.scope} at {candidate}; no study driver executed")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
