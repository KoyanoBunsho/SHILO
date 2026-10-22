"""Audit every locked package, including platform and Python-specific variants.

Run with: uv run --no-project --with pip-audit python scripts/audit_dependencies.py
"""

from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib


def main():
    lock = Path(__file__).resolve().parents[1] / "poetry.lock"
    packages = tomllib.loads(lock.read_text())["package"]
    # pip-audit rejects two versions of a name in one requirements file.
    # Audit alternate versions separately instead of dropping their markers.
    batches = []
    for package in packages:
        name = package["name"]
        line = f"{name}=={package['version']}"
        line += "".join(f" --hash={file['hash']}" for file in package["files"])
        batch = next((batch for batch in batches if name not in batch), None)
        if batch is None:
            batch = {}
            batches.append(batch)
        batch[name] = line

    status = 0
    with tempfile.TemporaryDirectory(prefix="shilo-audit-") as directory:
        for index, batch in enumerate(batches, 1):
            requirements = Path(directory) / f"requirements-{index}.txt"
            requirements.write_text("\n".join(batch.values()) + "\n")
            print(f"Auditing batch {index}/{len(batches)}: {len(batch)} locked packages", flush=True)
            result = subprocess.run(
                [sys.executable, "-m", "pip_audit", "--disable-pip",
                 "--require-hashes", "-r", str(requirements)],
                check=False,
            )
            if result.returncode:
                status = result.returncode
    return status


if __name__ == "__main__":
    raise SystemExit(main())
