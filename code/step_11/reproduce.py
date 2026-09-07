#!/usr/bin/env python3
"""Verify and rebuild Step11 v1.2 publication tables/figures from saved data.

No raw extraction, new thresholds, model runs, or inferential tests are used.
The two archive parts are combined and checked before safe member extraction.
"""
import argparse
import gzip
import hashlib
import io
import json
import lzma
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "data/derived/step11_v1_2"
ARCHIVE_SHA = "7ed39e3e12dfd9c4c0455cef8a705956ba1a11f9f92fa7ed44525cfa9d846303"
AUDIT_SHA = "9d85579f62eacc48947bd6c7323ab1e358a587bbf46ce797fd77ed0cb6635fab"


def verify(payload, expected, label):
    if hashlib.sha256(payload).hexdigest() != expected:
        raise RuntimeError("SHA-256 mismatch: " + label)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=ROOT / "reproduced/step11_publication")
    args = ap.parse_args()
    out = args.output.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    archive = b"".join((DATA / ("tables.tar.xz." + n)).read_bytes() for n in ("001", "002"))
    verify(archive, ARCHIVE_SHA, "saved measurement archive")
    manifest = json.loads((DATA / "PUBLIC_DATA_MANIFEST.json").read_text())
    source = lzma.decompress((HERE / "build_publication.py.xz").read_bytes())
    verify(source, AUDIT_SHA, "publication audit source")
    with tempfile.TemporaryDirectory(prefix="step11_publication_") as tmp:
        tmp = Path(tmp)
        inputs = tmp / "inputs"
        inputs.mkdir()
        found = set()
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:xz") as tf:
            for member in tf.getmembers():
                name = member.name
                if not member.isfile() or Path(name).name != name or name not in manifest or name in found:
                    raise RuntimeError("Unexpected archive member: " + name)
                payload = tf.extractfile(member).read()
                verify(payload, manifest[name]["sha256"], name)
                if len(payload) != manifest[name]["size_bytes"]:
                    raise RuntimeError("Size mismatch: " + name)
                target = manifest[name]["original_file"]
                if Path(target).name != target:
                    raise RuntimeError("Unsafe target name")
                if target.endswith(".gz"):
                    payload = gzip.compress(payload, mtime=0)
                (inputs / target).write_bytes(payload)
                found.add(name)
        if found != set(manifest):
            raise RuntimeError("Incomplete archive")
        script = tmp / "build_publication.py"
        script.write_bytes(source)
        subprocess.run([sys.executable, str(script), "--input", str(inputs), "--output", str(out)], check=True)
    (out / "public_archive_verification.json").write_text(json.dumps({
        "status": "PASS", "archive_sha256": ARCHIVE_SHA, "audit_source_sha256": AUDIT_SHA,
        "measurement_files_verified": len(found),
        "path_prefix_change": "/root/nmda/IV_NMDA/ -> IV_NMDA/",
        "note": "Temporary gzip encodings may differ; decompressed measurement bytes are verified."
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
