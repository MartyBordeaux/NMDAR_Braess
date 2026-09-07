#!/usr/bin/env python3
"""Run the exact, hash-verified Step11 v1.2 extraction source.

The XZ file is standard lossless compression of the original Python source.
Raw ABFs are required for extraction; --self-test uses synthetic traces only.
"""
import hashlib
import lzma
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SHA = "b994c232f779d1a2a04940530c359c9be17a9895b35617b03ddd2cfc34476419"

def main():
    source = lzma.decompress((HERE / "step11_experimental_controls_v1_2.py.xz").read_bytes())
    if hashlib.sha256(source).hexdigest() != SHA:
        raise RuntimeError("Frozen Step11 source hash mismatch")
    args = sys.argv[1:]
    if not any(x == "--output" or x.startswith("--output=") for x in args):
        args += ["--output", str(HERE.parents[1] / "reproduced/step11_raw_v1_2")]
    with tempfile.TemporaryDirectory(prefix="step11_frozen_") as tmp:
        script = Path(tmp) / "step11_experimental_controls_v1_2.py"
        script.write_bytes(source)
        subprocess.run([sys.executable, str(script), *args], check=True)

if __name__ == "__main__":
    main()
