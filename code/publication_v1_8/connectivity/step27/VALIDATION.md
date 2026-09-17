# Validation notes

The packaged code was checked with:

```text
python -m py_compile step27_curved_bridge_optimization.py step25_engine_snapshot.py
python step27_curved_bridge_optimization.py --self-test
```

The self-test verifies exact preservation of kinetic and auxiliary endpoint values under the curved sine-basis path parameterization and basic graph-union behavior.

A reduced end-to-end smoke run was also executed against the uploaded Step-26 result using a synthetic Step-25 accepted-broad geometry. That test exercised:

```text
input loading
Step-25 engine replay
native coarse-component reconstruction
reuse of Step-26 refined edges
new direct-path verification
curved native bridge optimization
final graph assembly
figure generation
scientific-decision output
```

The synthetic smoke result is not a scientific result and is not included in the release archive.
