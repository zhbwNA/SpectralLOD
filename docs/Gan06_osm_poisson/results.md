Reproduction target: Gander, "Optimized Schwarz Methods", for Poisson OSM convergence behavior and optimized Robin transmission.
Created: 2026-05-24
Updated: 2026-05-25
Verification entry point: `verify/verify_ddm_study.m`; extension: `verify/verify_ddm_extension.m`
Main utilities: `assembleStiffness2D`, `partitionMesh2D`, `optimizedSchwarzPoisson2D`, `optimizedSchwarzPoisson2D_overlap`, `assembleBoundaryMass2D`

# Gan06 OSM Poisson Reproduction Summary

## Verification Entry Points

Run the main OSM reproduction study with:

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_study.m');"
```

The main entry-point script is `verify/verify_ddm_study.m`. It prints both ASM and OSM tables; this report records the OSM/Gander target rows.

Run the 2D extension study with:

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_extension.m');"
```

The extension entry-point script is `verify/verify_ddm_extension.m`.

## Results (h=1/32..1/256 in 1D, h=1/12..1/24 in 2D, h=1/6..1/8 in 3D)

| Dimension | Method | nSub | ρ | Iters |
|---|---|---:|---:|---:|
| 1D | OSM | 2 | 0.15 | 9 |
| 1D | OSM | 4 | 0.46 | 18 |
| 1D | OSM | 8 | 0.91 | 100 |
| 2D | OSM | 2 | 0.94 | >100 |
| 2D | OSM | 3 | 0.99 | >100 |
| 3D | OSM | 2 | 0.96 | >100 |

## Findings

1. The tested one-level OSM converges in 1D with optimized Robin parameter `α`.
2. The same one-level method degrades in 2D/3D without a coarse space.
3. This file records the OSM/Gander reproduction target separately from the ASM/Toselli-Widlund checks.

## Extension Study

The extension study on 2D Poisson includes:
- OSM strip vs checkerboard (`ρ`, iterations)
- Overlapping OSM `δ` effect (`ρ` vs `δ`, `α`)
- Two-level OSM with coarse space

Extension results:
- Overlapping OSM: `δ=0.16` gives `ρ=0.04` (6 iterations), compared with `ρ=0.94` for non-overlapping.
- Two-level OSM: P1 coarse `H=1/6` gives `ρ=0.17` (9 iterations), compared with `ρ=0.996` for one-level.
