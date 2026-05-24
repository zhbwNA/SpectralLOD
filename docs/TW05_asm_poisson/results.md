Reproduction target: Table 3.1 in Toselli-Widlund, *Domain Decomposition Methods--Algorithms and Theory*, for Poisson additive Schwarz behavior.
Created: 2026-05-24
Updated: 2026-05-25
Verification entry point: `verify/verify_ddm_study.m`
Main utilities: `assembleStiffness2D`, `partitionMesh2D`, `additiveSchwarz`

# TW05 ASM Poisson Reproduction Summary

## Verification Entry Point

Run the reproduction study with:

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_study.m');"
```

The entry-point script is `verify/verify_ddm_study.m`. It prints both ASM and OSM tables; this report records the ASM/Toselli-Widlund target rows.

## Results (h=1/32..1/256)

| Dimension | Method | nSub | δ>0? | κ | ρ | Iters |
|---|---|---:|---|---:|---:|---:|
| 1D | ASM | 2 | yes | 1.0 | ~0 | 3 |
| 1D | ASM | 4 | yes | 1.2-1.3 | 0.002-0.004 | 5 |
| 1D | ASM | 8 | yes | 2.0-2.5 | 0.030-0.051 | 9 |
| 1D | ASM | any | **no** | 80-550 | 1.2-1.5 | **fails** |
| 2D | ASM | 2 | yes | 1.4-2.0 | 0.006-0.030 | 6-7 |
| 2D | ASM | 3 | yes | 2.6-5.2 | 0.056-0.151 | 9-13 |
| 2D | ASM | 4 | yes | 3.8-9.1 | 0.105-0.253 | 11-17 |
| 2D | ASM | any | **no** | 700-1600 | 1.1 | **fails** |
| 3D | ASM | 2 | yes | 2.0-2.5 | 0.029-0.050 | 7-8 |
| 3D | ASM | 3 | yes | 5.3 | 0.155 | 13 |
| 3D | ASM | any | **no** | 8000+ | 1.0+ | **fails** |

## Findings

1. ASM-PCG is effective across all tested dimensions when overlap is present.
2. Positive overlap is essential for the Dirichlet-inner-boundary ASM setup used here.
3. The measured condition-number behavior is consistent with dependence on `H/δ` rather than direct dependence on `h`.
