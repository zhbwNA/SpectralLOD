# MATLAB Function Usage Audit

Date: 2026-05-25

This is a conservative static scan of top-level MATLAB function files. A row means the function name was not found as a textual reference in another `.m` file in this repository. This does not prove the function is dead: scripts, external users, callbacks, generated commands, or manual MATLAB calls can still use it.

## No Static Callers

| Function | Path | Note |
|---|---|---|
| `orasSchwarz2D` | `src/DDM/Solvers/orasSchwarz2D.m` | Research solver/API candidate; no in-repo caller found. |
| `nedelec1_2D` | `src/FE/Nedelec/nedelec1_2D.m` | Basis evaluator/API candidate; may be used interactively or by future assemblers. |
| `nedelec1_3D` | `src/FE/Nedelec/nedelec1_3D.m` | Basis evaluator/API candidate; may be used interactively or by future assemblers. |
| `nedelec2_3D` | `src/FE/Nedelec/nedelec2_3D.m` | Basis evaluator/API candidate; may be used interactively or by future assemblers. |
| `rasPMLHelmholtz2D` | `src/Preconditioners/rasPMLHelmholtz2D.m` | Alias/API entry point; no in-repo caller found. |
| `shiftedLaplacianPreconditionerMaxwell2D` | `src/Preconditioners/shiftedLaplacianPreconditionerMaxwell2D.m` | Maxwell preconditioner wrapper; no in-repo caller found. |
| `prolongate_P2_P2` | `src/Utils/Transfer/prolongate_P2_P2.m` | Transfer utility; no in-repo caller found. |
| `plot_partition_asm2d` | `verify/plot_partition_asm2d.m` | Manual plotting entry point. |
| `plot_partition1d` | `verify/plot_partition1d.m` | Manual plotting entry point. |
| `plot_partition2d` | `verify/plot_partition2d.m` | Manual plotting entry point. |
| `plot_partition3d_slice` | `verify/plot_partition3d_slice.m` | Manual plotting entry point. |
| `verify_all` | `verify/verify_all.m` | Master verification entry point. |
| `verify_gggls_pml_decay_convergence` | `verify/verify_gggls_pml_decay_convergence.m` | Reproduction entry point. |
| `verify_gggls_ras_pml_reproduce` | `verify/verify_gggls_ras_pml_reproduce.m` | Reproduction entry point. |
| `verify_oras_largek_iterations` | `verify/verify_oras_largek_iterations.m` | Large-run verification entry point. |
| `verify_oras_reproduce_gggs` | `verify/verify_oras_reproduce_gggs.m` | Reproduction entry point. |
| `verify_theory_new_functions` | `verify/verify_theory_new_functions.m` | Interest-driven theory verification entry point. |
