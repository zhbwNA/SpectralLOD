Reproduction target: Gander-Graham-Spence shifted-Laplacian GMRES paper, Tables 4-5.
Created: 2026-05-25
Updated: 2026-05-25
Verification entry point: `verify/verify_ggs15_shifted_laplacian_tables45.m`
Main utilities: `squaremesh`, `assembleHelmholtz2D`, `shiftedLaplacianPDE`, `shiftedLaplacianPreconditioner2D`, MATLAB `gmres`

# GGS15 Shifted-Laplacian GMRES Tables 4-5

Source: Gander, Graham, and Spence, Numerische Mathematik 131 (2015), Tables 4-5. Author manuscript: <https://purehost.bath.ac.uk/ws/portalfiles/portal/109280714/GaGrSp14_revised.pdf>.

This report records a repo-local reproduction path for Tables 4 and 5 of the shifted-Laplacian GMRES paper. The paper reports the distance `d` from the origin to the numerical range of the left-preconditioned matrix and unrestarted GMRES iterations to relative residual `1e-6`.

## Mathematical Form

The model problem is the 2D interior impedance Helmholtz problem on `Omega=(0,1)^2`. The repo discretization uses continuous P1 finite elements on the uniform triangulation returned by `squaremesh([0,1,0,1],1/n)`. The unshifted matrix is

```text
A = K - k^2 M - i k M_boundary.
```

For a shift `epsilon > 0`, the shifted-Laplacian preconditioner matrix is

```text
A_epsilon = K - (k^2 + i epsilon) M - i sqrt(k^2 + i epsilon) M_boundary.
```

The verified linear system is the left-preconditioned system

```text
A_epsilon^{-1} A u = A_epsilon^{-1} b,
```

with `b = ones(N,1)`, zero initial guess, and unrestarted MATLAB `gmres` tolerance `1.0e-06`. This follows the paper table setup except for the numerical-range distance, which is computed here only as a dense angular approximation when `N <= GGS15_D_MAXDOF`.

## Paper Parameters

- Table 4 mesh rule: `n = 2*k`.
- Table 5 mesh rule: `n = ceil(k^(3/2))`.
- Shift columns: `epsilon = k/4, k/2, k, 2k, 4k, k^(3/2), k^2`.
- Default configured Table 4 rows: `[10 20]`.
- Default configured Table 5 rows: `10`.
- Default maximum DOFs: `20000`; dense `d` approximation maximum DOFs: `700`.

### Table 4 Paper Targets: fixed points per wavelength

| k | metric | k/4 | k/2 | k | 2k | 4k | k^(3/2) | k^2 |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 10 | d | 0.933 | 0.871 | 0.764 | 0.597 | 0.386 | 0.459 | 0.147 |
| 10 | it | 4 | 5 | 6 | 7 | 9 | 8 | 13 |
| 20 | d | 0.927 | 0.862 | 0.749 | 0.580 | 0.373 | 0.341 | 0.054 |
| 20 | it | 4 | 5 | 6 | 8 | 10 | 11 | 25 |
| 40 | d | 0.925 | 0.857 | 0.741 | 0.568 | 0.359 | 0.231 | 0.016 |
| 40 | it | 4 | 5 | 6 | 8 | 11 | 13 | 47 |
| 80 | d | 0.923 | 0.854 | 0.736 | 0.561 | 0.352 | 0.148 | 0.004 |
| 80 | it | 4 | 5 | 6 | 7 | 10 | 16 | 84 |
| 160 | d | 0.922 | 0.853 | 0.734 | 0.555 | 0.344 | 0.087 | 0.003 |
| 160 | it | 4 | 5 | 6 | 7 | 10 | 19 | 148 |

### Table 5 Paper Targets: fixed scaled points per wavelength

| k | metric | k/4 | k/2 | k | 2k | 4k | k^(3/2) | k^2 |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 10 | d | 0.932 | 0.871 | 0.763 | 0.594 | 0.381 | 0.455 | 0.143 |
| 10 | it | 4 | 5 | 6 | 7 | 9 | 8 | 13 |
| 20 | d | 0.926 | 0.860 | 0.746 | 0.575 | 0.370 | 0.337 | 0.052 |
| 20 | it | 4 | 5 | 6 | 8 | 11 | 11 | 24 |
| 40 | d | 0.923 | 0.854 | 0.736 | 0.561 | 0.353 | 0.228 | 0.015 |
| 40 | it | 4 | 5 | 6 | 8 | 11 | 14 | 48 |
| 80 | d | 0.920 | 0.849 | 0.728 | 0.550 | 0.342 | 0.144 | 0.006 |
| 80 | it | 4 | 5 | 6 | 8 | 10 | 16 | 86 |

## Repo Run Results

### Table 4 Configured Run

| k | n | DOFs | epsilon | paper d | repo d | paper it | repo it | flag | relres | note |
|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|:---|
| 10 | 20 | 441 | k/4 | 0.933 | 0.930 | 4 | 4 | 0 | 7.42e-08 | ran |
| 10 | 20 | 441 | k/2 | 0.871 | 0.866 | 5 | 5 | 0 | 4.99e-08 | ran |
| 10 | 20 | 441 | k | 0.764 | 0.755 | 6 | 6 | 0 | 8.66e-08 | ran |
| 10 | 20 | 441 | 2k | 0.597 | 0.585 | 7 | 7 | 0 | 4.00e-07 | ran |
| 10 | 20 | 441 | 4k | 0.386 | 0.374 | 9 | 9 | 0 | 2.07e-07 | ran |
| 10 | 20 | 441 | k^(3/2) | 0.459 | 0.448 | 8 | 8 | 0 | 3.63e-07 | ran |
| 10 | 20 | 441 | k^2 | 0.147 | 0.139 | 13 | 13 | 0 | 3.40e-07 | ran |
| 20 | 40 | 1681 | k/4 | 0.927 | -- | 4 | 4 | 0 | 1.79e-07 | ran |
| 20 | 40 | 1681 | k/2 | 0.862 | -- | 5 | 5 | 0 | 9.93e-08 | ran |
| 20 | 40 | 1681 | k | 0.749 | -- | 6 | 6 | 0 | 1.61e-07 | ran |
| 20 | 40 | 1681 | 2k | 0.580 | -- | 8 | 7 | 0 | 7.78e-07 | ran |
| 20 | 40 | 1681 | 4k | 0.373 | -- | 10 | 10 | 0 | 4.12e-07 | ran |
| 20 | 40 | 1681 | k^(3/2) | 0.341 | -- | 11 | 10 | 0 | 9.70e-07 | ran |
| 20 | 40 | 1681 | k^2 | 0.054 | -- | 25 | 24 | 0 | 5.98e-07 | ran |

### Table 5 Configured Run

| k | n | DOFs | epsilon | paper d | repo d | paper it | repo it | flag | relres | note |
|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|:---|
| 10 | 32 | 1089 | k/4 | 0.932 | -- | 4 | 4 | 0 | 6.38e-08 | ran |
| 10 | 32 | 1089 | k/2 | 0.871 | -- | 5 | 4 | 0 | 9.39e-07 | ran |
| 10 | 32 | 1089 | k | 0.763 | -- | 6 | 6 | 0 | 5.20e-08 | ran |
| 10 | 32 | 1089 | 2k | 0.594 | -- | 7 | 7 | 0 | 2.92e-07 | ran |
| 10 | 32 | 1089 | 4k | 0.381 | -- | 9 | 9 | 0 | 2.05e-07 | ran |
| 10 | 32 | 1089 | k^(3/2) | 0.455 | -- | 8 | 8 | 0 | 3.65e-07 | ran |
| 10 | 32 | 1089 | k^2 | 0.143 | -- | 13 | 13 | 0 | 2.21e-07 | ran |

## Re-run Controls

PowerShell examples:

```powershell
$env:GGS15_TABLE4_KVALS = '10 20 40'
$env:GGS15_TABLE5_KVALS = '10 20'
$env:GGS15_MAXDOF = '50000'
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ggs15_shifted_laplacian_tables45.m');"
```

Set `GGS15_COMPUTE_D=0` to skip dense numerical-range estimates. Set `GGS15_D_MAXDOF` and `GGS15_D_ANGLES` higher only for small matrices, because the approximation forms the dense preconditioned matrix and diagonalizes angular Hermitian parts.

## Status

The verifier now gives a self-contained reproduction harness for Tables 4-5 and stores the paper target values in code and in this report. The default run is a small executable subset of the paper tables; larger rows should be run explicitly after choosing the resource envelope.
