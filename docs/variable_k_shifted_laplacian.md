# Variable Wavenumber And Shifted-Laplacian Wrappers

The scalar Helmholtz assemblers now accept either the old constant `k` input or a PDE structure:

```matlab
pde = helmholtzPDE(@(x,y) 20 + 2*x, ...
    'source', @(x,y) exp(-40*((x-0.5).^2 + (y-0.5).^2)), ...
    'eta', 'k');
[A,b] = assembleHelmholtz2D(node, elem, bdFlag, pde, [], [], 1);
```

For shifted-Laplacian preconditioning, use

```matlab
pdeShift = shiftedLaplacianPDE(k, 'epsilon', 'quadratic', 'eta', 'sqrt');
Aeps = assembleHelmholtz2D(node, elem, bdFlag, pdeShift, [], [], degree);
```

This implements

```text
Aeps = K - int (k(x)^2 + i*epsilon(x)) phi_i phi_j dx
       - i int_Gamma eta(x) phi_i phi_j ds.
```

The supported named rules are:

| field | rule | meaning |
|---|---|---|
| `epsilon` | `'zero'` | no volume absorption |
| `epsilon` | `'linear'` | `epsilon = abs(k)` |
| `epsilon` | `'quadratic'` | `epsilon = abs(k)^2` |
| `eta` | `'k'` | `eta = k` |
| `eta` | `'sqrt'` | `eta = sqrt(k^2+i*epsilon)` |
| `eta` | `'zero'` | no boundary absorption |

Both `epsilon` and `eta` can also be numbers or function handles. Function handles may depend on `(x,y)`, `(x,y,k)`, or just `(k)`.

The shifted-Laplacian helper

```matlab
[applyPrecon,Aeps,pdeShift] = shiftedLaplacianPreconditioner2D(node, elem, bdFlag, k, degree, opts);
```

returns a left-preconditioner applying `Aeps\r` via sparse LU by default. The same interface is available for 3D Helmholtz as `shiftedLaplacianPreconditioner3D`. This follows the shifted form described by Gander, Graham, and Spence: `Aeps = S-(k^2+i*epsilon)M-i*eta*N`, with `eta` usually chosen as either `k` or `sqrt(k^2+i*epsilon)`.

For 2D Maxwell NE_1 experiments, use

```matlab
Amax = assembleMaxwell2D(node, elem, bdFlag, pdeShift);
[applyMax,AmaxShift] = shiftedLaplacianPreconditionerMaxwell2D(node, elem, bdFlag, k, opts);
```

The Maxwell wrapper assembles `curlcurl - (k^2+i*epsilon) mass - i eta boundaryMass` using the existing Nedelec utilities.
