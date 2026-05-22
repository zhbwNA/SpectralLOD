function parts = ggglsSubdomains2D(node, elem, gridSize, delta, kappa)
% GGGLSSUBDOMAINS2D  Paper-style overlapping subdomains plus PML boxes.

if isscalar(gridSize)
    nx = gridSize;
    ny = 1;
else
    nx = gridSize(1);
    ny = gridSize(2);
end
nSub = nx * ny;
parts = repmat(struct('elemIdx', [], 'nodeIdx', [], 'coreBox', [], ...
    'pmlBox', [], 'ij', [], 'weightFun', []), nSub, 1);

cent = (node(elem(:,1), :) + node(elem(:,2), :) + node(elem(:,3), :)) / 3;
idx = 0;
for j = 1:ny
    for i = 1:nx
        idx = idx + 1;
        x0 = (i - 1) / nx; x1 = i / nx;
        y0 = (j - 1) / ny; y1 = j / ny;
        core = [max(0, x0 - delta), min(1, x1 + delta), ...
                max(0, y0 - delta), min(1, y1 + delta)];
        pbox = [core(1) - kappa, core(2) + kappa, core(3) - kappa, core(4) + kappa];
        eIdx = find(cent(:,1) >= pbox(1) - 1e-12 & cent(:,1) <= pbox(2) + 1e-12 & ...
                    cent(:,2) >= pbox(3) - 1e-12 & cent(:,2) <= pbox(4) + 1e-12);
        parts(idx).elemIdx = eIdx;
        parts(idx).nodeIdx = unique(elem(eIdx, :));
        parts(idx).coreBox = core;
        parts(idx).pmlBox = pbox;
        parts(idx).ij = [i, j];
    end
end
end
