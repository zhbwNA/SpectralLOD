function parts = ggglsPartitionOfUnity2D(parts, gridSize, delta)
% GGGLSPARTITIONOFUNITY2D  Attach Section 8 GGGLS smooth PoU weights.

if isscalar(gridSize)
    nx = gridSize;
    ny = 1;
else
    nx = gridSize(1);
    ny = gridSize(2);
end

for s = 1:numel(parts)
    ij = parts(s).ij;
    parts(s).weightFun = @(x,y) normalizedWeight(x, y, ij(1), ij(2), nx, ny, delta);
end
end


function w = normalizedWeight(x, y, i, j, nx, ny, delta)
den = zeros(size(x));
num = rawTensorWeight(x, y, i, j, nx, ny, delta);
for jj = 1:ny
    for ii = 1:nx
        den = den + rawTensorWeight(x, y, ii, jj, nx, ny, delta);
    end
end
w = zeros(size(x));
good = den > 0;
w(good) = num(good) ./ den(good);
end


function w = rawTensorWeight(x, y, i, j, nx, ny, delta)
wx = rawOneDimWeight(x, i, nx, delta);
if ny == 1
    wy = ones(size(y));
else
    wy = rawOneDimWeight(y, j, ny, delta);
end
w = wx .* wy;
end


function w = rawOneDimWeight(x, i, n, delta)
coreL = (i - 1) / n;
coreR = i / n;
a = max(0, coreL - delta);
b = min(1, coreR + delta);
eta = 0.3 * delta;
w = zeros(size(x));

if n == 1
    w(:) = 1;
elseif i == 1
    right = b - eta;
    mask = x < right;
    denom = max(right - x(mask), eps);
    w(mask) = exp(-(b - a) ./ (2 * denom));
elseif i == n
    left = a + eta;
    mask = x > left;
    denom = max(x(mask) - left, eps);
    w(mask) = exp(-(b - a) ./ (2 * denom));
else
    left = a + eta;
    right = b - eta;
    mask = x > left & x < right;
    denom = 4 * (x(mask) - left) .* (right - x(mask));
    w(mask) = exp(-((b - a)^2) ./ max(denom, eps));
end
end
