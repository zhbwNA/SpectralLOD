function parts = smoothPartitionOfUnity2D(parts, bbox, gridSize, overlap)
% SMOOTHPARTITIONOFUNITY2D  Attach C-infinity bump POU weights to parts.
%
%   parts = SMOOTHPARTITIONOFUNITY2D(parts,bbox,[nx ny],overlap)
%
%   The raw weights are tensor products of smooth ramps over the overlap
%   regions. Consumers should normalize overlapping weights nodally so that
%   the final partition sums to one.

xmin = bbox(1); xmax = bbox(2);
ymin = bbox(3); ymax = bbox(4);
nx = gridSize(1); ny = gridSize(2);
Hx = (xmax - xmin) / nx;
Hy = (ymax - ymin) / ny;

for j = 1:ny
    for i = 1:nx
        s = (j - 1) * nx + i;
        xL = xmin + (i - 1) * Hx;
        xR = xmin + i * Hx;
        yB = ymin + (j - 1) * Hy;
        yT = ymin + j * Hy;

        parts(s).coreBox = [xL, xR, yB, yT];
        parts(s).extendedBox = [max(xmin, xL - overlap), min(xmax, xR + overlap), ...
                                max(ymin, yB - overlap), min(ymax, yT + overlap)];
        parts(s).overlap = overlap;
        parts(s).weightFun = @(x,y) smoothBoxWeight(x, y, ...
            xL, xR, yB, yT, xmin, xmax, ymin, ymax, overlap);
    end
end
end


function w = smoothBoxWeight(x, y, xL, xR, yB, yT, xmin, xmax, ymin, ymax, overlap)
w = smoothOneDimWeight(x, xL, xR, xmin, xmax, overlap) .* ...
    smoothOneDimWeight(y, yB, yT, ymin, ymax, overlap);
end


function w = smoothOneDimWeight(x, xL, xR, xmin, xmax, overlap)
if overlap <= 0
    w = double(x >= xL - 1e-12 & x <= xR + 1e-12);
    return;
end

w = ones(size(x));
if xL > xmin + 1e-12
    t = (x - (xL - overlap)) / (2 * overlap);
    w = min(w, smoothStep(t));
    w(x < xL - overlap) = 0;
end
if xR < xmax - 1e-12
    t = ((xR + overlap) - x) / (2 * overlap);
    w = min(w, smoothStep(t));
    w(x > xR + overlap) = 0;
end
w = max(0, min(1, w));
end


function y = smoothStep(t)
y = zeros(size(t));
y(t >= 1) = 1;
mid = t > 0 & t < 1;
if any(mid(:))
    tm = t(mid);
    a = exp(-1 ./ tm);
    b = exp(-1 ./ (1 - tm));
    y(mid) = a ./ (a + b);
end
end
