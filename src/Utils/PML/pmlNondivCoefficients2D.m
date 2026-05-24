function [d11, d22, beta1, beta2, s1, s2] = pmlNondivCoefficients2D(x, y, k, pml)
% PMLNONDIVCOEFFICIENTS2D  Coefficients for sum_l (1/s_l d_l)^2.

if nargin < 4 || isempty(pml)
    pml = struct();
end
if ~isfield(pml, 'physicalBox')
    error('pmlNondivCoefficients2D:missingBox', 'pml.physicalBox is required.');
end
if ~isfield(pml, 'sigmaMax') || isempty(pml.sigmaMax)
    pml.sigmaMax = k;
end
if ~isfield(pml, 'sigmaOrder') || isempty(pml.sigmaOrder)
    pml.sigmaOrder = 2;
end

box = pml.physicalBox(:).';
xmin = box(1); xmax = box(2);
ymin = box(3); ymax = box(4);

if isfield(pml, 'pmlBox') && ~isempty(pml.pmlBox)
    pbox = pml.pmlBox(:).';
else
    pbox = [min(x(:)), max(x(:)), min(y(:)), max(y(:))];
end

txL = max(xmin - pbox(1), eps);
txR = max(pbox(2) - xmax, eps);
tyB = max(ymin - pbox(3), eps);
tyT = max(pbox(4) - ymax, eps);

[sigma1, dsigma1] = oneDimSigmaDerivative(x, xmin, xmax, txL, txR, ...
    pml.sigmaMax, pml.sigmaOrder);
[sigma2, dsigma2] = oneDimSigmaDerivative(y, ymin, ymax, tyB, tyT, ...
    pml.sigmaMax, pml.sigmaOrder);

s1 = 1 + 1i * sigma1 / k;
s2 = 1 + 1i * sigma2 / k;
ds1 = 1i * dsigma1 / k;
ds2 = 1i * dsigma2 / k;

d11 = 1 ./ (s1.^2);
d22 = 1 ./ (s2.^2);
beta1 = -ds1 ./ (s1.^3);
beta2 = -ds2 ./ (s2.^3);
end


function [sigma, dsigma] = oneDimSigmaDerivative(x, xmin, xmax, tLeft, tRight, sigmaMax, sigmaOrder)
sigma = zeros(size(x));
dsigma = zeros(size(x));

left = x < xmin;
if any(left(:))
    rRaw = (xmin - x(left)) ./ tLeft;
    active = rRaw > 0 & rRaw < 1;
    r = min(1, max(0, rRaw));
    sigma(left) = sigmaMax * r.^sigmaOrder;
    tmp = zeros(size(rRaw));
    if sigmaOrder > 0
        tmp(active) = -sigmaMax * sigmaOrder * rRaw(active).^(sigmaOrder - 1) ./ tLeft;
    end
    dsigma(left) = tmp;
end

right = x > xmax;
if any(right(:))
    rRaw = (x(right) - xmax) ./ tRight;
    active = rRaw > 0 & rRaw < 1;
    r = min(1, max(0, rRaw));
    sigma(right) = sigmaMax * r.^sigmaOrder;
    tmp = zeros(size(rRaw));
    if sigmaOrder > 0
        tmp(active) = sigmaMax * sigmaOrder * rRaw(active).^(sigmaOrder - 1) ./ tRight;
    end
    dsigma(right) = tmp;
end
end
