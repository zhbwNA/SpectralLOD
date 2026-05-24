function A = assembleNondivStiffness2D(node, elem, degree, coef, opts)
% ASSEMBLENONDIVSTIFFNESS2D  Assemble int D grad u.grad v + int v beta.grad u.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end
if nargin < 5 || isempty(opts), opts = struct(); end

[diffCoef, betaCoef] = splitNondivCoefficients(coef);
A = assembleDiffusion2D(node, elem, degree, diffCoef, opts) + ...
    assembleAdvection2D(node, elem, degree, betaCoef, opts);
end


function [diffCoef, betaCoef] = splitNondivCoefficients(coef)
if isnumeric(coef) || isa(coef, 'function_handle')
    diffCoef = coef;
    betaCoef = [0, 0];
    return;
end

if ~isstruct(coef)
    error('assembleNondivStiffness2D:coef', 'Coefficient must be scalar, function handle, or struct.');
end

diffCoef = struct();
for name = {'d11', 'd12', 'd21', 'd22'}
    key = name{1};
    if isfield(coef, key)
        diffCoef.(key) = coef.(key);
    end
end

betaCoef = struct();
for name = {'beta1', 'beta2'}
    key = name{1};
    if isfield(coef, key)
        betaCoef.(key) = coef.(key);
    end
end
end
