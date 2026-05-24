function A = assembleNondivCurlCurl2D(node, elem, coef, opts)
% ASSEMBLENONDIVCURLCURL2D  Assemble int coef curl w_i curl w_j for NE_1 2D.

if nargin < 3 || isempty(coef), coef = 1; end
if nargin < 4 || isempty(opts), opts = struct(); end
A = assembleWeightedCurlCurl2D(node, elem, coef, opts);
end
