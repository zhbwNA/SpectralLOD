function [A, b, freeDof, bdDof, info] = assembleHelmholtzPML2D(node, elem, k, pml, f)
% ASSEMBLEHELMHOLTZPML2D  Assemble P1 Helmholtz operator with non-divergence PML.
%
%   [A,b,freeDof,bdDof] = ASSEMBLEHELMHOLTZPML2D(node,elem,k,pml,f)
%
%   The bilinear form is the weak form of -sum_l (1/s_l partial_l)^2 - k^2:
%       int d_l partial_l u partial_l v + int v beta_l partial_l u - k^2 int u v,
%   with d_l = s_l^(-2) and beta_l = -s_l'/s_l^3.
%
%   Homogeneous Dirichlet data are imposed by the caller on bdDof, usually
%   by solving A(freeDof,freeDof) u = b(freeDof).

if nargin < 5
    f = [];
end
if size(elem, 2) ~= 3
    error('assembleHelmholtzPML2D:p1Only', 'PML assembly currently supports P1 triangles only.');
end
if nargin < 4 || isempty(pml)
    pml = struct();
end
if ~isfield(pml, 'physicalBox') || isempty(pml.physicalBox)
    pml.physicalBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
if ~isfield(pml, 'pmlBox') || isempty(pml.pmlBox)
    pml.pmlBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
if ~isfield(pml, 'quadOrder') || isempty(pml.quadOrder)
    pml.quadOrder = 4;
end

N = size(node, 1);
coef = struct();
coef.d11 = @(x, y) pmlCoefField(x, y, k, pml, 'd11');
coef.d22 = @(x, y) pmlCoefField(x, y, k, pml, 'd22');
coef.beta1 = @(x, y) pmlCoefField(x, y, k, pml, 'beta1');
coef.beta2 = @(x, y) pmlCoefField(x, y, k, pml, 'beta2');

K = assembleNondivStiffness2D(node, elem, 1, coef, struct('quadOrder', pml.quadOrder));
M = assembleMass2D(node, elem, 1);
A = K - k^2 * M;

if nargout > 1
    b = assembleLoadVector(node, elem, pml.quadOrder, f);
end

bdDof = outerBoundaryNodes2D(node, pml.pmlBox);
freeDof = setdiff((1:N).', bdDof(:));

info = struct();
info.form = 'nondivergence';
info.physicalBox = pml.physicalBox;
info.pmlBox = pml.pmlBox;
info.dirichletDof = bdDof;
info.freeDof = freeDof;
end


function b = assembleLoadVector(node, elem, quadOrder, f)
N = size(node, 1);
if isempty(f)
    b = zeros(N, 1);
    return;
end

[lambda, weight] = quadtriangle(quadOrder);
x1 = node(elem(:,1), 1);  y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);  y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);  y3 = node(elem(:,3), 2);
area = abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1)) / 2;
rhs = zeros(size(elem, 1), 3);

for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1)*x1 + lq(2)*x2 + lq(3)*x3;
    yq = lq(1)*y1 + lq(2)*y2 + lq(3)*y3;
    wphys = 2 * weight(q) * area;
    fq = evalLoad(f, xq, yq);
    rhs = rhs + wphys .* fq .* lq;
end

b = accumarray(reshape(elem, [], 1), reshape(rhs, [], 1), [N, 1], @sum, 0);
end


function val = pmlCoefField(x, y, k, pml, name)
[d11, d22, beta1, beta2] = pmlNondivCoefficients2D(x, y, k, pml);
switch name
    case 'd11'
        val = d11;
    case 'd22'
        val = d22;
    case 'beta1'
        val = beta1;
    case 'beta2'
        val = beta2;
    otherwise
        error('assembleHelmholtzPML2D:coefName', 'Unknown PML coefficient field.');
end
end


function val = evalLoad(f, x, y)
if isnumeric(f)
    val = f * ones(size(x));
else
    val = f(x, y);
end
end


function bd = outerBoundaryNodes2D(node, box)
tol = 100 * eps(max(1, max(abs(box))));
onX = abs(node(:,1) - box(1)) <= tol | abs(node(:,1) - box(2)) <= tol;
onY = abs(node(:,2) - box(3)) <= tol | abs(node(:,2) - box(4)) <= tol;
bd = find(onX | onY);
end
