function val = evalPDECoefficient(coef, x, y, z, kval)
% EVALPDECOEFFICIENT  Evaluate scalar/nodal/function PDE coefficient data.

if nargin < 4, z = []; end
if nargin < 5, kval = []; end

if isnumeric(coef)
    if isscalar(coef)
        val = coef * ones(size(x));
    elseif numel(coef) == numel(x)
        val = reshape(coef, size(x));
    else
        error('evalPDECoefficient:size', 'Numeric coefficient size is incompatible with query points.');
    end
    return;
end

if isa(coef, 'function_handle')
    if isempty(z)
        try
            val = coef(x, y, kval);
        catch
            try
                val = coef(x, y);
            catch
                val = coef(kval);
            end
        end
    else
        try
            val = coef(x, y, z, kval);
        catch
            try
                val = coef(x, y, z);
            catch
                val = coef(kval);
            end
        end
    end
    if isscalar(val), val = val * ones(size(x)); end
    val = reshape(val, size(x));
    return;
end

error('evalPDECoefficient:type', 'Unsupported coefficient type.');
end
