function val = vap_coolprop(output, spec1, val1, spec2, val2, fluid)
%VAP_COOLPROP Strict SI query. Errors retain their cause and requested state.
if nargin < 6 || isempty(fluid), fluid = 'parahydrogen'; end
validateattributes(val1, {'numeric'}, {'real','finite','scalar'});
validateattributes(val2, {'numeric'}, {'real','finite','scalar'});
try
    val = double(py.CoolProp.CoolProp.PropsSI(output, spec1, val1, spec2, val2, fluid));
    if ~isscalar(val) || ~isfinite(val)
        error('vap:NonfiniteProperty', 'CoolProp returned a nonfinite scalar.');
    end
catch ME
    problem = MException('vap:PropertyFailure', ...
        'CoolProp %s(%s=%g, %s=%g, %s): %s', ...
        output, spec1, val1, spec2, val2, fluid, ME.message);
    problem = addCause(problem, ME);
    throw(problem);
end
end
