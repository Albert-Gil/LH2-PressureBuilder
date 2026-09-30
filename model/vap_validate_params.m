function p = vap_validate_params(p)
%VAP_VALIDATE_PARAMS Resolve derived geometry and reject inconsistent inputs.
if ~isstruct(p) || ~isscalar(p), error('vap:Parameters','p must be a scalar struct.'); end
defaults = vap_default_params();
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(p,names{k}), p.(names{k}) = defaults.(names{k}); end
end
if ~(ischar(p.fluid) && isrow(p.fluid)) && ~(isstring(p.fluid) && isscalar(p.fluid))
    error('vap:Parameters','fluid must be a scalar name.');
end
p.fluid = char(p.fluid);
positive = {'g','D','D_tank','vol_L','mdot_min','mdot_max','mdot_hard_max', ...
    'pressure_abs_tol','pressure_rel_tol'};
for k = 1:numel(positive)
    validateattributes(p.(positive{k}),{'numeric'},{'scalar','real','finite','positive'},'',positive{k});
end
if isempty(p.H_elev), p.H_elev = p.D_tank; end
if isempty(p.L_vap), p.L_vap = p.H_elev; end
nonnegative = {'H_elev','L_liq','L_vap','K_L','K_g','roughness'};
for k = 1:numel(nonnegative)
    validateattributes(p.(nonnegative{k}),{'numeric'},{'scalar','real','finite','nonnegative'},'',nonnegative{k});
end
validateattributes(p.bracket_points,{'numeric'},{'scalar','integer','>=',8,'<=',10000});
if ~isempty(p.T_out_max)
    validateattributes(p.T_out_max,{'numeric'},{'scalar','real','finite','positive'});
end
if p.L_vap < p.H_elev || p.roughness >= p.D
    error('vap:Geometry','Return length must cover its elevation and roughness must be below diameter.');
end
if p.mdot_min >= p.mdot_max || p.mdot_max > p.mdot_hard_max
    error('vap:Bounds','Require 0 < mdot_min < mdot_max <= mdot_hard_max.');
end
% Ignore stale derived fields in legacy structs.
p.A = pi*p.D^2/4;
p.tank_length = 4*(p.vol_L/1000)/(pi*p.D_tank^2);
end
