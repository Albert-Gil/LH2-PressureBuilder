function th = vap_thermo(P, T_L, fluid)
%VAP_THERMO Consistent inlet and saturation properties, cached per state.
persistent cacheKey cacheValue
validateattributes(P,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(T_L,{'numeric'},{'scalar','real','finite','positive'});
key = {P,T_L,char(fluid)};
if isequal(key,cacheKey), th = cacheValue; return; end
cp = @(o,s,v) vap_coolprop(o,'P',P,s,v,fluid);
th.P = P; th.T_L = T_L; th.fluid = char(fluid);
th.T_min = cp('Tmin','T',T_L);
if T_L < th.T_min
    error('vap:InletTemperature','T_L=%g K is below the fluid EOS minimum (%g K).',T_L,th.T_min);
end
th.T_sat = cp('T','Q',0);
if T_L > th.T_sat + 1e-7
    error('vap:InletPhase','T_L=%g K exceeds saturation (%g K); a liquid inlet is required.',T_L,th.T_sat);
end
th.h_f = cp('H','Q',0); th.h_g = cp('H','Q',1);
th.rho_f = cp('D','Q',0); th.rho_g = cp('D','Q',1);
th.mu_f = cp('V','Q',0); th.mu_g = cp('V','Q',1);
if abs(T_L-th.T_sat) <= 1e-7
    th.h_L = th.h_f; th.rho_L = th.rho_f; th.mu_L = th.mu_f;
else
    th.h_L = cp('H','T',T_L);
    th.rho_L = cp('D','T',T_L);
    th.mu_L = cp('V','T',T_L);
end
th.h_fg = th.h_g-th.h_f;
th.T_max = cp('Tmax','T',T_L);
th.h_max = cp('H','T',th.T_max);
if th.h_fg <= 0 || th.h_L > th.h_f+1e-4
    error('vap:InletPhase','No valid subcritical liquid/vapor state at the requested inlet.');
end
cacheKey = key; cacheValue = th;
end
