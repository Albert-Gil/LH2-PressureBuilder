function st = vap_flow_state(mdot, H, Qdot, p, th)
%VAP_FLOW_STATE Trial flow in a horizontal heater / adiabatic vertical pipe.
% p must be normalized by vap_validate_params. At constant P the return pipe
% has uniform enthalpy, density and viscosity equal to the heater exit.
validateattributes(mdot,{'numeric'},{'scalar','real','finite','positive'});
h = th.h_L+Qdot/mdot;
if h < th.h_f
    st.x_out = 0; st.regime = 'liquid';
    st.T_out = vap_coolprop('T','P',th.P,'H',h,th.fluid);
    st.rho_out = vap_coolprop('D','P',th.P,'H',h,th.fluid);
    st.mu_out = vap_coolprop('V','P',th.P,'H',h,th.fluid);
elseif h <= th.h_g
    st.x_out = (h-th.h_f)/th.h_fg;
    st.regime = 'two_phase'; st.T_out = th.T_sat;
    st.rho_out = 1/((1-st.x_out)/th.rho_f+st.x_out/th.rho_g);
    % Explicit homogeneous no-slip viscosity closure, continuous at x=0,1.
    st.mu_out = 1/((1-st.x_out)/th.mu_f+st.x_out/th.mu_g);
else
    st.x_out = 1; st.regime = 'dry_vapor';
    st.T_out = vap_coolprop('T','P',th.P,'H',h,th.fluid);
    st.rho_out = vap_coolprop('D','P',th.P,'H',h,th.fluid);
    st.mu_out = vap_coolprop('V','P',th.P,'H',h,th.fluid);
end
if any(~isfinite([st.T_out,st.rho_out,st.mu_out])) || st.rho_out<=0 || st.mu_out<=0
    error('vap:InvalidState','Trial state has invalid thermodynamic properties.');
end
st.mdot = mdot; st.h_out = h;
st.rho_m = st.rho_out; % compatibility: actual uniform return-pipe density
st.Re_L = mdot*p.D/(p.A*th.mu_L);
st.Re_g = mdot*p.D/(p.A*st.mu_out);
st.f_L = friction_factor(st.Re_L,p.roughness/p.D);
st.f_g = friction_factor(st.Re_g,p.roughness/p.D);
G2 = (mdot/p.A)^2;
st.dP_drive = p.g*(th.rho_L*H-st.rho_out*p.H_elev);
st.dP_friction = G2/2*(st.f_L*p.L_liq/(p.D*th.rho_L) + st.f_g*p.L_vap/(p.D*st.rho_out));
st.dP_minor = G2/2*(p.K_L/th.rho_L+p.K_g/st.rho_out);
st.dP_acc = G2*(1/st.rho_out-1/th.rho_L);
st.dP_loss = st.dP_friction+st.dP_minor+st.dP_acc;
st.pressure_residual = st.dP_drive-st.dP_loss;
st.energy_residual = mdot*(h-th.h_L)-Qdot;
st.dT_superheat = max(st.T_out-th.T_sat,0);
end

function f = friction_factor(Re,rel)
% Darcy factor; smoothly blend laminar and Haaland through transition.
laminar = 64/Re;
if Re <= 2300, f = laminar; return; end
turbulent = (-1.8*log10((rel/3.7)^1.11+6.9/Re))^(-2);
if Re >= 4000, f = turbulent; return; end
s = (Re-2300)/1700; w = s*s*(3-2*s);
f = (1-w)*laminar+w*turbulent;
end
