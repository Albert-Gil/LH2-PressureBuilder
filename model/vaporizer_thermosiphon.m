function out = vaporizer_thermosiphon(P, T_L, H, Qdot, p)
%VAPORIZER_THERMOSIPHON Quasi-steady bottom-heater / vertical-return ROM.
% Inputs: P [Pa], liquid T_L [K], bottom takeoff head H [m], absorbed Qdot [W].
% Positive flow requires finite properties and a verified pressure balance.
% Failed/ambiguous solves have ok=false and mdot=NaN. No-flow success is
% reserved for empty inventory or zero heat without sufficient head.
out = empty_out(); out.P=P; out.T_L=T_L; out.H=H; out.Qdot=Qdot;
if nargin<5 || isempty(p), p=vap_default_params(); end
try
    validateattributes(P,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(T_L,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(H,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(Qdot,{'numeric'},{'scalar','real','finite','nonnegative'});
    p=vap_validate_params(p);
    if H>p.D_tank, error('vap:Height','H exceeds the tank diameter.'); end
catch ME
    out.status='invalid_input'; out.message=ME.message; return
end
try
    th=vap_thermo(P,T_L,p.fluid);
    out.rho_L=th.rho_L; out.mu_L=th.mu_L; out.h_L=th.h_L;
    out.h_f=th.h_f; out.h_g=th.h_g; out.h_fg=th.h_fg; out.T_sat=th.T_sat;
    out.mdot_Q_ceiling=Qdot/(th.h_g-th.h_L);
    out.A=p.A; out.H_elev=p.H_elev;
    if H==0
        out=no_flow(out,'empty_tank','No liquid inventory at the takeoff.'); return
    end
    if Qdot==0 && H<=p.H_elev
        out=no_flow(out,'zero_heat','Zero absorbed heat and insufficient liquid head.'); return
    end
    Tmax=th.T_max;
    if ~isempty(p.T_out_max), Tmax=min(Tmax,p.T_out_max); end
    if Tmax<=th.T_sat, error('vap:TemperatureLimit','T_out_max must exceed saturation temperature.'); end
    hmax=th.h_max;
    if Tmax<th.T_max, hmax=vap_coolprop('H','P',P,'T',Tmax,p.fluid); end
    % Thermodynamic limits bound the search; they never imply stagnation.
    lo=max(p.mdot_min,Qdot/(hmax-th.h_L)*(1+1e-9));
    hi=p.mdot_max;
    out.T_out_limit=Tmax;
    if lo>=p.mdot_hard_max
        out.status='no_root_in_domain'; out.message='Required enthalpy exceeds the temperature/search limits.'; return
    end
    hi=max(hi,min(p.mdot_hard_max,2*lo));
    residual=@(m) get_residual(m,H,Qdot,p,th);
    while residual(hi)>0 && hi<p.mdot_hard_max
        hi=min(2*hi,p.mdot_hard_max);
    end
    grid=logspace(log10(lo),log10(hi),p.bracket_points);
    % Explicit phase boundaries also resolve narrow near-dryout brackets.
    special=Qdot./[th.h_g-th.h_L,th.h_f-th.h_L];
    special=special(isfinite(special)&special>lo&special<hi);
    grid=unique([lo,grid,special,hi]);
    values=arrayfun(residual,grid);
    rootsFound=grid(values==0);
    options=optimset('Display','off','TolX',1e-12);
    for k=1:numel(grid)-1
        if sign(values(k))*sign(values(k+1))<0
            [m,~,exitflag]=fzero(residual,[grid(k),grid(k+1)],options);
            if exitflag<=0, error('vap:Convergence','Root search did not converge.'); end
            st=vap_flow_state(m,H,Qdot,p,th);
            tol=p.pressure_abs_tol+p.pressure_rel_tol*max(abs([st.dP_drive,st.dP_loss]));
            if abs(st.pressure_residual)>tol
                error('vap:Residual','Root candidate fails momentum balance: %g Pa.',st.pressure_residual);
            end
            rootsFound(end+1)=m; %#ok<AGROW>
        end
    end
    rootsFound=sort(rootsFound);
    if numel(rootsFound)>1
        rootsFound=rootsFound([true,diff(rootsFound)>1e-7*max(rootsFound(2:end),1e-12)]);
    end
    out.root_candidates=rootsFound; out.root_count=numel(rootsFound);
    out.search_bounds=[lo,hi];
    if isempty(rootsFound)
        out.status='no_root_in_domain';
        out.message='No pressure-balanced root in the valid property/flow domain; not a confirmed no-flow state.'; return
    elseif numel(rootsFound)>1
        out.status='multiple_roots'; out.message='Multiple steady roots; no stability/branch model selects one.'; return
    end
    st=vap_flow_state(rootsFound(1),H,Qdot,p,th);
    tol=p.pressure_abs_tol+p.pressure_rel_tol*max(abs([st.dP_drive,st.dP_loss]));
    if abs(st.pressure_residual)>tol || st.T_out>Tmax+1e-6 || ...
            abs(st.energy_residual)>1e-8*max(Qdot,1)
        error('vap:Residual','Final solution fails pressure, energy or temperature acceptance.');
    end
    names=fieldnames(st);
    for k=1:numel(names), out.(names{k})=st.(names{k}); end
    out.ok=true; out.status='solved'; out.message='ok';
catch ME
    out.ok=false; out.mdot=NaN; out.regime='invalid';
    out.status='calculation_failed'; out.message=[ME.identifier ': ' ME.message];
end
end

function r=get_residual(m,H,Qdot,p,th)
st=vap_flow_state(m,H,Qdot,p,th); r=st.pressure_residual;
end

function out=no_flow(out,status,message)
out.ok=true; out.mdot=0; out.regime='no_flow'; out.status=status; out.message=message;
% No steady outlet state or absorbed-heat balance is assigned at zero flow.
end

function out=empty_out()
out=struct('ok',false,'status','invalid_input','message','','regime','invalid', ...
    'P',NaN,'T_L',NaN,'H',NaN,'Qdot',NaN,'mdot',NaN,'mdot_Q_ceiling',NaN, ...
    'x_out',NaN,'T_out',NaN,'h_out',NaN,'h_L',NaN,'h_f',NaN,'h_g',NaN,'h_fg',NaN, ...
    'rho_L',NaN,'mu_L',NaN,'rho_out',NaN,'mu_out',NaN,'rho_m',NaN,'T_sat',NaN, ...
    'dP_drive',NaN,'dP_loss',NaN,'dP_acc',NaN,'dP_friction',NaN,'dP_minor',NaN, ...
    'pressure_residual',NaN,'energy_residual',NaN,'Re_L',NaN,'Re_g',NaN, ...
    'f_L',NaN,'f_g',NaN,'dT_superheat',NaN,'A',NaN,'H_elev',NaN, ...
    'root_count',0,'root_candidates',[],'search_bounds',[],'T_out_limit',NaN);
end
