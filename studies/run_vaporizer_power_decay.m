function S = run_vaporizer_power_decay(p, outputDir, makePlots)
%RUN_VAPORIZER_POWER_DECAY Quasi-steady samples with liquid/vapor bookkeeping.
% The imposed fixed P,T are reservoir conditions, not a closed-tank transient.
if nargin<1 || isempty(p), p=vap_default_params(); end
if nargin<2 || isempty(outputDir), outputDir=vap_results_dir(); end
if nargin<3, makePlots=usejava('desktop'); end
p=vap_validate_params(p);
validateattributes(p.t_end,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.dt,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.Qdot0,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(p.Qdot_tau,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.H0,{'numeric'},{'scalar','real','finite','>=',0,'<=',p.D_tank});
validateattributes(p.deplete_level,{'logical','numeric'},{'scalar','binary'});
S.p=p; S.t=unique([0:p.dt:p.t_end,p.t_end])';
S.Qdot=p.Qdot0*exp(-S.t/p.Qdot_tau);
n=numel(S.t);
for name={'H','mdot','x_out','T_out','mass_L','pressure_residual','energy_residual'}
    S.(name{1})=nan(n,1);
end
S.ok=false(n,1); S.regime=repmat({'invalid'},n,1);
S.status=repmat({'not_evaluated'},n,1); S.message=repmat({''},n,1);
S.mass_vaporized=zeros(n,1); % cumulative transfer, not loss from total tank
S.vapor_transfer_rate=nan(n,1);
S.interval_vaporized=nan(n-1,1);
S.interval_active_fraction=nan(n-1,1);
th=vap_thermo(p.P,p.T_L,p.fluid);
S.T_sat=th.T_sat; S.rho_L=th.rho_L;
S.H(1)=p.H0;
S.mass_L(1)=th.rho_L*vap_height_to_volume(p.H0,p.D_tank,p.tank_length);
for k=1:n
    if ~isfinite(S.H(k))
        S.status{k}='upstream_failure';
        S.message{k}='Inventory is unknown after a prior solver failure.';
        if k<n, S.mass_vaporized(k+1)=NaN; end
        continue
    end
    r=vaporizer_thermosiphon(p.P,p.T_L,S.H(k),S.Qdot(k),p);
    S.ok(k)=r.ok; S.regime{k}=r.regime; S.status{k}=r.status; S.message{k}=r.message;
    if r.ok
        for name={'mdot','x_out','T_out','pressure_residual','energy_residual'}
            S.(name{1})(k)=r.(name{1});
        end
        if r.mdot==0, rate=0; else, rate=r.mdot*r.x_out; end
        S.vapor_transfer_rate(k)=rate;
    end
    if k==n, continue; end
    if ~p.deplete_level
        % Externally maintained level, independent of whether this solve fails.
        S.H(k+1)=S.H(k); S.mass_L(k+1)=S.mass_L(k);
        if r.ok
            dm=rate*(S.t(k+1)-S.t(k));
            S.interval_vaporized(k)=dm; S.interval_active_fraction(k)=1;
            S.mass_vaporized(k+1)=S.mass_vaporized(k)+dm;
        else
            S.mass_vaporized(k+1)=NaN;
        end
    elseif ~r.ok
        S.mass_vaporized(k+1)=NaN; % never advance with a fabricated zero flow
    else
        step=S.t(k+1)-S.t(k);
        [S.mass_L(k+1),dm]=vap_update_inventory(S.mass_L(k),r.mdot,r.x_out,step);
        S.interval_vaporized(k)=dm;
        if rate>0, S.interval_active_fraction(k)=dm/(rate*step); else, S.interval_active_fraction(k)=double(S.mass_L(k)>0); end
        S.mass_vaporized(k+1)=S.mass_vaporized(k)+dm;
        S.H(k+1)=vap_SOC_to_height(S.mass_L(k+1)/(th.rho_L*p.vol_L/1000),p.D_tank,p.tank_length,p.vol_L/1000);
    end
end
S.failed_count=nnz(~S.ok);
if S.failed_count>0, warning('vap:FailedCases','Power decay has %d failed/unevaluated samples.',S.failed_count); end
S.metadata=vap_result_metadata();
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,'power_decay_results.mat'),'-struct','S');
fprintf('Power decay: %d samples, %d failures; H %.6f -> %.6f m, vapor transfer %.4f kg\n', ...
    n,S.failed_count,S.H(1),S.H(end),S.mass_vaporized(end));
if makePlots, vap_plot_results(S,'power_decay',fullfile(outputDir,'figures')); end
end
