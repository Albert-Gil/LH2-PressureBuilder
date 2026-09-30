function S = run_vaporizer_sweep(p, outputDir, makePlots)
%RUN_VAPORIZER_SWEEP Evaluate H x absorbed power, retaining failure diagnostics.
if nargin<1 || isempty(p), p=vap_default_params(); end
if nargin<2 || isempty(outputDir), outputDir=vap_results_dir(); end
if nargin<3, makePlots=usejava('desktop'); end
p=vap_validate_params(p);
S.p=p; S.H_list=p.H_list; S.Q_list=p.Qdot_list;
[Qg,Hg]=meshgrid(S.Q_list,S.H_list);
fprintf('Sweep: %d heights x %d powers\n',numel(S.H_list),numel(S.Q_list));
R=vap_run_cases(p.P,p.T_L,Hg,Qg,p);
S.mdot_map=R.mdot; S.x_map=R.x_out; S.Tout_map=R.T_out;
S.ok=R.ok; S.regime_map=R.regime; S.status=R.status; S.message=R.message;
S.pressure_residual=R.pressure_residual; S.energy_residual=R.energy_residual;
S.failed_count=R.failed_count;
S.sp=vaporizer_thermosiphon(p.P,p.T_L,p.H_ref,p.Qdot_ref,p);
[~,S.iHref]=min(abs(S.H_list-p.H_ref));
[~,S.iQref]=min(abs(S.Q_list-p.Qdot_ref));
S.mdot_vs_Q=S.mdot_map(S.iHref,:); S.x_vs_Q=S.x_map(S.iHref,:);
S.Tout_vs_Q=S.Tout_map(S.iHref,:);
S.mdot_vs_H=S.mdot_map(:,S.iQref); S.Tout_vs_H=S.Tout_map(:,S.iQref);
th=vap_thermo(p.P,p.T_L,p.fluid);
S.ceiling_vs_Q=S.Q_list/(th.h_g-th.h_L);
S.metadata=vap_result_metadata();
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,'sweep_results.mat'),'-struct','S');
fprintf('Reference: ok=%d, %.6f kg/min, %.3f K; sweep failures %d/%d\n', ...
    S.sp.ok,60*S.sp.mdot,S.sp.T_out,S.failed_count,numel(Qg));
if makePlots, vap_plot_results(S,'sweep',fullfile(outputDir,'figures')); end
end
