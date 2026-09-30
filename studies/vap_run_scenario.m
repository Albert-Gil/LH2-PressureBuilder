function S = vap_run_scenario(kind,p,outputDir,makePlots)
%VAP_RUN_SCENARIO Shared scenario implementation, including failure records.
p=vap_validate_params(p);
S=vap_scenario_schedule(p); S.p=p;
S.D=p.D_tank; S.Length=p.tank_length; S.vol_m3=p.vol_L/1000;
S.SOC0=p.scenario_SOC0; S.SOC1=p.scenario_SOC1;
if strcmp(kind,'SOC')
    S.H=vap_SOC_to_height(S.SOC,S.D,S.Length,S.vol_m3);
elseif strcmp(kind,'constH')
    S.H_const=p.H_ref; S.H=repmat(p.H_ref,size(S.t));
    S=rmfield(S,{'SOC','SOC0','SOC1'});
else
    error('vap:Scenario','Unknown scenario kind.');
end
R=vap_run_cases(p.P,p.T_L,S.H,S.Qdot,p);
names=fieldnames(R);
for k=1:numel(names), S.(names{k})=R.(names{k}); end
S.metadata=vap_result_metadata();
S.paired_count=zeros(size(S.labels));
S.mean_mdot=nan(size(S.labels)); S.mean_T_out=nan(size(S.labels));
for k=1:numel(S.labels)
    % Both averages use exactly the same valid, defined outlet-state samples.
    paired=S.ok(:,k)&isfinite(S.mdot(:,k))&isfinite(S.T_out(:,k));
    S.paired_count(k)=nnz(paired);
    if any(paired)
        S.mean_mdot(k)=mean(S.mdot(paired,k));
        S.mean_T_out(k)=mean(S.T_out(paired,k));
    end
    fprintf('%s | %s: paired mean %.4f kg/min, %.2f K (%d/%d samples), failures=%d\n', ...
        kind,S.labels{k},60*S.mean_mdot(k),S.mean_T_out(k),nnz(paired),numel(paired),nnz(~S.ok(:,k)));
end
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,[kind '_scenario_results.mat']),'-struct','S');
if makePlots, vap_plot_results(S,kind,fullfile(outputDir,'figures')); end
end
