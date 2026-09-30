function s = vap_scenario_schedule(p)
%VAP_SCENARIO_SCHEDULE Shared drivers for controlled SOC/constant-H comparison.
validateattributes(p.scenario_t_end,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.scenario_dt,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.scenario_Q0,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(p.scenario_Q_end,{'numeric'},{'vector','real','finite','nonnegative','nonempty'});
validateattributes(p.scenario_SOC0,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(p.scenario_SOC1,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
s.t=unique([0:p.scenario_dt:p.scenario_t_end,p.scenario_t_end])';
s.tmin=s.t/60;
s.SOC=p.scenario_SOC0+(p.scenario_SOC1-p.scenario_SOC0)*s.t/p.scenario_t_end;
s.Q0=p.scenario_Q0; s.Q_end_list=p.scenario_Q_end(:)';
s.Qdot=s.Q0+(s.t/p.scenario_t_end)*(s.Q_end_list-s.Q0);
s.labels=cell(numel(s.Q_end_list),1);
for k=1:numel(s.labels)
    if s.Q_end_list(k)==s.Q0
        s.labels{k}=sprintf('Q = %g kW (const)',s.Q0/1e3);
    else
        s.labels{k}=sprintf('Q: %g -> %g kW',s.Q0/1e3,s.Q_end_list(k)/1e3);
    end
end
end
