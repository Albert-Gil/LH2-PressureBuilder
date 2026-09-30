function [candidate, report] = calibrate_vaporizer_loss(target_kg_min, p)
%CALIBRATE_VAPORIZER_LOSS One-point effective K_g fit; never changes defaults.
% Uses the stated approximate flow as calibration, NOT independent validation.
if nargin<1, target_kg_min=0.30; end
if nargin<2 || isempty(p), p=vap_default_params(); end
validateattributes(target_kg_min,{'numeric'},{'scalar','real','finite','positive'});
candidate=vap_validate_params(p); th=vap_thermo(p.P,p.T_L,p.fluid);
m=target_kg_min/60;
baseline=vap_flow_state(m,p.H_ref,p.Qdot_ref,candidate,th);
scale=(m/candidate.A)^2/(2*baseline.rho_out);
newK=candidate.K_g+baseline.pressure_residual/scale;
if ~isfinite(newK) || newK<0
    error('vap:Calibration','Target requires a negative/nonfinite K_g; review geometry and measurements.');
end
candidate.K_g=newK;
report.original_K_g=p.K_g; report.fitted_K_g=newK;
report.target_kg_min=target_kg_min;
report.solution=vaporizer_thermosiphon(p.P,p.T_L,p.H_ref,p.Qdot_ref,candidate);
if ~report.solution.ok || abs(report.solution.mdot-m)>max(1e-10,1e-6*m)
    error('vap:Calibration','Fitted state failed independent root verification.');
end
fprintf('One-point effective K_g: %.6g -> %.6g; %.6f kg/min, %.3f K. Defaults unchanged.\n', ...
    p.K_g,newK,report.solution.mdot*60,report.solution.T_out);
end
