function S = run_vaporizer_constH_scenario(p, outputDir, makePlots)
%RUN_VAPORIZER_CONSTH_SCENARIO Fixed height with the same powers as the SOC run.
if nargin<1 || isempty(p), p=vap_default_params(); end
if nargin<2 || isempty(outputDir), outputDir=vap_results_dir(); end
if nargin<3, makePlots=usejava('desktop'); end
S=vap_run_scenario('constH',p,outputDir,makePlots);
end
