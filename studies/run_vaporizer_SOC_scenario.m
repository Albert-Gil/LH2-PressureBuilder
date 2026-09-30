function S = run_vaporizer_SOC_scenario(p, outputDir, makePlots)
%RUN_VAPORIZER_SOC_SCENARIO Prescribed SOC using the shared power schedules.
if nargin<1 || isempty(p), p=vap_default_params(); end
if nargin<2 || isempty(outputDir), outputDir=vap_results_dir(); end
if nargin<3, makePlots=usejava('desktop'); end
S=vap_run_scenario('SOC',p,outputDir,makePlots);
end
