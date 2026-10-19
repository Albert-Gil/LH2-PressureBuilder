%MAIN_VAPORIZER Single entry point for the LH2 vaporizer studies.
% Switch studies on/off in the RUN block, adjust options and parameters,
% then run this script. Every study uses the same parameter struct p.

%% Simulation switches
RUN.sweep        = true;   % height x power map              -> sweep_results.mat
RUN.power_decay  = true;   % exponential power decay         -> power_decay_results.mat
RUN.SOC          = true;   % prescribed SOC drawdown         -> SOC_scenario_results.mat
RUN.constH       = true;   % constant height, same powers    -> constH_scenario_results.mat
RUN.calibration  = false;  % one-point K_g fit; defaults and saved results unchanged
RUN.python_plots = false;  % batch PNGs via python/plot_vaporizer_results.py

%% Options
makePlots = true;          % MATLAB figures saved to <outputDir>/figures
outputDir = '';            % '' means <project root>/results
calibration_target_kg_min = 0.30*3.6; % [kg/min] scaled target = 1.08

%% Parameters (shared by every study)
projectRoot = fileparts(mfilename('fullpath'));
codeDirs = fullfile(projectRoot,{'model','geometry','studies','output'});
addpath(codeDirs{:});
p = vap_default_params();
% p.D = 0.025;             % example override; area is recomputed automatically

%% Run selected studies
if isempty(outputDir), outputDir = vap_results_dir(); end
results = struct();
if RUN.sweep,       results.sweep       = run_vaporizer_sweep(p,outputDir,makePlots); end
if RUN.power_decay, results.power_decay = run_vaporizer_power_decay(p,outputDir,makePlots); end
if RUN.SOC,         results.SOC         = run_vaporizer_SOC_scenario(p,outputDir,makePlots); end
if RUN.constH,      results.constH      = run_vaporizer_constH_scenario(p,outputDir,makePlots); end
if RUN.calibration
    [p_fit, fit] = calibrate_vaporizer_loss(calibration_target_kg_min,p);
    results.calibration = struct('p',p_fit,'report',fit);
end
if RUN.python_plots
    % Use the same Python interpreter that MATLAB uses for CoolProp.
    pe = pyenv;
    [status, cmdout] = system(sprintf('"%s" "%s" --directory "%s"', pe.Executable, ...
        fullfile(projectRoot,'python','plot_vaporizer_results.py'), outputDir));
    fprintf('%s', cmdout);
    if status ~= 0, warning('vap:PythonPlots','Python plotting failed (exit %d).',status); end
end
