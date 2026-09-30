function d = vap_results_dir()
%VAP_RESULTS_DIR Default output folder: <project root>/results.
d = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
end
