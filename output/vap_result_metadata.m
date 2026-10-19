function metadata = vap_result_metadata()
%VAP_RESULT_METADATA Save enough provenance to distinguish model generations.
metadata.model_version='2.0-horizontal-heater';
metadata.schema_version=1;
metadata.generated_utc=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'));
metadata.matlab_version=version;
metadata.coolprop_version=char(py.getattr(py.importlib.import_module('CoolProp'),'__version__'));
metadata.geometry='Horizontal lumped heater; constant-enthalpy vertical return';
metadata.heat_definition='Net absorbed heat in the horizontal heater [W]';
% Git is optional (e.g. a downloaded source archive). Unknown is not clean.
metadata.git_commit='unknown';
metadata.git_dirty=NaN;
projectRoot=fileparts(fileparts(mfilename('fullpath')));
[status,commit]=system(sprintf('git -C "%s" rev-parse HEAD',projectRoot));
if status==0
    metadata.git_commit=strtrim(commit);
    % Exclude generated outputs: saving one study must not dirty the next
    % study's source provenance. Untracked source files are still included.
    [status,changes]=system(sprintf( ...
        'git -C "%s" status --porcelain --untracked-files=normal -- . ":(exclude)results"',projectRoot));
    if status==0, metadata.git_dirty=~isempty(strtrim(changes)); end
end
end
