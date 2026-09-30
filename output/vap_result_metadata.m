function metadata = vap_result_metadata()
%VAP_RESULT_METADATA Save enough provenance to distinguish model generations.
metadata.model_version='2.0-horizontal-heater';
metadata.generated_utc=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'));
metadata.matlab_version=version;
metadata.coolprop_version=char(py.getattr(py.importlib.import_module('CoolProp'),'__version__'));
metadata.geometry='Horizontal lumped heater; constant-enthalpy vertical return';
metadata.heat_definition='Net absorbed heat in the horizontal heater [W]';
end
