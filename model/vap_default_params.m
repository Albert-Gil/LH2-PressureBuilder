function p = vap_default_params()
%VAP_DEFAULT_PARAMS SI parameters for a horizontal heater and vertical return.
% Qdot is net heat absorbed in the bottom horizontal fin-and-plate heater.
% The bare return pipe is approximated as adiabatic, at constant pressure.
p.fluid = 'parahydrogen';
p.g = 9.81;
p.P = 8.273e5;              % [Pa] representative ullage/loop pressure
p.T_L = 27.0;               % [K] subcooled inlet
p.D_tank = 2.49;            % [m] horizontal cylindrical tank diameter
p.vol_L = 68000;            % [litres] geometric tank capacity
p.H_ref = 1.54;             % [m] approximately 65% volume fill

p.D = 0.020;               % [m] circular supply/return pipe ID
% A is derived from D by vap_validate_params, never independently specified.
p.L_liq = 1.5;             % [m] unheated liquid supply pipe
p.H_elev = [];             % [m] [] means return to tank top: D_tank
p.L_vap = [];              % [m] [] means straight vertical return: H_elev
p.roughness = 15e-6;       % [m] pipe wall roughness

% Historical effective losses, retained, NOT refitted after repair.
% Lump fittings and unrepresented heater passage resistance, referenced
% to inlet (K_L) and outlet (K_g) velocities using the circular pipe area.
% Plate passage geometry is unknown; this is not a detailed heater model.
p.K_L = 7.5;
p.K_g = 6.0;
p.Qdot_ref = 6e3;          % [W] documented approximate field comparison
p.STM_mdot_annotation = 0.30/60; % [kg/s] comparison only, not solver input

p.H_list = linspace(0.01, 2.0, 32);
p.Qdot_list = linspace(0.1, 10, 50)*1e3;
% Identical power schedules isolate the effect of height between scenarios.
p.scenario_t_end = 90*60;
p.scenario_dt = 30;
p.scenario_Q0 = 8e3;
p.scenario_Q_end = [8,7,6,5,4,3]*1e3;
p.scenario_SOC0 = 0.90;
p.scenario_SOC1 = 0.30;

p.t_end = 30*60;
p.dt = 30;
p.Qdot0 = 12e3;
p.Qdot_tau = 15*60;
p.H0 = 1.2;
p.deplete_level = true;

p.mdot_min = 1e-8;         % [kg/s] search limit, not a stagnation test
p.mdot_max = 0.5;          % [kg/s] initial search upper bound
p.mdot_hard_max = 5;       % [kg/s] maximum search upper bound
p.bracket_points = 48;
p.pressure_abs_tol = 1e-5; % [Pa] final momentum-balance acceptance
p.pressure_rel_tol = 1e-8;
p.T_out_max = [];          % [K] [] uses the fluid EOS maximum temperature
end
