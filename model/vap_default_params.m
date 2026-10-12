function p = vap_default_params()
%VAP_DEFAULT_PARAMS SI parameters for a horizontal heater and vertical return.
% Qdot is net heat absorbed in the bottom horizontal fin-and-plate heater.
% The bare return pipe is approximated as adiabatic, at constant pressure.
p.fluid = 'parahydrogen';
p.g = 9.81;
p.P = 8.273e5;              % [Pa] representative ullage/loop pressure
p.T_L = 25.0;               % [K] subcooled inlet
p.D_tank = 2.49;            % [m] horizontal cylindrical tank diameter
p.vol_L = 68000;            % [litres] geometric tank capacity
p.H_ref = 1.54;             % [m] approximately 65% volume fill

% Scaled design: 3.6x the flow of D = 0.020 m at the same outlet temperature
% and liquid heights, with all powers also scaled by 3.6.
p.D = 0.03592;             % [m] circular supply/return pipe ID
% A is derived from D by vap_validate_params, never independently specified.
p.L_liq = 1.5;             % [m] unheated liquid supply pipe
p.H_elev = [];             % [m] [] means return to tank top: D_tank
p.L_vap = [];              % [m] [] means straight vertical return: H_elev
p.roughness = 15e-6;       % [m] pipe wall roughness

% Effective losses lump fittings and unrepresented heater passage resistance,
% referenced to inlet (K_L) and outlet (K_g) velocities using the circular
% pipe area. Plate passage geometry is unknown; this is not a detailed heater
% model. Historical values at D = 0.020 m were 7.5 and 6.0; they are scaled
% by the fully turbulent friction factor ratio fT(D)/fT(0.020) = 0.875
% (Crane-style K = n*fT), assuming the heater flow area scales with power.
p.K_L = 6.565;
p.K_g = 5.252;
p.Qdot_ref = 21.6e3;       % [W] 3.6 x 6 kW reference point
p.STM_mdot_annotation = 1.1/60; % [kg/s] comparison only, not solver input

p.H_list = linspace(0.01, 2.0, 32);
p.Qdot_list = linspace(0.1, 10, 50)*1e3*3.6;
% Identical power schedules isolate the effect of height between scenarios.
p.scenario_t_end = 90*60;
p.scenario_dt = 30;
p.scenario_Q0 = 8e3*3.6;
p.scenario_Q_end = [8,7,6,5,4,3]*1e3*3.6;
p.scenario_SOC0 = 0.90;
p.scenario_SOC1 = 0.30;

p.t_end = 30*60;
p.dt = 30;
p.Qdot0 = 42e3;            % [W] start below the 250 K outlet cap at H0
p.Qdot_tau = 15*60;
p.H0 = 1.2;
p.deplete_level = true;

p.mdot_min = 1e-8;         % [kg/s] search limit, not a stagnation test
p.mdot_max = 0.5;          % [kg/s] initial search upper bound
p.mdot_hard_max = 5;       % [kg/s] maximum search upper bound
p.bracket_points = 48;
p.pressure_abs_tol = 1e-5; % [Pa] final momentum-balance acceptance
p.pressure_rel_tol = 1e-8;
p.T_out_max = 250;         % [K] outlet temperature limit; [] uses the fluid EOS maximum
end
