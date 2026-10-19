function tests = test_vaporizer
%TEST_VAPORIZER MATLAB regressions; property tests require Python CoolProp.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath=path;
addpath(fullfile(root,'model'),fullfile(root,'geometry'), ...
    fullfile(root,'studies'),fullfile(root,'output'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testGeometryRoundTrip(testCase)
p=vap_validate_params(vap_default_params());
fractions=[0,1e-12,1e-8,.01,.25,.5,.75,.99,1-1e-8,1-1e-12,1];
h=vap_SOC_to_height(fractions,p.D_tank,p.tank_length,p.vol_L/1000);
v=vap_height_to_volume(h,p.D_tank,p.tank_length);
verifyEqual(testCase,v/(p.vol_L/1000),fractions,'AbsTol',1e-12);
verifyEqual(testCase,h([1,end]),[0,p.D_tank]);
end

function testInventoryAccounting(testCase)
[remaining,transferred]=vap_update_inventory(10,.2,.25,30);
verifyEqual(testCase,[remaining,transferred],[8.5,1.5],'AbsTol',1e-12);
[remaining,transferred]=vap_update_inventory(1,.2,1,30);
verifyEqual(testCase,[remaining,transferred],[0,1]);
[remaining,transferred]=vap_update_inventory(1,0,NaN,30);
verifyEqual(testCase,[remaining,transferred],[1,0]);
end

function testInvalidInputIsUnknownFlow(testCase)
p=vap_default_params();
r=vaporizer_thermosiphon(p.P,p.T_L,-1,p.Qdot_ref,p);
verifyFalse(testCase,r.ok);
verifyEqual(testCase,r.status,'invalid_input');
verifyTrue(testCase,isnan(r.mdot));
end

function testReferenceBalance(testCase)
p=vap_default_params();
r=vaporizer_thermosiphon(p.P,p.T_L,p.H_ref,p.Qdot_ref,p);
verifyTrue(testCase,r.ok,r.message);
verifyEqual(testCase,60*r.mdot,1.60928222448,'AbsTol',1e-6);
verifyEqual(testCase,r.T_out,61.6333127716,'AbsTol',1e-5);
verifyLessThan(testCase,abs(r.pressure_residual),1e-3);
verifyLessThan(testCase,abs(r.energy_residual),1e-6);
end

function testNoFlowHasNoOutletState(testCase)
p=vap_default_params();
r=vaporizer_thermosiphon(p.P,p.T_L,0,p.Qdot_ref,p);
verifyTrue(testCase,r.ok,r.message);
verifyEqual(testCase,r.status,'empty_tank');
verifyEqual(testCase,r.mdot,0);
verifyTrue(testCase,isnan(r.T_out));
r=vaporizer_thermosiphon(p.P,p.T_L,p.H_ref,0,p);
verifyTrue(testCase,r.ok,r.message);
verifyEqual(testCase,r.status,'zero_heat');
verifyEqual(testCase,r.mdot,0);
end

function testTemperatureCapFailureIsUnknownFlow(testCase)
p=vap_default_params();
r=vaporizer_thermosiphon(p.P,p.T_L,.01,36e3,p);
verifyFalse(testCase,r.ok);
verifyEqual(testCase,r.status,'no_root_in_domain');
verifyTrue(testCase,isnan(r.mdot));
end

function testPhaseBoundaryContinuity(testCase)
p=vap_validate_params(vap_default_params());
th=vap_thermo(p.P,p.T_L,p.fluid);
m=.02;
for h=[th.h_f,th.h_g]
    left=vap_flow_state(m,p.H_ref,m*(h-th.h_L-.01),p,th);
    right=vap_flow_state(m,p.H_ref,m*(h-th.h_L+.01),p,th);
    verifyEqual(testCase,left.rho_out,right.rho_out,'RelTol',1e-5);
    verifyEqual(testCase,left.mu_out,right.mu_out,'RelTol',1e-5);
    verifyEqual(testCase,left.T_out,right.T_out,'AbsTol',1e-3);
end
end

function testDefaultCalibration(testCase)
p=vap_default_params();
[candidate,report]=calibrate_vaporizer_loss();
verifyEqual(testCase,report.target_kg_min,.30*3.6,'AbsTol',1e-12);
verifyTrue(testCase,report.solution.ok,report.solution.message);
verifyEqual(testCase,60*report.solution.mdot,1.08,'AbsTol',1e-6);
verifyEqual(testCase,report.solution.T_out,94.1517558349,'AbsTol',1e-5);
verifyEqual(testCase,candidate.K_g,9.2633744238,'AbsTol',1e-5);
verifyEqual(testCase,candidate.K_L,p.K_L);
verifyEqual(testCase,vap_default_params(),p);
end
