function [mass_next, mass_vaporized] = vap_update_inventory(mass, mdot, x_out, dt)
%VAP_UPDATE_INVENTORY Returned liquid stays in the tank; only vapor transfers.
validateattributes(mass,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(mdot,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(dt,{'numeric'},{'scalar','real','finite','nonnegative'});
if mdot==0
    mass_vaporized=0; mass_next=mass; return
end
validateattributes(x_out,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
mass_vaporized=min(mass,x_out*mdot*dt);
mass_next=mass-mass_vaporized;
end
