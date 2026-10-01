function vap_plot_results(S,kind,figDir)
%VAP_PLOT_RESULTS Desktop MATLAB plots. Python provides the batch equivalent.
if ~isfolder(figDir), mkdir(figDir); end
if strcmp(kind,'sweep')
    [Qg,Hg]=meshgrid(S.Q_list/1e3,S.H_list);
    for pair={{'mdot_map_H_Q',S.mdot_map*60,'Mass flow [kg/min]'}, ...
              {'T_out_map_H_Q',S.Tout_map,'Outlet temperature [K]'}}
        item=pair{1}; f=figure('Color','w');
        contourf(Qg,Hg,item{2},16,'LineColor','none'); colorbar;
        xlabel('Absorbed heater power [kW]'); ylabel('Liquid height [m]'); title(item{3});
        savefig_png(f,figDir,item{1});
    end
    data={S.Q_list/1e3,S.mdot_vs_Q*60,'mdot_vs_Q','Absorbed heater power [kW]','Mass flow [kg/min]'; ...
          S.H_list,S.mdot_vs_H*60,'mdot_vs_H','Liquid height [m]','Mass flow [kg/min]'; ...
          S.Q_list/1e3,S.x_vs_Q,'x_out_vs_Q','Absorbed heater power [kW]','Outlet vapor fraction'; ...
          S.Q_list/1e3,S.Tout_vs_Q,'T_out_vs_Q','Absorbed heater power [kW]','Outlet temperature [K]'};
    for k=1:size(data,1)
        f=figure('Color','w'); plot(data{k,1},data{k,2},'-o','LineWidth',1.5); hold on;
        if k==1
            plot(S.Q_list/1e3,S.ceiling_vs_Q*60,'--','LineWidth',1.3);
            legend('Circulation','Complete-vaporization boundary Q/(h_g-h_L)','Location','best');
        elseif k==4
            yline(S.sp.T_sat,'--','Saturation temperature');
        end
        if k==2, title(sprintf('Absorbed power = %.3g kW',S.Q_list(S.iQref)/1e3));
        else, title(sprintf('Liquid height = %.3f m',S.H_list(S.iHref))); end
        xlabel(data{k,4}); ylabel(data{k,5}); grid on;
        savefig_png(f,figDir,data{k,3});
    end
elseif strcmp(kind,'power_decay')
    f=figure('Color','w','Position',[100 100 900 1000]); tiledlayout(4,1);
    nexttile; yyaxis left; plot(S.t/60,S.Qdot/1e3); ylabel('Absorbed power [kW]');
    yyaxis right; plot(S.t/60,S.H); ylabel('Height [m]'); grid on;
    for pair={{S.mdot*60,'Mass flow [kg/min]'}, {S.T_out,'Outlet temperature [K]'}, {S.x_out,'Outlet vapor fraction'}}
        item=pair{1}; nexttile; plot(S.t/60,item{1},'LineWidth',1.5); ylabel(item{2}); xlabel('Time [min]'); grid on;
    end
    sgtitle('Quasi-steady power decay; imposed pressure and inlet temperature');
    savefig_png(f,figDir,'power_decay_transient');
else
    f=figure('Color','w','Position',[100 100 900 750]); tiledlayout(2,1);
    nexttile;
    if strcmp(kind,'SOC')
        yyaxis left; plot(S.tmin,100*S.SOC); ylabel('Volume fill [%]');
        yyaxis right; plot(S.tmin,S.H); ylabel('Height [m]');
    else
        plot(S.tmin,S.H); ylabel('Height [m]');
    end
    xlabel('Time [min]'); grid on; title('Prescribed liquid level');
    nexttile; plot(S.tmin,S.Qdot/1e3,'LineWidth',1.5); grid on;
    xlabel('Time [min]'); ylabel('Absorbed heater power [kW]'); legend(S.labels,'Location','best');
    sgtitle([kind ' scenario drivers']);
    savefig_png(f,figDir,[kind '_scenario_drivers']);
    f=figure('Color','w','Position',[100 100 900 750]); tiledlayout(2,1);
    nexttile; plot(S.tmin,60*S.mdot,'LineWidth',1.5); ylabel('Mass flow [kg/min]');
    xlabel('Time [min]'); grid on; legend(S.labels,'Location','best');
    nexttile; plot(S.tmin,S.T_out,'LineWidth',1.5); ylabel('Outlet temperature [K]');
    xlabel('Time [min]'); grid on; legend(S.labels,'Location','best');
    sgtitle(sprintf('%s scenario (failed points: %d)',kind,S.failed_count));
    savefig_png(f,figDir,[kind '_scenario_mdot_Tout']);
end
end
function savefig_png(f,dir,name)
fontsize(f,14,'points');
print(f,fullfile(dir,[name '.png']),'-dpng','-r150'); close(f);
end
