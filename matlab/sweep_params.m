function [] = sweep_params(paramToSweep,lowerBound,upperBound,equalPermSwitch,kcRatio,N,params)
% run a sweep in some parameter for the numerical solution
%also allows kcA~=kcP (based on equalPermSwitch, kcRatio)

% add paths
add_paths
changeplot

% define baseline parameters
p = params;
numberofsims= N;
sweep = {paramToSweep,logspace(lowerBound,upperBound, numberofsims)};

for ii = 1:length(sweep{1,2})
    startValue=get(PduParams_MCP,sweep{1,1});
    set(p, sweep{1,1},sweep{1,2}(ii)*startValue);
    
    if equalPermSwitch
        p.kcP = kcRatio*p.kcA; %keep kcX the same
    end

    % run the simulation
    exec = FullPduModelExecutor(p);
    res = exec.RunAnalytical();
    analytical = ConstantMCPAnalyticalSolution(p);

    %concs in MCP and cytosol
    Pfull(ii)=analytical.p_full_uM/1000;
    Afull(ii)=analytical.a_full_uM/1000;
    Pcenter(ii)=res.p_MCP_uM(end,1)/1000;
    Acenter(ii)=res.a_MCP_uM(end,1)/1000;
    Pedge(ii)=res.p_MCP_uM(end,p.xnum)/1000;
    Aedge(ii)=res.a_MCP_uM(end,p.xnum)/1000;

    
end

% plot results
loglog(sweep{1,2}, Pedge, 'o','Color',[43 172 226]./256,'MarkerSize',4) %edge of MCP
hold on
plot(sweep{1,2}, Aedge, 'o','Color',[248 149 33]./256,'MarkerSize',4)
plot(sweep{1,2}, Pcenter, 'v','Color',[43 172 226]./256,'MarkerSize',4) %center of MCP
plot(sweep{1,2}, Acenter, 'v','Color',[248 149 33]./256,'MarkerSize',4)
%analytical solution assuming constant conc in MCP
plot(sweep{1,2}, Pfull, '-','Color',[43 172 226]./256)
plot(sweep{1,2}, Afull, '-','Color',[248 149 33]./256)
xlabel(['parameter: ' sweep{1,1}])
ylabel('mM')
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KCDE/1000 p.KCDE/1000], 'Color', [43 172 226]./256,'LineWidth',2) %saturation halfmax conc of 1,2-PD for PduCDE
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KPQ/1000 p.KPQ/1000], 'Color',[248 149 33]./256,'LineWidth',2) %saturation halfmax conc of propanal for PduPQ
line([1 1],ylim,'Color','k','LineStyle',':','LineWidth',1.5)
xlim([sweep{1,2}(1) sweep{1,2}(end)])
ax=gca;
set(ax,'FontSize',10)
ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
end
ax.XTickLabel=xLabels;
ax.XTickLabelRotation = 45;
axis square

    
    