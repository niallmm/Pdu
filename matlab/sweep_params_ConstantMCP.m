function [] = sweep_params_ConstantMCP(paramToSweep,lowerBound,upperBound,equalPermSwitch,kcRatio,N,params)
% run a sweep in some parameter for the solution assuming constant
% concentrations in the MCP
%also allows kcA~=kcP (based on equalPermSwitch, kcRatio)

% add paths
add_paths
changeplot

% define baseline parameters
p = params;

numberofsims= N;

sweep = {paramToSweep,logspace(lowerBound,upperBound, numberofsims)};

startValue=get(params,sweep{1,1});

for ii = 1:length(sweep{1,2})
    set(p, sweep{1,1},sweep{1,2}(ii)*startValue);
    
    if equalPermSwitch
        p.kcP = kcRatio*p.kcA; %keep kcX the same
    end
    
    % run the simulation
    res = ConstantMCPAnalyticalSolution(p);

    Pfull(ii)=res.p_full_uM/1000;
    Afull(ii)=res.a_full_uM/1000;
    Psatsat(ii)=res.p_satsat_uM/1000;
    Punsatunsat(ii)=res.p_unsatunsat_uM/1000;
    Aunsatunsat(ii)=res.a_unsatunsat_uM/1000;
    Asatunsat(ii)=res.a_satunsat_uM/1000;
    Pcyto(ii)=res.p_cyto_uM/1000;
    Acyto(ii)=res.a_cyto_uM/1000;
    
end

% plot results
loglog(sweep{1,2},Pfull,'--','Color',[43 172 226]./256,'LineWidth',1.5)
hold on
plot(sweep{1,2},Pcyto, '-.','Color',[43 172 226]./256,'LineWidth',1.5)
xlim([sweep{1,2}(1) sweep{1,2}(end)])
ylim([10^-8 10^4])
plot(sweep{1,2},Afull,'--','Color',[248 149 33]./256,'LineWidth',1.5)
plot(sweep{1,2},Acyto, '-.','Color',[248 149 33]./256,'LineWidth',1.5) 
xlabel(['parameter: ' sweep{1,1}])
ylabel('mM')
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KCDE/1000 p.KCDE/1000], 'Color',[43 172 226]./256,'LineWidth',1.5) %saturation halfmax conc of 1,2-PD for PduCDE
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KPQ/1000 p.KPQ/1000], 'Color',[248 149 33]./256,'LineWidth',1.5) %saturation halfmax conc of propanal for PduPQ
line([get(params,sweep{1,1}) get(params,sweep{1,1})],[10^-8 10^4],'Color','k','LineStyle',':','LineWidth',1.5)
xlim([sweep{1,2}(1) sweep{1,2}(end)])
ylim([10^-8 10^4])
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
    
    