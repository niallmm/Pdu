function [] = sweepPerformance(paramToSweep,lowerBound,upperBound,equalPermSwitch,plotSwitch,N,multiPlot,plotRow)
%calculate MCP performance ratio for a range of paramToSweep
%calculates absolute and relative PduP/Q flux and leakage values
%plotSwitch controls which plots (absolute, relative) to display

p1=PduParams_MCP;
p2=PduParams_NoMCP;

numberofsims= N;

sweep = {paramToSweep,logspace(lowerBound,upperBound, numberofsims)};

relPQ=zeros(numberofsims,1);
relFlux=zeros(numberofsims,1);
AcytoMCP=zeros(numberofsims,1);
AcytoNoMCP=zeros(numberofsims,1);
AMCPMCP=zeros(numberofsims,1);
AMCPNoMCP=zeros(numberofsims,1);

startValue1=get(PduParams_MCP,sweep{1,1});
startValue2=get(PduParams_NoMCP,sweep{1,1});

for ii = 1:length(sweep{1,2})
    
    set(p1, sweep{1,1},sweep{1,2}(ii)*startValue1);
    if equalPermSwitch
        p1.kcP = p1.kcA; %keep kcX the same
    end
    
    set(p2, sweep{1,1},sweep{1,2}(ii)*startValue2);
    if equalPermSwitch
        p2.kcP = p2.kcA; %keep kcX the same
    end
    
    [rCDE1(ii), rPQ1(ii), fluxA1(ii), Acyto1(ii), AMCP1(ii)]=GradientsAcrossCell(p1,0);
    
    [rCDE2(ii), rPQ2(ii), fluxA2(ii), Acyto2(ii), AMCP2(ii)]=GradientsAcrossCellNoMCP(p2,0);
    
    relPQ(ii)=rPQ1(ii)/rPQ2(ii);
    relFlux(ii)=fluxA2(ii)/fluxA1(ii);
    AcytoMCP(ii)=Acyto1(ii);
    AcytoNoMCP(ii)=Acyto2(ii);
    AMCPMCP(ii)=AMCP1(ii);
    AMCPNoMCP(ii)=AMCP2(ii);
    
end

if plotSwitch==1
    subplot(1,2,1)
    loglog(sweep{1,2},AcytoMCP,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},AcytoNoMCP,'Color',[102 25 255]./256,'LineWidth',2)
    plot(sweep{1,2},AMCPMCP,'--','Color',[38 230 0]./256,'LineWidth',2)
    legend('ACytoMCP','ACytoNoMCP','AMCPMCP','Location','SouthEast')
    xlabel(['parameter: ' sweep{1,1}])
    ylabel({'Absolute Concentrations' '(mM)'})
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',10)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square

    subplot(1,2,2)
    loglog(sweep{1,2},rPQ1,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},rPQ2,'Color',[102 25 255]./256,'LineWidth',2)
    plot(sweep{1,2},fluxA1,'--','Color',[38 230 0]./256,'LineWidth',2)
    plot(sweep{1,2},fluxA2,'--','Color',[102 25 255]./256,'LineWidth',2)
    legend('Flux w MCPs','Flux w/o MCPs','Leakage w MCPs','Leakage w/o MCPs','Location','SouthEast')
    xlabel(['parameter: ' sweep{1,1}])
    ylabel({'Absolute Fluxes' '\mu mol/cell-s'})
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',12)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square
end

if plotSwitch==2
    loglog(sweep{1,2},relPQ,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},relFlux,'--','Color',[38 230 0]./256,'LineWidth',2)
    legend('Relative Carbon Flux','Relative Leakage','Location','SouthEast')
    xlabel(['parameter: ' sweep{1,1}])
    ylabel('Relative Flux')
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',12)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square
end

if plotSwitch==3
    subplot(multiPlot,3,(plotRow-1)*3+1)
    loglog(sweep{1,2},AcytoMCP,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},AcytoNoMCP,'Color',[102 25 255]./256,'LineWidth',2)
    plot(sweep{1,2},AMCPMCP,'--','Color',[38 230 0]./256,'LineWidth',2)
    xlabel(['parameter: ' sweep{1,1}])
    ylabel({'Absolute Concentrations' '(mM)'})
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',12)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square

    subplot(multiPlot,3,(plotRow-1)*3+2)
    loglog(sweep{1,2},rPQ1,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},rPQ2,'Color',[102 25 255]./256,'LineWidth',2)
    plot(sweep{1,2},fluxA1,'--','Color',[38 230 0]./256,'LineWidth',2)
    plot(sweep{1,2},fluxA2,'--','Color',[102 25 255]./256,'LineWidth',2)
    xlabel(['parameter: ' sweep{1,1}])
    ylabel({'Absolute Fluxes' '\mu mol/cell-s'})
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',12)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square
    
    subplot(multiPlot,3,(plotRow-1)*3+3)
    loglog(sweep{1,2},relPQ,'Color',[38 230 0]./256,'LineWidth',2)
    hold on
    plot(sweep{1,2},relFlux,'--','Color',[38 230 0]./256,'LineWidth',2)
    xlabel(['parameter: ' sweep{1,1}])
    ylabel('Relative Flux')
    line([1 1],ylim,'LineStyle','--','Color','k')
    xlim([sweep{1,2}(1) sweep{1,2}(end)])
    ax=gca;
    set(ax,'FontSize',12)
    ax.XTick=(sweep{1,2}(1:round(numberofsims/5):numberofsims));
    for i=1:round(numberofsims/5):numberofsims
        labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
        xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
    end
    ax.XTickLabel=xLabels;
    ax.XTickLabelRotation = 45;
    axis square
    
end