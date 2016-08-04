%calculate MCP performance ratio for a range of Pout

p1=PduParams_MCP;
p2=PduParams_NoMCP;

numberofsims= 25;

sweep = {'Pout',logspace(-8,2, numberofsims)};

% first entry is the name of the parameter as defined in the class
% (PduParams)

relPQ=zeros(numberofsims,1);
relFlux=zeros(numberofsims,1);
AcytoMCP=zeros(numberofsims,1);
AcytoNoMCP=zeros(numberofsims,1);
AMCPMCP=zeros(numberofsims,1);
AMCPNoMCP=zeros(numberofsims,1);

for ii = 1:length(sweep{1,2})
    startValue1=get(PduParams_MCP,sweep{1,1});
    set(p1, sweep{1,1},sweep{1,2}(ii)*startValue1);
    %p1.kcP = p1.kcA; %keep kcX the same
    
    startValue2=get(PduParams_NoMCP,sweep{1,1});
    set(p2, sweep{1,1},sweep{1,2}(ii)*startValue2);
    %p2.kcP = p2.kcA; %keep kcX the same
    
    [rCDE1(ii), rPQ1(ii), fluxA1(ii), Acyto1(ii), AMCP1(ii)]=GradientsAcrossCell(p1,0);
    
    [rCDE2(ii), rPQ2(ii), fluxA2(ii), Acyto2(ii), AMCP2(ii)]=GradientsAcrossCellNoMCP(p2,0);
    
    relPQ(ii)=rPQ1(ii)/rPQ2(ii);
    relFlux(ii)=fluxA2(ii)/fluxA1(ii);
    AcytoMCP(ii)=Acyto1(ii);
    AcytoNoMCP(ii)=Acyto2(ii);
    AMCPMCP(ii)=AMCP1(ii);
    AMCPNoMCP(ii)=AMCP2(ii);
    
end

figure
subplot(1,3,1)
loglog(sweep{1,2},AcytoMCP,'Color',[43 172 226]./256,'LineWidth',2)
hold on
plot(sweep{1,2},AcytoNoMCP,'Color',[248 149 33]./256,'LineWidth',2)
plot(sweep{1,2},AMCPMCP,'--','Color',[43 172 226]./256,'LineWidth',2)
legend('ACytoMCP','ACytoNoMCP','AMCPMCP','Location','SouthEast')
xlabel(['log ' sweep{1,1}])
ylabel('mM')
title('Absolute Concentrations')
axis square

subplot(1,3,2)
loglog(sweep{1,2},rPQ1,'Color',[43 172 226]./256,'LineWidth',2)
hold on
plot(sweep{1,2},rPQ2,'--','Color',[43 172 226]./256,'LineWidth',2)
plot(sweep{1,2},fluxA1,'Color',[248 149 33]./256,'LineWidth',2)
plot(sweep{1,2},fluxA2,'--','Color',[248 149 33]./256,'LineWidth',2)
legend('Flux w MCPs','Flux w/o MCPs','Leakage w MCPs','Leakage w/o MCPs','Location','East')
xlabel(['log ' sweep{1,1}])
ylabel(['\mu' 'mol/cell-s'])
title('Absolute Fluxes')
axis square

subplot(1,3,3)
loglog(sweep{1,2},relPQ,'Color',[43 172 226]./256,'LineWidth',2)
hold on
plot(sweep{1,2},relFlux,'Color',[248 149 33]./256,'LineWidth',2)
legend('Relative Carbon Flux','Relative Leakage','Location','East')
xlabel(['log ' sweep{1,1}])
ylabel('Relative Flux')
title('Relative Fluxes')
axis square


