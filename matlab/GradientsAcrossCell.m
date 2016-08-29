function [rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(params,plotSwitch)
%plot the intracellular concentrations of 1,2-PD and propionaldehyde
%assuming constant concentrations inside the MCP

% add paths
add_paths;
changeplot

% define parameters
p = params;

% run the simulation
res = ConstantMCPAnalyticalSolution(p);

%concs in MCP and cytosol (mM)
A=res.a_full_uM/1000;
P=res.p_full_uM/1000;
Acyto=res.a_cyto_rad_uM/1000;
Pcyto=res.p_cyto_rad_uM/1000;

AcytoMean=mean(Acyto);
AMCPMean=mean(A);
rPQ=2*p.VPQ*A*1000/(p.KPQ+A*1000)*(p.VMCP/1000);
rCDE=p.VCDE*P*1000/(p.KCDE+P*1000)*(p.VMCP/1000);
fluxA=p.kmA*(Acyto(end)-p.Aout/1000)*1000*p.SAcell/1000;
%fprintf('rCDE is %4.2e umol/cell-s, rPQ is %4.2e umol/cell-s, rCDE-rPQ is %4.2e umol/cell-s and fluxA is %4.2e umol/cell-s \r',rCDE,rPQ,rCDE-rPQ,fluxA)
%fprintf('rCDE is %4.2e molecules/cell-s, rPQ is %4.2e molecules/cell-s, rCDE-rPQ is %4.2e molecules/cell-s and fluxA is %4.2e molecules/cell-s \r',rCDE*6.02*10^17,rPQ*6.02*10^17,(rCDE-rPQ)*6.02*10^17,fluxA*6.02*10^17)

% get the concentration across the cell (as a function of the radius)
rb = linspace(p.Rc, p.Rb, 100);

%plot
if plotSwitch
    semilogy(rb, Acyto,'Color',[248 149 33]./256)
    hold on
    plot(rb, Pcyto,'Color',[43 172 226]./256)
    plot([0 p.Rc],[A A],'Color',[248 149 33]./256)
    plot([0 p.Rc], [P P],'Color',[43 172 226]./256)
    xlabel('Radius (cm)')
    ylabel('Concentration (mM)')
    line([p.Rc p.Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 1)
    xlim([0 5e-5])
    ylim([1e-4 100])
    plot([0 p.Rb],[p.KCDE/1000 p.KCDE/1000], '--','Color',[43 172 226]./256)
    plot([0 p.Rb],[p.KPQ/1000 p.KPQ/1000], '--','Color',[248 149 33]./256)
    box off

axis square
end
