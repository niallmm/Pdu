function [rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCellNoMCP(params,plotSwitch)
%plot the intracellular concentrations of 1,2-PD and propionaldehyde
%assuming absence of MCP

add_paths;

changeplot
% define parameters
p = params;

% run the simulation
exec = FullPduModelExecutor(p);
res = exec.RunAnalytical();

%concs in MCP and cytosol (mM)
A=res.a_MCP_uM(end,:)/1000;
P=res.p_MCP_uM(end,:)/1000;

AcytoMean=mean(A);
AMCPMean=mean(A);
rPQ=p.VPQ*mean(A)*1000/(p.KPQ+mean(A)*1000)*(p.Vcell/1000);
rCDE=p.VCDE*mean(P)*1000/(p.KCDE+mean(P)*1000)*(p.Vcell/1000);
fluxA=p.kmA*(A(end)-p.Aout/1000)*1000*p.SAcell/1000;
%if desired, also print output
%fprintf('rCDE is %4.2e umol/cell-s, rPQ is %4.2e umol/cell-s, rCDE-rPQ is %4.2e umol/cell-s and fluxA is %4.2e umol/cell-s \r',rCDE,rPQ,rCDE-rPQ,fluxA)
%fprintf('rCDE is %4.2e molecules/cell-s, rPQ is %4.2e molecules/cell-s, rCDE-rPQ is %4.2e molecules/cell-s and fluxA is %4.2e molecules/cell-s \r',rCDE*6.02*10^17,rPQ*6.02*10^17,(rCDE-rPQ)*6.02*10^17,fluxA*6.02*10^17)

% get the concentration across the cell (as a function of the radius)
rb = linspace(0, p.Rb, 100);
%plot
if plotSwitch
    semilogy(rb, A,'Color',[248 149 33]./256)
    hold on
    plot(rb, P,'Color',[43 172 226]./256)
    xlabel('Radius (cm)')
    ylabel('Concentration (mM)')
    xlim([0 5e-5])
    ylim([1e-4 100])
    plot([0 p.Rb],[p.KCDE/1000 p.KCDE/1000], '--','Color',[43 172 226]./256)
    plot([0 p.Rb],[p.KPQ/1000 p.KPQ/1000], '--','Color',[248 149 33]./256)
box off
axis square
end
