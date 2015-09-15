% Check Gradients across cell

% add paths

% define parameters
p = PduParams_MCP;
p.jc = 5;
p.kcA = 1e-4;
p.kcP = p.kcP;
p.alpha =0;

% run the simulation
exec = FullPduModelExecutor(p);
res = exec.RunAnalytical();

% get the concentrations in the carboxysome (which is assumed to be
% constant)
AMCP = res.a_MCP_mM;
PMCP = res.p_MCP_mM;

% get the concentration across the cell (as a function of the radius)
rb = linspace(p.Rc, p.Rb, 100);
Acyto = res.a_cyto_rad_uM/10^3;
Pcyto = res.p_cyto_rad_uM/10^3;

%plot
figure(5)
semilogy(rb, Acyto, 'r')
hold on
plot(rb, Pcyto, 'b')
line([0 p.Rc],[AMCP AMCP], 'Color','r')
line([0 p.Rc], [PMCP PMCP], 'Color', 'b')
xlabel('Cell radius (cm)')
ylabel('Concentration (mM)')
%line([p.Rc p.Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 3)
%xlim([0 5e-5])
%ylim([1e-5 100])
%set(gca,'XTick',[0 1e-5 2e-5 3e-5 4e-5 5e-5])
%set(gca,'YTick',[1e-5  1e-3  1e-1  10])
box off
