function [] = sweep_params_ConstantMCP(paramToSweep)
% run a sweep in some parameter for the numerical solution

% add paths
add_paths
changeplot

% define a path for saving your results
saveLocationRoot = '/Users/niallmm/Dropbox/GitHub/Pdu/matlab/CCMtesting/';
% saveLocationRoot = 'C:\Users\groupadmin\Dropbox\Berkeley\Lab\pdumodeling\Pdu\matlab\testing\';
%saveLocationRoot = '/Users/chrisjakobson/Dropbox/Berkeley/Lab/pdumodeling/Pdu/matlab/testing/';

% define baseline parameters
p = PduParams_MCP;
%p.kcA = 1e-3;
%p.kcP = p.kcA;
p.alpha =0;

% =========================================================================
% define parameter you want to change
% =========================================================================

% Define parameter sweeps in cell array
numberofsims= 50;
sweep = {paramToSweep,logspace(-10,6, numberofsims)};
% first entry is the name of the parameter as defined in the class
% (PduParams)

figure
for ii = 1:length(sweep{1,2})
    startValue=get(PduParams_MCP,sweep{1,1});
    set(p, sweep{1,1},sweep{1,2}(ii)*startValue);
    p.kcP = p.kcA; %keep kcX the same

     % save location for this particular parameter combination run
        
        savefolder1 = savefolder(); % creates a folder name with date and time to use as save location
        saveLocation = ([saveLocationRoot savefolder1 '/']);
        mkdir(saveLocation)
        
        save([saveLocation 'p.mat'], 'p');
        
    % run the simulation
    res = ConstantMCPAnalyticalSolution(p);

    % save results
    save([saveLocation 'res.mat'], 'res');
    
    % plot results
    %concs in MCP
    loglog(sweep{1,2}(ii), res.p_full_uM/1000, '-og')
    hold on
    plot(sweep{1,2}(ii), res.a_full_uM/1000, '-om')
    plot(sweep{1,2}(ii), res.p_satsat_uM/1000, 'xb')
    plot(sweep{1,2}(ii), res.a_satsat_uM/1000, 'xr')
    plot(sweep{1,2}(ii), res.p_unsatunsat_uM/1000, '+b')
    plot(sweep{1,2}(ii), res.a_unsatunsat_uM/1000, '+r')
    plot(sweep{1,2}(ii), res.a_unsatsat_uM/1000, 'sr')
    plot(sweep{1,2}(ii), res.a_satunsat_uM/1000, 'vr')
    
    legend('full P','full A','sat P','satsatA','unsat P','unsatunsatA',...
        'unsatsatA','satunsatA','Location','southeast')

    
end

xlabel(['parameter: ' sweep{1,1}])
ylabel('A and P concentration in compartment and cytosol [mM]')
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KCDE/1000 p.KCDE/1000], 'Color', 'b') %saturation halfmax conc of 1,2-PD for PduCDE
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KPQ/1000 p.KPQ/1000], 'Color','r') %saturation halfmax conc of propanal for PduPQ

    
    