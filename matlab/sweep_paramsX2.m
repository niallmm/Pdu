% run a sweep in some parameter for the numerical solution

% add paths
add_paths
changeplot

% define a path for saving your results
%saveLocationRoot = '/Users/niallmangan/Dropbox/CCMtesting/';
saveLocationRoot = 'C:\Users\groupadmin\Dropbox\Berkeley\Lab\pdumodeling\Pdu\matlab\testing\';

% define baseline parameters
p = PduParams_MCP;
p.kcA = 1e-4;
p.kcP = p.kcA;
p.alpha =0;

% =========================================================================
% define parameter you want to change
% =========================================================================

% Define parameter sweeps in cell array
numberofsims= 3;


sweep = {'jc',logspace(-3,3, numberofsims)
        'kcA', logspace(-3, 3, numberofsims)};
% first entry is the name of the parameter as defined in the class
% (CCMParams)

a = figure;

for ii = 1:length(sweep{1,2})
        set(p, sweep{1,1},sweep{1,2}(ii));
    for jj = 1:length(sweep{2,2})
        set(p, sweep{2,1}, sweep{2,2}(jj));
        % set the compartment permeabilities equal !!! take out if you want
        % them to be different 
        p.kcP = p.kcA; %
   

     % save location for this particular parameter combination run
        
        savefolder1 = savefolder(); % creates a folder name with date and time to use as save location
        saveLocation = ([saveLocationRoot savefolder1 '/']);
        mkdir(saveLocation)
        
        save([saveLocation 'p.mat'], 'p');
        
    % run the simulation
    exec = FullPduModelExecutor(p);
    res = exec.RunAnalytical();

    % save results
    save([saveLocation 'res.mat'], 'res');
    
    Acyto = res.a_cyto_rad_uM/10^3;
    Pcyto = res.p_cyto_rad_uM/10^3;

    % plot results -- currently does not differentiate different
    % permeabilities
    loglog(sweep{1,2}(ii), res.p_MCP_mM, 'or')
    hold on
    plot(sweep{1,2}(ii), res.a_MCP_mM, 'ob')
    plot(sweep{1,2}(ii), mean(Pcyto), 'xb')
    plot(sweep{1,2}(ii), mean(Acyto), 'xr')
    end
    
end

xlabel(['parameter: ' sweep{1,1}])
ylabel('A and P concentration in compartment [mM]')
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KCDE/1000 p.KCDE/1000], 'Color', 'b') %saturation halfmax conc of 1,2-PD for PduCDE
line([sweep{1,2}(1) sweep{1,2}(end)],[p.KPQ/1000 p.KPQ/1000], 'Color','r') %saturation halfmax conc of propanal for PduPQ
    
    