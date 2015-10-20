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


sweep = {'jc',logspace(-2,4, numberofsims)
        'kcA', logspace(-6, 1, numberofsims)};
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
    exec = PduModelExecutor(p);
    res = exec.RunNumerical();

    % save results
    save([saveLocation 'res.mat'], 'res');
    

    % plot results -- currently does not differentiate different
    % permeabilities
    loglog(sweep{1,2}(ii), res.p_MCP_mM, 'or')
    hold on
    plot(sweep{1,2}(ii), res.a_MCP_mM, 'ob')
    end
    
end

xlabel(['parameter: ' sweep{1,1}])
ylabel('A and P concentration in compartment [mM]')
    
    