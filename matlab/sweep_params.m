% run a sweep in some parameter for the numerical solution

% add paths
add_paths

% define a path for saving your results
saveLocationRoot = '/Users/niallmangan/Dropbox/CCMtesting/';

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
sweep = {'jc',logspace(-2,4, numberofsims)};
% first entry is the name of the parameter as defined in the class
% (CCMParams)

a = figure;
for ii = 1:length(sweep{1,2})
    
    set(p, sweep{1,1},sweep{1,2}(ii));

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
    % plot results
    loglog(sweep{1,2}(ii), res.p_MCP_mM, 'or')
    hold on
    plot(sweep{1,2}(ii), res.a_MCP_mM, 'ob')
end

xlabel(['parameter: ' sweep{1,1}])
ylabel('A and P concentration in compartment [mM]')
    
    