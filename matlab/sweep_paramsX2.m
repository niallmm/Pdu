% run a sweep in some parameter for the analytical solution

% add paths
add_paths
changeplot

% define a path for saving your results
%saveLocationRoot = '/Users/niallmangan/Dropbox/CCMtesting/';
saveLocationRoot = 'C:\Users\groupadmin\Dropbox\Berkeley\Lab\pdumodeling\Pdu\matlab\testing\';
%saveLocationRoot = '/Users/chrisjakobson/Dropbox/Berkeley/Lab/pdumodeling/Pdu/matlab/testing/';

% define baseline parameters
p = PduParams_MCP;
%p.kcA = 1e-4;
%p.kcP = p.kcA;
%p.alpha =0;

% =========================================================================
% define parameter you want to change
% =========================================================================

% Define parameter sweeps in cell array
numberofsims= 50;


sweep = {'kmA',logspace(-8,8, numberofsims)
        'kmP', logspace(-8,8, numberofsims)};
% first entry is the name of the parameter as defined in the class
% (CCMParams)

a = figure;

A=zeros(length(sweep{1,2}),length(sweep{2,2}));
P=zeros(length(sweep{1,2}),length(sweep{2,2}));
Acyto=zeros(length(sweep{1,2}),length(sweep{2,2}));
Pcyto=zeros(length(sweep{1,2}),length(sweep{2,2}));

for ii = 1:length(sweep{1,2})
        startValue1=get(PduParams_MCP,sweep{1,1});
        set(p, sweep{1,1},sweep{1,2}(ii)*startValue1);
    for jj = 1:length(sweep{2,2})
        startValue2=get(PduParams_MCP,sweep{2,1});
        set(p, sweep{2,1},sweep{2,2}(jj)*startValue2);
        % set the compartment permeabilities equal !!! take out if you want
        % them to be different 
        %p.kcP = p.kcA; %
   

     % save location for this particular parameter combination run
        
        savefolder1 = savefolder(); % creates a folder name with date and time to use as save location
        saveLocation = ([saveLocationRoot savefolder1 '/']);
        %mkdir(saveLocation)
        
        %save([saveLocation 'p.mat'], 'p');
        
    % run the simulation
    res = ConstantMCPAnalyticalSolution(p);

    % save results
    %save([saveLocation 'res.mat'], 'res');
    
    A(ii,jj)=res.a_full_uM/1000;
    P(ii,jj)=res.p_full_uM/1000;
    Acyto(ii,jj)=res.a_cyto_uM/1000;
    Pcyto(ii,jj)=res.p_cyto_uM/1000;

    end
    
end

%plot results

subplot(2,2,1)
surf(log10(sweep{1,2}),log10(sweep{2,2}),log10(A));
ylabel(['parameter: ' sweep{1,1}])
xlabel(['parameter: ' sweep{2,1}])
zlabel('A concentration in compartment [mM]')

subplot(2,2,2)
surf(log10(sweep{1,2}),log10(sweep{2,2}),log10(P));
ylabel(['parameter: ' sweep{1,1}])
xlabel(['parameter: ' sweep{2,1}])
zlabel('P concentration in compartment [mM]')

subplot(2,2,3)
surf(log10(sweep{1,2}),log10(sweep{2,2}),log10(Acyto));
ylabel(['parameter: ' sweep{1,1}])
xlabel(['parameter: ' sweep{2,1}])
zlabel('A concentration in cytosol [mM]')

subplot(2,2,4)
surf(log10(sweep{1,2}),log10(sweep{2,2}),log10(Pcyto));
ylabel(['parameter: ' sweep{1,1}])
xlabel(['parameter: ' sweep{2,1}])
zlabel('P concentration in cytosol [mM]')


%line([sweep{1,2}(1) sweep{1,2}(end)],[p.KCDE/1000 p.KCDE/1000], 'Color', 'b') %saturation halfmax conc of 1,2-PD for PduCDE
%line([sweep{1,2}(1) sweep{1,2}(end)],[p.KPQ/1000 p.KPQ/1000], 'Color','r') %saturation halfmax conc of propanal for PduPQ
    
    