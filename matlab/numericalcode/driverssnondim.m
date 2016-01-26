function [r1, p1, a1, fintime, t] = driverssnondim(xnum, pdu_params, initv)
%==========================================================================
% Non-linear Diffusion Equation Solver
% uses matlab ode functions
% calls setgrid.m, initold.m, variablechangesep.m
% cytosol.m calculates solutions outside csome
% checks error
% plots solution
% for two concentrations
% 
% concentrations returned are non-dimensional and need to be converted.\
% pass parameters to ODEs (spherediffssnondim) in PduParams struct -CMJ
%==========================================================================
tic;
%==========================================================================
% Control parameters
%==========================================================================

finaltime = 100000000000;        % total simulation time
time = linspace(0,finaltime,100);
%time = [0 0.0001];
abstol = 1e-12;          % error tolerance
reltol = 1e-13;


%==========================================================================
% call setgrid.m to intialize the grid x. where u(x).
% it will return the grid used by the solver, even if that means its in a 
% different coordinate frame than solution should be ploted in
%==========================================================================

[x, dx] = setgridcsome(xnum, 1);

pdu_params.x=x; %define numerical integration settings in pdu_params
pdu_params.dx=dx;
pdu_params.xnum=xnum;

%==========================================================================
% call to set an initial condition only in the carboxysome
%==========================================================================

u0 = initold(xnum, initv);

%==========================================================================
% call ode solver
%==========================================================================

options = odeset('RelTol',reltol, 'AbsTol', abstol);
length(u0);

[t, u]=ode15s(@spherediffssnondim, time, u0, options, pdu_params); 

[fintime, junk] = size(t);

%==========================================================================
% Change variables back
%==========================================================================

[r1,p1,a1] = variablechangesep(x,u);

 toc;