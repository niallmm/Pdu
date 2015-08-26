function [r1, h1, c1, fintime, t] = driverssnondim(xnum, ccm_params, initv)
%==========================================================================
% Non-linear Diffusion Equation Solver
% uses matlab ode functions
% calls setgrid.m, initold.m, variablechangesep.m
% cytosol.m calculates solutions outside csome
% checks error
% plots solution
% for two concentrations
% 
% concentrations returned are non-dimensional and need to be converted.
%==========================================================================
tic;
%==========================================================================
% Control parameters
%==========================================================================

p = ccm_params;
parameters = [xnum p.xi p.gamma p.kappa p.beta_h ...
    p.epsilon_h p.beta_c p.epsilon_c p.beta_c2 p.Vmax ...
    p.Km p.Vba p.Kca];

finaltime = 100000000000;        % total simulation time
time = linspace(0,finaltime,100);
%time = [0 0.0001];
abstol = 1e-12;          % error tolerance
reltol = 1e-13;

%xnum =400;             % number of mesh points
xnum = parameters(1);

%==========================================================================
% call setgrid.m to intialize the grid x. where u(x).
% it will return the grid used by the solver, even if that means its in a 
% different coordinate frame than solution should be ploted in
%==========================================================================

[x, dx] = setgridcsome(xnum, 1);

param = [parameters dx x]; 

%param = [xnum xi gamma kappa beta_h epsilon_h beta_c epsilon_c beta_c2 dx x];
% xi = param(2);
% param2 = [xnum Rc Rb D k Jc Vmax Km Vba Kba Vca Kca kmH kmC Hout Cout alpha dx x];

%==========================================================================
% call to set an initial condition only in the carboxysome
%==========================================================================

u0 = initold(xnum, initv);

%==========================================================================
% call ode solver
%==========================================================================

options = odeset('RelTol',reltol, 'AbsTol', abstol);
length(u0);

[t, u]=ode15s(@spherediffssnondim, time, u0, options, param); 

[fintime, junk] = size(t);

%==========================================================================
% Change variables back
%==========================================================================

[r1,h1,c1] = variablechangesep(x,u);

 toc;