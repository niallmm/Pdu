function dy = spherediffssnondim(t,y,param)
% pass parameters in PduParams struct -CMJ
% the equations in x have been discretized using finite-difference method

%numerical integration parameters
% call setgrid same function as in driver to get lattice
dx = param.dx;
x = param.x;
xnum = param.xnum;

dp = zeros(xnum, 1);        %derivative of p
da = zeros(xnum, 1);        %derivative of a

% the vector y contains both species, where y(1) = p(1), y(2) = a(1),
% y(3) = p(2), etc... we want to convert this to two vectors p = 1,2-PD, and
% a = propanal
p = y(1:2:xnum*2)';
a = y(2:2:xnum*2)';

% boundary condition at the center (symmetry)
beta = param.xi/(dx^2); % grouping of parameters
RPQ = a(1)/(1 + a(1)); %non-dimensional PduP/PduQ reaction
RCDE = p(1)/(1 + p(1)); % non-dimensional PduCDE reaction
dp(1) = beta*((3*(x(1)+x(2))/2)^(4/3))*(p(2)-p(1)) - 1/param.kappa*RCDE; % change in 1,2-PD
da(1) = beta*((3*(x(1)+x(2))/2)^(4/3))*(a(2)-a(1)) + RCDE - param.gamma*RPQ; % change in propanal

% inside MCP

for i = 2:xnum-1
    xhp = (3*(x(i)+x(i+1))/2)^(4/3); % factor from finite difference scheme in spherical coordinates
    xhm = (3*(x(i)+x(i-1))/2)^(4/3);
    RPQ = a(i)/(1 + a(i)); %non-dimensional PduP/PduQ reaction
    RCDE = p(i)/(1 + p(i)); % non-dimensional PduCDE reaction
    dp(i) = beta*(xhp*(p(i+1)-p(i))-xhm*(p(i)-p(i-1)))- 1/param.kappa*RCDE; %change in 1,2-PD
    da(i) = beta*(xhp*(a(i+1)-a(i))-xhm*(a(i)-a(i-1)))+ RCDE - param.gamma*RPQ; %change in propanal
end

%at MCP boundary

i = xnum;
    
xhp = (3*(x(i)+x(i)+0.5*dx)/2)^(2/3); % 
xhm = (3*(x(i)+x(i-1))/2)^(4/3);
dp(i) = beta*(xhp*(param.beta_p*p(i)+param.epsilon_p)*dx - xhm*(p(i)-p(i-1)));
da(i) = beta*(xhp*(param.beta_a*a(i)+param.epsilon_a)*dx - xhm*(a(i)-a(i-1)));
    

% now we need to convert the derivative vectors dp and da into one vector
dy = zeros(xnum*2,1);
dy(1:2:xnum*2)= dp;
dy(2:2:xnum*2)= da;
