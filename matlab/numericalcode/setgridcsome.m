%==========================================================================
% Make the discretized grid to solve equations on.
% For solving the inner sphere only.
% Grid spacing is for flux, concentration sits in between grid spacings.
% This is why there are xnum+1 grid spacings and xnum concentration pts.
%==========================================================================


function [x, dx1] = setgridcsome(xnum,Rc)
r1 = Rc;

%switch to m coordinates
m1 = (r1^3)/3;                      

% grid spacing 
dx1 = m1/(xnum);                 

% make grid
x = linspace(0, m1, xnum);

