% takes 3 vectors x, y, and z. finds unique values of x and y.
% converst z to a matrix who's positions correspond to the values of x and
% y and makes a contour plot of z with x and y as the dependent variables 

function [xvec, yvec, zmatrix] = xyzcontour(x,y,z)

F = scatteredInterpolant(x,y,z, 'natural', 'nearest');

[xvec, i1, junk] = unique(x);
[yvec, i1, junk] = unique(y);


    [x2,y2] = meshgrid(xvec, yvec);


    zmatrix = F(x2,y2);
