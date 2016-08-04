%script to generate Figure S7
%CMJ 20160804
clear variables

M=3;

figure

sweepPerformance('kcA',-5,2,1,3,M);

figure

sweepPerformance('kcA',-5,2,0,3,M);

figure

sweepPerformance('kcP',-5,2,0,3,M);