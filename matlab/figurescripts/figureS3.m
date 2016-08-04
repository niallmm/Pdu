%script to generate Figure S3
%CMJ 20160804
clear variables

M=50;

figure('units','normalized','outerposition',[0 0 1 1])

sweepPerformance('kmA',-5,2,1,3,M);
