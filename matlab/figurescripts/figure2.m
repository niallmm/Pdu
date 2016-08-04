%script to generate Figure 2
%CMJ 20160804
clear variables

figure

M=50;

subplot(1,2,1)
doubling_time_ConstantMCP('kcA',1,M);

subplot(1,2,2)
doubling_time_ConstantMCP('Pout',1,M);