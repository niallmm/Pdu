%script to generate Figure 8
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=100;

parameters=PduParams_MCP;
sweep_paramsX2('kcA',-5,5,'kcP',-5,5,0,1,M,parameters);