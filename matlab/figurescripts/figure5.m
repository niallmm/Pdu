%script to generate Figure 5
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=100;

parameters=PduParams_MCP;
subplot(2,2,1)
sweep_paramsX2('kcA',-5,5,'Pout',-5,2,1,0,M,parameters);

parameters=PduParams_MCP;
subplot(2,2,2)
sweep_paramsX2('kmP',-6,3,'Pout',-5,2,1,0,M,parameters);

parameters=PduParams_MCP;
subplot(2,2,3)
sweep_paramsX2('kmA',-6,3,'Pout',-5,2,1,0,M,parameters);