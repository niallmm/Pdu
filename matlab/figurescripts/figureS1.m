%script to generate Figure S1
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=25;

parameters=PduParams_MCP;
subplot(2,3,1)
sweep_params('kcA',-5,5,1,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,2)
sweep_params('Pout',-5,2,1,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,4)
sweep_params('kmA',-6,3,1,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,5)
sweep_params('kmP',-6,3,1,1,M,parameters);