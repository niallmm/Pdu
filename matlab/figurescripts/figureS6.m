%script to generate Figure S6
%CMJ 20160804
clear variables

figure

M=100;

parameters=PduParams_MCP;
subplot(1,2,1)
sweep_params_ConstantMCP('kcA',-5,5,1,0.1,M,parameters);

parameters=PduParams_MCP;
subplot(1,2,2)
sweep_params_ConstantMCP('kcA',-5,5,1,10,M,parameters);