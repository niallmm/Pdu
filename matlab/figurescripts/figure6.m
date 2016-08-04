%script to generate Figure 6
%CMJ 20160804
clear variables

figure

M=100;

parameters=PduParams_MCP;
parameters.jc=1;
subplot(1,2,1)
sweep_paramsX2('kcA',-5,5,'jc',-8,5,1,0,M,parameters);

parameters=PduParams_MCP;
parameters.jc=1;
subplot(1,2,2)
sweep_paramsX2('kmP',-8,3,'jc',-8,5,1,0,M,parameters);