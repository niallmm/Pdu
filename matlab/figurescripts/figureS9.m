%script to generate Figure S7
%CMJ 20160804
clear variables

M=50;

figure('units','normalized','outerposition',[0 0 1 1])

parameters=PduParams_MCP;
parameters.Aout=1;
subplot(2,3,1)
sweep_paramsX2('kcA',-5,5,'kcP',-5,5,0,1,M,parameters);

parameters=PduParams_MCP;
parameters.Aout=1;
parameters.jc=1;
subplot(2,3,4)
sweep_paramsX2('kcA',-5,5,'jc',-8,5,1,0,M,parameters);


parameters=PduParams_MCP;
parameters.Aout=1e3;
subplot(2,3,2)
sweep_paramsX2('kcA',-5,5,'kcP',-5,5,0,1,M,parameters);

parameters=PduParams_MCP;
parameters.Aout=1e3;
parameters.jc=1;
subplot(2,3,5)
sweep_paramsX2('kcA',-5,5,'jc',-8,5,1,0,M,parameters);


parameters=PduParams_MCP;
parameters.Aout=10e3;
subplot(2,3,3)
sweep_paramsX2('kcA',-5,5,'kcP',-5,5,0,1,M,parameters);

parameters=PduParams_MCP;
parameters.Aout=10e3;
parameters.jc=1;
subplot(2,3,6)
sweep_paramsX2('kcA',-5,5,'jc',-8,5,1,0,M,parameters);