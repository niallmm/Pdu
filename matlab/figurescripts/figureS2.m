%script to generate Figure S2
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=3;

parameters=PduParams_NoMCP;
subplot(2,2,1)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCellNoMCP(parameters,1);

parameters=PduParams_MCP;
parameters.kcA=10e10;
parameters.kcP=10e10;
subplot(2,2,2)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);

parameters=PduParams_MCP;
subplot(2,2,3)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);

parameters=PduParams_MCP;
parameters.kcA=10e-10;
parameters.kcP=10e-10;
subplot(2,2,4)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);