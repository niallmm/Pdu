%plot gradients across cell for various org cases

changeplot
p1=PduParams_NoMCP;
p2=PduParams_MCP;
p3=PduParams_MCP;
p4=PduParams_MCP;

p2.kcA=1e10;
p2.kcP=p2.kcA;
p4.kcA=1e-10;
p4.kcP=p4.kcA;

subplot(2,3,1)
[rCDE(1), rPQ(1), fluxA(1), AcytoMean(1), AMCPMean(1)] = GradientsAcrossCellNoMCP(p1,1);

subplot(2,3,2)
[rCDE(2), rPQ(2), fluxA(2), AcytoMean(2), AMCPMean(2)] = GradientsAcrossCell(p2,1);

subplot(2,3,4)
[rCDE(3), rPQ(3), fluxA(3), AcytoMean(3), AMCPMean(3)] = GradientsAcrossCell(p3,1);

subplot(2,3,5)
[rCDE(4), rPQ(4), fluxA(4), AcytoMean(4), AMCPMean(4)] = GradientsAcrossCell(p4,1);
