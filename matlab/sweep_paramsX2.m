function [] = sweep_paramsX2(paramToSweep1,lowerBound1,upperBound1,paramToSweep2,lowerBound2,upperBound2,equalPermSwitch,paritySwitch,N,params)
% run a sweep in two parameters for the analytical solution

% add paths
add_paths
changeplot

% define baseline parameters
p = params;

numberofsims= N;


sweep = {paramToSweep1,logspace(lowerBound1,upperBound1, numberofsims)
        paramToSweep2,logspace(lowerBound2,upperBound2, numberofsims)};

startValue1=get(params,sweep{1,1});
startValue2=get(params,sweep{2,1});

A=zeros(length(sweep{1,2}),length(sweep{2,2}));
P=zeros(length(sweep{1,2}),length(sweep{2,2}));
Acyto=zeros(length(sweep{1,2}),length(sweep{2,2}));
Pcyto=zeros(length(sweep{1,2}),length(sweep{2,2}));
fluxA=zeros(length(sweep{1,2}),length(sweep{2,2}));

for ii = 1:length(sweep{1,2})
        set(p, sweep{1,1},sweep{1,2}(ii)*startValue1);
    for jj = 1:length(sweep{2,2})
        set(p, sweep{2,1},sweep{2,2}(jj)*startValue2);
        if equalPermSwitch
            p.kcP = p.kcA; 
        end
   
    res = ConstantMCPAnalyticalSolution(p);
    
    A(ii,jj)=res.a_full_uM/1000;
    P(ii,jj)=res.p_full_uM/1000;
    Acyto(ii,jj)=res.a_cyto_uM/1000;
    Pcyto(ii,jj)=res.p_cyto_uM/1000;
    fluxA(ii,jj)=p.kmA*Acyto(ii,jj)*1000*p.SAcell/1000;
    
    end
    
end



%plot results
contourf(log10(sweep{2,2}),log10(sweep{1,2}),(P*1000/p.KCDE>1)+(A*1000/p.KPQ>1),[-0.1 0.9 1.9]);
colormap([211 211 211;248 149 33;43 172 226]./256)
hold on
contour(log10(sweep{2,2}),log10(sweep{1,2}),Acyto*1000,[0 0.01],'--r','LineWidth',1.5);
contour(log10(sweep{2,2}),log10(sweep{1,2}),Acyto*1000,[0 1],'-r','LineWidth',1.5);

line([0 0],[min(log10(sweep{1,2})) max(log10(sweep{1,2}))],'Color','k','LineStyle',':','LineWidth',1.5)
line([min(log10(sweep{2,2})) max(log10(sweep{2,2}))],[0 0],'Color','k','LineStyle',':','LineWidth',1.5)

if paritySwitch
    line([min(log10(sweep{1,2})) max(log10(sweep{1,2}))], [min(log10(sweep{2,2})) max(log10(sweep{2,2}))],'Color','g')
end

ylabel([sweep{1,1}])
xlabel([sweep{2,1}])
ax=gca;
set(ax,'FontSize',10)
ax.YTick=log10(sweep{1,2}(1:round(numberofsims/5):numberofsims));
for i=1:round(numberofsims/5):numberofsims
    labelString=num2str(startValue1.*sweep{1,2}(i),'%1.2e');
    yLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
end
ax.YTickLabel=yLabels;
ax.XTick=log10(sweep{2,2}(1:round(numberofsims/5):numberofsims));
for i=1:round(numberofsims/5):numberofsims
    labelString=num2str(startValue2.*sweep{2,2}(i),'%1.2e');
    xLabels{(i-1)/round(numberofsims/5)+1}=[labelString(1:3) 'x10^{', labelString(6:end), '}'];
end
ax.XTickLabel=xLabels;
ax.XTickLabelRotation = 45;
axis square

    
    