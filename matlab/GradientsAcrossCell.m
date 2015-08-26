% Check Gradients across cell
p = CCMParams_Csome;
p.jc = 1e-4;
p.kcC = 1e-4;
p.kcH = p.kcC;
p.pH = 8;
p.alpha =0;


exec = CCMModelExecutor(p);
res = exec.RunNumerical();

Ccsome = res.c_csome_mM;
Hcsome = res.h_csome_mM;


rb = res.rb;
Ccyto = res.c_cyto_mM;
Hcyto = res.h_cyto_mM;


figure(5)
semilogy(rb, Ccyto, 'r')
hold on
plot(rb, Hcyto, 'b')
line([0 p.Rc],[Ccsome Ccsome], 'Color','r')
line([0 p.Rc], [Hcsome Hcsome], 'Color', 'b')
xlabel('Cell radius (cm)')
ylabel('Concentration (mM)')
line([p.Rc p.Rc], [1e-6 1e6], 'Color', [0.5 0.5 0.5], 'LineStyle', '-.', 'LineWidth', 3)
xlim([0 5e-5])
ylim([1e-5 100])
set(gca,'XTick',[0 1e-5 2e-5 3e-5 4e-5 5e-5])
set(gca,'YTick',[1e-5  1e-3  1e-1  10])
box off
