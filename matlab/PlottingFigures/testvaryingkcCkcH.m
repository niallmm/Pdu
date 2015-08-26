% vary kcC and kcH separately and see what happens

p = CCMParams_Csome;
    p.kRub = 11.6; % rxns/s maximum reaction rate at single active site
    p.Km_8 = 340;    % half max reaction rate of RuBisCO, uM
	p.S_sat = 43;  % specificity ratio
    p.KO = 972; 
% kcHvec = [1e-3 1e-4 1e-5];
kcCvec = logspace(1, -5, 10);
p.pH = 8;
p.jc = 0.001;
% for j = 1:length(kcHvec)
j =1;
for i = 1:length(kcCvec)
    p.kcC = kcCvec(i);
    p.kcH = 10*p.kcC;
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    Keq_vec(i,j) = p.Keq;
    Ccsome(i,j) = res.c_csome_mM;
    Hcsome(i,j) = res.h_csome_mM;
end
% end
% 
% figure(11)
% semilogy(pHvec, Keq_vec, 'ok')
% xlabel('pH')
% ylabel('Keq')

figure(12)
loglog(kcCvec, Ccsome, 'r')
hold on
loglog(kcCvec, Hcsome, 'b')
xlabel('kcC')
ylabel('Ccsome')