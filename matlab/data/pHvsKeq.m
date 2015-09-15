% test Keq for a range of pH

p = CCMParams_Csome;

pHvec = linspace(2,11);
p.jc = 0.01;
for i = 1:length(pHvec)
    p.pH = pHvec(i);
    exec = FullCCMModelExecutor(p);
    res = exec.RunAnalytical();
    Keq_vec(i) = p.Keq;
    Ccsome(i) = res.c_csome_mM;
    Hcsome(i) = res.h_csome_mM;
end

figure(11)
semilogy(pHvec, Keq_vec, 'ok')
xlabel('pH')
ylabel('Keq')

figure(12)
plot(pHvec, Ccsome, 'r')
hold on
plot(pHvec, Hcsome, 'b')
xlabel('pH')
ylabel('Ccsome')