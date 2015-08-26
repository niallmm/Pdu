%==========================================================================
% Compares the results of the analytical and numerical models in the case
% that there is no carboxysome - the enzymes are dispersed through the
% cell.
%==========================================================================

xnum = 400; % number of mesh points to discritize x
initv = zeros(2, xnum); % initialize vectors for CO2 and HCO3- concentrations

% parameters for calculation.
p_ncsome = CCMParams_NoCsome();

% Optimal jc for the no csome case.
Hmax = 30000; % 30 mM in uM
jc_opt_nocsome = p_ncsome.CalcOptimalJc(Hmax); 
p_ncsome.jc = jc_opt_nocsome;

% Numerical solution for no csome case.
executor_ncsome = NoCsomeModelExecutor(p_ncsome);
num_res_ncsome = executor_ncsome.RunNumerical(xnum);
h_ncsome_mM = num_res_ncsome.h_mM;
c_ncsome_mM = num_res_ncsome.c_mM;
r_ncsome = num_res_ncsome.r;

% Analytical solution from the no csome case
ana_res_ncsome = executor_ncsome.RunAnalytical();
Hcyto_ncsome_mM = ana_res_ncsome.h_cyto_mM;
Ccyto_ncsome_mM = ana_res_ncsome.c_cyto_mM;

% Params for the no CCM case (no CA, no active transport, no cbsome)
p_nccm = CCMParams_NoCsome();
p_nccm.NCA = 1e-6;
p_nccm.jc = 0;

% Numerical solution for no CCM case
executor_nccm = NoCsomeModelExecutor(p_nccm);
num_res_nccm = executor_nccm.RunNumerical(xnum);
h_nccm_mM = num_res_nccm.h_mM;
c_nccm_mM = num_res_nccm.c_mM;
r_nccm = num_res_nccm.r;

% Analytical solution from the no csome case
ana_res_nccm = executor_nccm.RunAnalytical();
Hcyto_nccm_mM = ana_res_nccm.h_cyto_mM;
Ccyto_nccm_mM = ana_res_nccm.c_cyto_mM;

figure(1);
semilogy(r_ncsome, h_ncsome_mM(end,:), '-r', 'LineWidth', 3);
hold on;
semilogy(r_ncsome, c_ncsome_mM(end,:), '-g', 'LineWidth', 3);
semilogy(r_ncsome, ones(size(r_ncsome)) * Hcyto_ncsome_mM, '-.c', 'LineWidth', 3);
semilogy(r_ncsome, ones(size(r_ncsome)) * Ccyto_ncsome_mM, '-.m', 'LineWidth', 3);
legend('hco3- numerical', 'co2 numerical', 'hco3- analytical', 'co2 analytical');
xlabel('radius');
ylabel('concentration uM');
title('No Carboxysome');

figure(2);
semilogy(r_nccm, h_nccm_mM(end,:), '-r', 'LineWidth', 3);
hold on;
semilogy(r_nccm, c_nccm_mM(end,:), '-g', 'LineWidth', 3);
semilogy(r_ncsome, ones(size(r_ncsome)) * Hcyto_nccm_mM, '-.c', 'LineWidth', 3);
semilogy(r_ncsome, ones(size(r_ncsome)) * Ccyto_nccm_mM, '-.m', 'LineWidth', 3);
legend('hco3- numerical', 'co2 numerical', 'hco3- analytical', 'co2 analytical');
xlabel('radius');
ylabel('concentration uM');
title('No CCM');