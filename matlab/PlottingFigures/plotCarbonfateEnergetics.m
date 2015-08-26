%Run calcCarbonfate_vs_pH.m first

% Cost of transport in the 3 cases (per fixation)
cost_h = Hin ./ CratewO;
cost_cy_h = Hin_cy ./ CratewO_cy;
cost_just_rbc_h = Hin_just_rbc ./ CratewO_just_rbc;
cost_scaffold_h = Hin_scaffold ./ CratewO_scaffold;

% Cost of fixation (per fixation)
% 3 ATP, 2 NADPH per CO2 fixed in Calvin Cycle.
% 4 H+ per ATP, 10 H+ per NADPH based on assumptions of
% ATP synthase stoichiometry and ETC H+ pumping.
cost_fixation = 4*3.0 + 10*2.0;
all_ones = ones(size(pH));

% Cost of 2PG recovery through the C2 pathway is:
% -- 1/2 carbon fixation (one decarboxylation for every two 2PG).
% -- 1 ATP (one dephosphorylation for every 2PG).
% -- 1 NADH (one oxidation for every 2PG).
% -- 1/2 amination (one net deamination for every two 2PG).
% 4 H+ per ATP, 10 H+ per NADPH based on assumptions of
% ATP synthase stoichiometry and ETC H+ pumping.
% Cost of re-aminating assumed to be 1 ATP + 1 NADPH based on the
% assumption that glutamine synthetase/glutamate synthase pathway is the
% primary pathway of ammonium incorporation as in E. coli and that
% glutamate is the product of that pathway.
cost_per_amination = 4 + 10;
cost_per_2pg = (10 + 4 + 0.5*cost_fixation + 0.5*cost_per_amination);

% Fraction of RuBisCO reactions that are oxygenations. <--- not fraction,
% but the number of oxygenations normalized by carboxylations ie
% oxygenations/carboxylations
% Specificity ratio as defined by Savir et al. 2010.
% S = (RC/RO)([O2]/[CO2]) => RO/RC = [O2]/(S [CO2])
% TODO: check that RuBisCO is saturated in all cases here.
<<<<<<< HEAD
% frac_oxygenation_ccm = ccm_params.O ./ (Ccsome .* ccm_params.S_sat);
% frac_oxygenation_cy = ccm_params_cell.O ./ (Ccyto_cy .* ccm_params_cell.S_sat);
% frac_oxygenation_just_rbc = ccm_params_just_rbc.O ./ (Ccyto_just_rbc .* ccm_params_just_rbc.S_sat);

% use the oxygenations and carboxylation rates calculated in the analytic
% and model executor files (calc_Carbonfate_vs_pH.m calls them)
frac_oxygenation_ccm = OratewC./CratewO;
frac_oxygenation_cy = OratewC_cy./CratewO_cy;
frac_oxygenation_just_rbc = OratewC_just_rbc./CratewO_just_rbc;

=======
oxygenation_rate_ratio_ccm = ccm_params.O ./ (Ccsome .* ccm_params.S_sat);
oxygenation_rate_ratio_cy = ccm_params_cell.O ./ (Ccyto_cy .* ccm_params_cell.S_sat);
oxygenation_rate_ratio_just_rbc = ccm_params_just_rbc.O ./ (Ccyto_just_rbc .* ccm_params_just_rbc.S_sat);
oxygenation_rate_ratio_scaffold = ccm_params_scaffold.O ./ (Ccyto_scaffold .* ccm_params_scaffold.S_sat);
>>>>>>> FETCH_HEAD

% Cost of 2PG recovery on a per-fixation basis in each case.
cost_2pg_ccm_h = cost_per_2pg * oxygenation_rate_ratio_ccm;
cost_2pg_cy_h = cost_per_2pg * oxygenation_rate_ratio_cy;
cost_2pg_just_rbc_h = cost_per_2pg * oxygenation_rate_ratio_just_rbc;
cost_2pg_scaffold = cost_per_2pg * oxygenation_rate_ratio_scaffold;

% Total cost of the system in H+/fixation for each case.
total_cost_ccm_h = cost_2pg_ccm_h + cost_h + cost_fixation;
total_cost_cy_h = cost_2pg_cy_h + cost_cy_h + cost_fixation;
total_cost_just_rbc_h = cost_2pg_just_rbc_h + cost_just_rbc_h + cost_fixation;
<<<<<<< HEAD
% =========================================================================
% Plot cummulative plot
% =========================================================================
%%
plot1 = 0;
if plot1 == 1
    % Cost breakdown of the ccm model as a function of pH.
    figure(5)
    subplot(1,3,1);
    area(pH, total_cost_ccm_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
    hold on;
    area(pH, cost_fixation + cost_2pg_ccm_h, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
    area(pH, cost_2pg_ccm_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
    plot(pH, total_cost_ccm_h, '-k', 'LineWidth', 3);
    set(gca, 'Yscale', 'log');
    axis([pH(1) pH(end) 1e-1 2e8]);
    ylabel('Energetic Cost H^+ / (CO_2 fixed)');
    xlabel('Cytoplasmic pH');
    title('Full CCM');
    
    % Cost breakdown of the cytosolic enzymes model as a function of pH.
    subplot(1,3,2);
    area(pH, total_cost_cy_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
    hold on;
    area(pH, cost_fixation + cost_2pg_cy_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
    area(pH, all_ones .* cost_fixation, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
    plot(pH, total_cost_cy_h, '-k', 'LineWidth', 3);
    set(gca, 'Yscale', 'log');
    axis([pH(1) pH(end) 1e-1 2e8]);
    set(gca, 'YTick', []);
    xlabel('Cytoplasmic pH');
    title('Reactions in Cytosol');
    
    % Cost breakdown of the cytosolic enzymes model as a function of pH.
    subplot(1,3,3);
    area(pH, total_cost_just_rbc_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
    hold on;
    area(pH, cost_fixation + cost_2pg_just_rbc_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
    area(pH, all_ones .* cost_fixation, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
    plot(pH, total_cost_just_rbc_h, '-k', 'LineWidth', 3);
    set(gca, 'Yscale', 'log');
    axis([pH(1) pH(end) 1e-1 2e8]);
    set(gca, 'YTick', []);
    xlabel('Cytoplasmic pH');
    title('RuBisCO Alone');
    legend('Transport', 'Photorespiration', 'Fixation', 'Total Cost');
end
%%
% =========================================================================
% Plot cost of transport, CCB cycle and photorespiration per carbon fixation
% =========================================================================
=======
total_cost_scaffold_h = cost_2pg_scaffold + cost_scaffold_h + cost_fixation;

% Fraction oxygenation as a function of pH. 
figure(222)
semilogy(pH, oxygenation_rate_ratio_ccm, 'b');
hold on;
semilogy(pH, oxygenation_rate_ratio_cy, 'g');
semilogy(pH, oxygenation_rate_ratio_just_rbc, 'r');
semilogy(pH, oxygenation_rate_ratio_scaffold, 'k');
title('Ratio of Oxygenations/Carboxylations');
legend('CCM', 'cytosolic enzymes', 'RuBisCO alone', 'scaffolded enzymes');

% Cost breakdown of the ccm model as a function of pH.
figure(5)
subplot(2,2,1);
area(pH, total_cost_ccm_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
hold on;
area(pH, cost_fixation + cost_2pg_ccm_h, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
area(pH, cost_2pg_ccm_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
plot(pH, total_cost_ccm_h, '-k', 'LineWidth', 3);
plot(pH, 3e3 .* all_ones, '--g', 'LineWidth', 2);
set(gca, 'Yscale', 'log');
axis([pH(1) pH(end) 1e-1 2e8]);
ylabel('Energetic Cost H^+ / (CO_2 fixed)');
xlabel('Cytoplasmic pH');
title('Full CCM');

% Cost breakdown of the scaffolded enzymes model as a function of pH.
subplot(2,2,2);
area(pH, total_cost_scaffold_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
hold on;
area(pH, cost_fixation + cost_2pg_scaffold, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
area(pH, all_ones .* cost_fixation, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
plot(pH, total_cost_scaffold_h, '-k', 'LineWidth', 3);
plot(pH, 3e3 .* all_ones, '--g', 'LineWidth', 2);
set(gca, 'Yscale', 'log');
axis([pH(1) pH(end) 1e-1 2e8]);
xlabel('Cytoplasmic pH');
title('Scaffolded Enzymes');
legend('Transport', 'Photorespiration', 'Fixation', 'Total Cost', 'Original Model');

% Cost breakdown of the cytosolic enzymes model as a function of pH.
subplot(2,2,3);
area(pH, total_cost_cy_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
hold on;
area(pH, cost_fixation + cost_2pg_cy_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
area(pH, all_ones .* cost_fixation, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
plot(pH, total_cost_cy_h, '-k', 'LineWidth', 3);
plot(pH, 3e3 .* all_ones, '--g', 'LineWidth', 2);
set(gca, 'Yscale', 'log');
axis([pH(1) pH(end) 1e-1 2e8]);
set(gca, 'YTick', []);
xlabel('Cytoplasmic pH');
title('Reactions in Cytosol');

% Cost breakdown of the cytosolic enzymes model as a function of pH.
subplot(2,2,4);
area(pH, total_cost_just_rbc_h, 'FaceColor', [.337,.165,.447], 'BaseValue', 1e-1);
hold on;
area(pH, cost_fixation + cost_2pg_just_rbc_h, 'FaceColor', [.667, .518, .224], 'BaseValue', 1e-1);
area(pH, all_ones .* cost_fixation, 'FaceColor', [.251, .557, .184], 'BaseValue', 1e-1);
plot(pH, total_cost_just_rbc_h, '-k', 'LineWidth', 3);
plot(pH, 3e3 .* all_ones, '--g', 'LineWidth', 2);
set(gca, 'Yscale', 'log');
axis([pH(1) pH(end) 1e-1 2e8]);
set(gca, 'YTick', []);
xlabel('Cytoplasmic pH');
title('RuBisCO Alone');

>>>>>>> FETCH_HEAD
figure(2)
% plot cost of transport/CO2 fixation with Cell membrane permeability
loglog(kmH, cost_h, 'Color', 'k', 'LineWidth', 3);
hold on
<<<<<<< HEAD
=======
text(10*kmH(end), 1.7*cost_h(end), 'CCM cost', 'Color', 'r');

loglog(kmH, cost_cy_h, 'Color', 'k', 'LineWidth', 3);
plot([1e-7 1e-7], [1e-1 2e7], 'Color', 'm', 'LineWidth', 20)
text(1.7e-7, 5e4, 'HCO_3^- permeability', 'Color', 'm')

plot([3e-3 3e-3], [1e-1 2e7], 'Color', 'c', 'LineWidth', 20)
text(7e-5, 5e4, 'H_2CO_3 permeability', 'Color', 'c')


plot([kmH(1) kmH(end)], [cost_fixation cost_fixation], '--b');
%text(kmH(1)/1100.0, 50, 'Calvin Cycle Fixation', 'Color', 'b');

plot([3e-4 3e-4], [1e-1 3e3], '--k');
plot([kmH(1) kmH(end)], [3e3 3e3], '--k');
%text(kmH(1)/1100.0, 5e3, 'Original Model', 'Color', 'k');

% cost is 1 nadph, 1 atp and 1 co2 times specificity ratio
% when RuBisCO is saturated: S = vc/vo ~ 13;
% not saturated: S = (vc/Kc)/(vo/Ko)~ 50 %SavirMilo paper
% Is rubisco saturated?
oxygenation_rate_ratio_cell = ccm_params.O ./ (Ccyto_cy .* ccm_params.S_sat);
cost_2pg_cell_h = oxygenation_rate_ratio_cell .* (10 + 4 + 1.5*4 + 10);
total_cost_cell_h = cost_2pg_cell_h + cost_cy_h + cost_fixation;

plot(kmH, total_cost_cell_h, '--g');
text(kmH(1)/1100.0, 1.5*total_cost_cell_h(end), 'No CCM', 'Color', 'g');

oxygenation_rate_ratio_noca = ccm_params.O ./ (Ccyto_just_rbc .* ccm_params.S_sat);
cost_2pg_noca_h = oxygenation_rate_ratio_noca .* (10 + 4 + 1.5*4 + 10);
total_cost_noca_h = cost_2pg_noca_h + cost_just_rbc_h + cost_fixation;

plot(kmH, total_cost_noca_h, '--m');
%text(kmH(1)/1100.0, 1.5*total_cost_noca_h(end), 'Just RuBisCO', 'Color', 'g');

>>>>>>> FETCH_HEAD
haxes1 = gca; % handle to axes
set(haxes1,'XColor','k','YColor','k')
axis([kmH(end) kmH(1) 1e-1 2e7])
xlabel('k_m^H, Effective Membrane Permeability of HCO_3^* (cm/s)')
haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes2 = axes('Position',haxes1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Yscale', 'log', ...
    'Color','none');
hold on
% plot cost of transport/CO2 fixation with pH
loglog(pH, cost_h,'Parent', haxes2,'Color','r', 'LineWidth', 3);
set(haxes2, 'XDir','Reverse')
axis([pH(1) pH(end) 1e-1 2e7])
ylabel('Energetic Cost H^+ / (CO_2 fixed)')
xlabel('Cytoplasmic pH')
text(5, 5e2, 'HCO3^- transport cost', 'Color', 'r')
% plot cost of fixation/CO2 (flat)
plot([pH(1) pH(end)], [cost_fixation cost_fixation], '--b');
text(5, 15, 'Calvin Cycle Fixation', 'Color', 'b');
% plot cost of photorespiration
loglog(pH, cost_2pg_ccm_h, 'Parent', haxes2,'Color','c', 'LineWidth', 3);
text(7.5, 1, 'Cost of photorespiration', 'Color', 'c')

%
% loglog(kmH, cost_cy_h, 'Color', 'k', 'LineWidth', 3);
% plot([1e-7 1e-7], [1e-1 2e7], 'Color', 'm', 'LineWidth', 20)
% text(1.7e-7, 5e4, 'HCO_3^- permeability', 'Color', 'm')
%
% plot([3e-3 3e-3], [1e-1 2e7], 'Color', 'c', 'LineWidth', 20)
% text(7e-5, 5e4, 'H_2CO_3 permeability', 'Color', 'c')
%
%
% plot([kmH(1) kmH(end)], [cost_fixation cost_fixation], '--b');
% text(kmH(1)/1100.0, 50, 'Calvin Cycle Fixation', 'Color', 'b');
%
% plot([3e-4 3e-4], [1e-1 3e3], '--k');
% plot([kmH(1) kmH(end)], [3e3 3e3], '--k');
% %text(kmH(1)/1100.0, 5e3, 'Original Model', 'Color', 'k');
%
% % cost is 1 nadph, 1 atp and 1 co2 times specificity ratio
% % when RuBisCO is saturated: S = vc/vo ~ 13;
% % not saturated: S = (vc/Kc)/(vo/Ko)~ 50 %SavirMilo paper
% % Is rubisco saturated?
% frac_oxygenation_cell = ccm_params.O ./ (Ccyto_cy .* ccm_params.S_sat);
% cost_2pg_cell_h = frac_oxygenation_cell .* (10 + 4 + 1.5*4 + 10);
% total_cost_cell_h = cost_2pg_cell_h + cost_cy_h + cost_fixation;
%
% plot(kmH, total_cost_cell_h, '--g');
% text(kmH(1)/1100.0, 1.5*total_cost_cell_h(end), 'No CCM', 'Color', 'g');
%
% frac_oxygenation_noca = ccm_params.O ./ (Ccyto_just_rbc .* ccm_params.S_sat);
% cost_2pg_noca_h = frac_oxygenation_noca .* (10 + 4 + 1.5*4 + 10);
% total_cost_noca_h = cost_2pg_noca_h + cost_just_rbc_h + cost_fixation;
%
% plot(kmH, total_cost_noca_h, '--m');
% text(kmH(1)/1100.0, 1.5*total_cost_noca_h(end), 'Just RuBisCO', 'Color', 'g');



% =========================================================================
% Plot CO2 and HCO3^- concentration in carboxysome as a function of pH
% =========================================================================
plot3 = 0;
if plot3 ==1
    figure(111)
    plot(pH, Ccsome, 'r')
    hold on
    plot(pH, Hcsome, 'b')
    xlabel('pH')
    ylabel('inorganic carbon')
end

