% plot bar plot of cost of varying conditions

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

oxygenation_rate_ratio_ccm = ccm_params.O ./ (Ccsome .* ccm_params.S_sat);
oxygenation_rate_ratio_cy = ccm_params_cell.O ./ (Ccyto_cy .* ccm_params_cell.S_sat);
oxygenation_rate_ratio_just_rbc = ccm_params_just_rbc.O ./ (Ccyto_just_rbc .* ccm_params_just_rbc.S_sat);
oxygenation_rate_ratio_scaffold = ccm_params_scaffold.O ./ (Ccyto_scaffold .* ccm_params_scaffold.S_sat);

% Cost of 2PG recovery on a per-fixation basis in each case.
cost_2pg_ccm_h = cost_per_2pg * oxygenation_rate_ratio_ccm;
cost_2pg_cy_h = cost_per_2pg * oxygenation_rate_ratio_cy;
cost_2pg_just_rbc_h = cost_per_2pg * oxygenation_rate_ratio_just_rbc;
cost_2pg_scaffold = cost_per_2pg * oxygenation_rate_ratio_scaffold;


% Total cost of the system in H+/fixation for each case.
total_cost_ccm_h = cost_2pg_ccm_h + cost_h + cost_fixation;
total_cost_cy_h = cost_2pg_cy_h + cost_cy_h + cost_fixation;
total_cost_scaffold = cost_2pg_scaffold+cost_scaffold_h+cost_fixation;
total_cost_just_rbc_h = cost_2pg_just_rbc_h + cost_just_rbc_h + cost_fixation;

bar([total_cost_ccm_h total_cost_cy_h total_cost_scaffold total_cost_just_rbc_h])
% low inorganic carbon (15 uM)

% high inorganic carbon (15 mM)

%external pH for freshwater: 6.5 -7.5

% external pH for salt water: 8.1 (use to be 8.2)

