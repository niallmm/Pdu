%run sweeps through set of desired parameters

params={'kmP' 'kmA' 'kcA' 'Pout'};

for i=1:length(params)
    
    sweep_params(params{i});
    
end
