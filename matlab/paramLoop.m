%run sweeps through set of desired parameters

params={'kmP' 'kmA' 'kcA' 'kcP' 'Pout' 'NCDE' 'NPQ'};

for i=1:length(params)
    
    sweep_params(params{i});
    
end
