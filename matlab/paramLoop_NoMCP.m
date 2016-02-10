%run sweeps through set of desired parameters

params={'kmP' 'kmA' 'kcA' 'Pout' 'NCDE' 'NPQ' 'jc'};

for i=1:length(params)
    
    sweep_params_NoMCP(params{i});
    
end
