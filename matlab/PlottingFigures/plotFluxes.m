% =========================================================================
% Plot CO2 and HCO3^- fluxes as a function of pH
% =========================================================================

    figure
    semilogy(pH, abs(Hin), 'r')
    hold on
    
    plot(pH, abs(Hleak), 'k')
    plot(pH, abs(Cleak), 'g')
        semilogy(pH, abs(CratewO), 'b')
    xlabel('Cytoplasmic pH')
    ylabel('inorganic carbon fluxes (picomoles/s)')
    legend('HCO_3^- transport', 'HCO_3^- leakage', 'CO_2 leakage', 'carboxylation')
    
    figure
    semilogy(pH, abs(CratewO), 'b')
    hold on
    plot(pH, abs(OratewC), 'c')
    xlabel('Cytoplasmic pH')
    ylabel('RuBisCO reactions (picomoles/s)')
    legend('carboxylation', 'oxygenation')
