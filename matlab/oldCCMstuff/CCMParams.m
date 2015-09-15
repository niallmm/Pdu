classdef CCMParams
    % Object defining CCM parameters - encapsulates various dependent
    % calculations of rates and volumes. 
    
    % Mutable properties that may depend on the model context.
    properties
        jc = 0.6;        % active uptake rate of HCO3- (cm/s)
        kcC = 1e-4;        % permeability of carboxysome to CO2 (cm/s)
        kcH = 1e-4;      % permeability of carboxysome to HCO3- cm/s
        Rc = 5e-6;       % radius of c-some (cm)
        Rb = 5e-5;       % radius of cell (cm)
        D = 1e-5;        % diffusion constant (cm^2/s)

        kmC = 0.3;       % cm/s permiability of outer membrane to CO2
         alpha = 0;       % reaction rate of conversion of CO2 to HCO3- at the cell membrane (cm/s)
        Cout = 15;          % microM from Henry's law.
        % values at pH 7.8 -- will be used to re-scale pH dependence in
        % function
        kRub = 11.6;           % rxns/s maximum reaction rate at single active site
        NRub = 2000;         % number of RuBisCO active sites
        Km7_8 = 340;            % half max reaction rate of RuBisCO, uM
        KO7_8 = 972;           % uM
        O = 300;             % uM, calculated from ambient O2
        S_sat = 43;          % Specificity ratio when RuBisCO is saturated

        kCAH = 4.6e4*0.9;    % maximum rate of bicarbonate dehydration by CA /active site @ pH = 8
        kCAC = 8e4;          % maximum rate of carbon dioxide hydration by CA /active site
        NCA = 100;            % number of active sites
        Kca = 3.2*1e3;       % uM half max reaction rate for carbon dioxide hydration
        Kba = 9.3*1e3;       % uM half max reaction rate for bicarbonate dehydration
        
        pH = 8;
        kmH_base = 3e-3; % cm/s this is the permeability of the membrane to pure H2CO3
        h_cyto_exp = 3000;   %uM of inorganic carbon expected in the cytosol
        
        pH_out = 7;         % extracellular pH, 7 is around freshwater
        
        salt = 0;           % is this 1, salt water (or 0, freshwater) 

    end
    
    % Values that cannot be edited by client code since they are physical
    % constants.
    properties (Constant)
        Na = 6.022e23;       % Avogadro's number is constant, of course.
        RT = 2.4788;           % (R = 8.314e-3 kJ/(K*mol))*(298.15 K)
        
        % pKa's for Ci pool in the cytosol, calculated assuming I = 0.2
        pKa1_cyto = 3.16;   % pKa for the H2CO3 to HCO3- equilibrium
        pKa2_cyto = 9.8;    % pKa for CO3(-2) to HCO3- equilibrium
        pKa_eff_cyto = 6.1; % effective pKa for CO2 to HCO3- equilibrium
        

    end
    
    properties (Dependent)
        Vcell   % volume of cell
        Vcsome  % volume of carboxysome
        SAcell  % surface area of the cell
        
        % Dependent paramters for the case that CA & RuBisCO are uniformly 
        % co-localized to the carboxysome.
        VmaxCsome    % uM/s RuBisCO max reaction rate/concentration
        VmaxCsome_pH8 % uM/s RuBisCO max reaction rate/concentration at pH 8
        % needed to calculate the oxygen reaction rate
        VbaCsome     % maximum rate of bicarbonate dehydration by CA
        VcaCsome     % maximum rate of carbon dioxide hydration by CA
        
        % Dependent paramters for the case that CA & RuBisCO are uniformly 
        % distributed through the cytoplasm.
        VmaxCell    % uM/s RuBisCO max reaction rate/concentration
        VmaxCell_pH8 %uM/s RuBisCO max reaction rate/concentration at pH 8
        VbaCell     % maximum rate of bicarbonate dehydration by CA 
        VcaCell     % maximum rate of carbon dioxide hydration by CA 
        
        %pH dependent things
        kmH         % cm/s permiability of outer membrane to HCO3- ph dependent
        kmH_out     % same as above, but with external
        Keq         % pH dependent equilibrium constant of carbonic anhydrase
        kRub_pH     % pH dependent RuBisCO reaction rate/s at single reaction site
        Km    % pH dependent RuBisCO 1/2 max concentration
        KO          % pH dependent RuBisCO 1/2 max concentration for oxygen, 
                    % interpreted from pH dependence of Km
        
        Hout
       % Cout < ---- now in first properties section (10 uM from Henry's
       % law)
       
        % pKa's for Ci pool ouside the cell, can be for either freshwater or
        % saltwater
        pKa1_out
        pKa2_out
        pKa_eff_out
    end
    
    properties (Abstract)
        % Non-dimensional params
        % See supplementary material and pdf document NonDimEqns2 for how these
        % Formulas contain the cytosol solutions.
        xi      % ratio of rate of diffusion across cell to rate of 
                %dehydration reaction of carbonic anhydrase (D/Rb^2)/(Vba/Kba)
        gamma   % ratio of forward and backward carbonic anhydrase rates (Vba/Vca)
        kappa   % ratio of 1/2 max forward and 1/2 max backward carbonic anhydrase rates
        beta_c  % dc/d\rho = beta_c*c - epsilon_c 
        beta_c2 % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        beta_h  % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        epsilon_c % dc/dp = beta_c*c - eps_c 
        epsilon_h % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        GC  % grouped params = D/(Rc^2 kcC) + 1/Rc - 1/Rb [1/cm]
        GH  % grouped params = D/(Rc^2 kcH) + 1/Rc - 1/Rb [1/cm]
        
        % Calculated appropriate to the volume in which the enzymes are
        % contained which depends on the situation (in cbsome or not).
        Vmax    % uM/s RuBisCO max reaction rate/concentration
        Vba     % maximum rate of bicarbonate dehydration by CA 
        Vca     % maximum rate of carbon dioxide hydration by CA 
    end
    
    methods
        function value = get.Vcell(obj)
            value = 4*pi*obj.Rb^3/3;
        end
        function value = get.Vcsome(obj)
            value = 4*pi*obj.Rc^3/3;
        end
        
        function value = get.SAcell(obj)
            value = 4*pi*obj.Rb^2;
        end
        function value = get.kRub_pH(obj)
            rxn = load('pH_Vmax.mat');
            VmaxpH7_8 = interp1(rxn.pH, rxn.Vmax, 7.8, 'pchip'); % find value of Vmax at pH 8 from data
            if (obj.pH > 8.36) || (obj.pH < 6)
                warning('pH is out of range of experimental RuBisCO measurements')
                if (obj.pH < 6)
                    temp = interp1(rxn.pH, rxn.Vmax, 6); % find value of Vmax at pH we want from data
                    value = temp*obj.kRub/VmaxpH7_8;
                end
                if (obj.pH>8)
                   temp = interp1(rxn.pH, rxn.Vmax, 6);
                   value = temp*obj.kRub/VmaxpH7_8;
                end
            else
            temp = interp1(rxn.pH, rxn.Vmax, obj.pH, 'pchip'); % find value of Vmax at pH we want from data
            value = temp*obj.kRub/VmaxpH7_8; % scale the interpolated value by our known kRub rate at pH 8
            % need to scale because the data was not in the right units.
            % now it is scaled to rxns/s per active site
            end
        end
        function value = get.Km(obj)
            rxn = load('pH_Km.mat');
            KmpH7_8 = interp1(rxn.pH, rxn.Km, 7.8, 'pchip'); % find value of Vmax at pH 8 from data
            if (obj.pH > 8.36) || (obj.pH < 6)
                warning('pH is out of range of experimental RuBisCO measurements')
                
                if (obj.pH <6)
                    value = interp1(rxn.pH, rxn.Km*obj.Km7_8/KmpH7_8, 6, 'pchip'); % Km was in mM
                end
                if obj.pH>8.36
                    value = interp1(rxn.pH, rxn.Km*obj.Km7_8/KmpH7_8, 8, 'pchip');
                end
            else
                value = interp1(rxn.pH, rxn.Km*obj.Km7_8/KmpH7_8, obj.pH,'pchip'); % Km was in mM
            end
        end
        
        function value = get.KO(obj)
            rxn = load('pH_Km.mat');
            KmpH7_8 = interp1(rxn.pH, rxn.Km, 7.8, 'pchip'); % find value of Vmax at pH 8 from data
            if (obj.pH > 8.36) || (obj.pH < 6)
                warning('pH is out of range of experimental RuBisCO measurements')
                
                if (obj.pH <6)
                    value = interp1(rxn.pH, rxn.Km*obj.KO7_8/KmpH7_8, 6, 'pchip'); % Km was in mM
                end
                if obj.pH>8.36
                    value = interp1(rxn.pH, rxn.Km*obj.KO7_8/KmpH7_8, 8, 'pchip');
                end
            else
                value = interp1(rxn.pH, rxn.Km*obj.KO7_8/KmpH7_8, obj.pH, 'pchip'); % Km was in mM
            end
        end
        
        function value = get.VmaxCsome(obj)
            value = obj.kRub_pH * obj.NRub*1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        function value = get.VmaxCsome_pH8(obj)
            value = obj.kRub*obj.NRub*1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        function value = get.VbaCsome(obj)
            value = obj.kCAH * obj.NCA * 1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        function value = get.VcaCsome(obj)
            value = obj.kCAC * obj.NCA * 1e6/(obj.Vcsome * obj.Na * 1e-3);
        end
        
        function value = get.VmaxCell(obj)
            value = obj.kRub_pH * obj.NRub*1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VmaxCell_pH8(obj)
            value = obj.kRub * obj.NRub*1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VbaCell(obj)
            value = obj.kCAH * obj.NCA * 1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VcaCell(obj)
            value = obj.kCAC * obj.NCA * 1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.pKa1_out(obj)
            if obj.salt == 0
                value = 3.26;
            elseif obj.salt == 1
                value = 3.1;
            else 
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        function value = get.pKa2_out(obj)
            if obj.salt == 0
                value = 10;
            elseif obj.salt == 1
                value = 8.91;
            else 
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        function value = get.pKa_eff_out(obj)
            if obj.salt == 0
                value = 6.2;
            elseif obj.salt == 1
                value = 5.86;
            else 
                warning('salt is set to a value other than 0 or 1 in CCMParams.m')
            end
        end
        
        function value = get.kmH(obj)
            value = obj.kmH_base*10^(obj.pKa - obj.pH);
        end
        
        function value = get.Keq(obj)
            value = exp(-obj.delGHC/obj.RT - log(10)*obj.pH);
            % second term is deltaH*log(10)*pH, where deltaH = 1 from
            % HCO3- -> H20+CO2
        end

        function value = get.Hout(obj)
            delGp0CO2 =obj.delG0CO2 + obj.delG0water + log(10)*obj.RT*obj.pH_out*2;
            delGp0HCO3 = obj.delG0HCO3 + log(10)*obj.RT*obj.pH_out;
            delGp0H2CO3 = obj.delG0H2CO3 + log(10)*obj.RT*obj.pH_out*2;
            prop_H = exp(-delGp0HCO3/obj.RT)/(exp(-delGp0HCO3/obj.RT)+...
                exp(-delGp0CO2/obj.RT) + exp(-delGp0H2CO3/obj.RT));
            prop_H2 = exp(-delGp0H2CO3/obj.RT)/(exp(-delGp0HCO3/obj.RT)+...
                exp(-delGp0CO2/obj.RT) + exp(-delGp0H2CO3/obj.RT));
            B = prop_H2/(1-prop_H2);
%             value =obj.Ci_tot*prop_H;
            value = prop_H*(1+ B)*obj.Cout/(1-prop_H-B);
        end

    end
    
    methods (Abstract)
        % Calculate the optimal jc value (active bicarbonate influx).
        % Calculated differently depending on the model.
        CalcOptimalJc(obj)
    end
end
