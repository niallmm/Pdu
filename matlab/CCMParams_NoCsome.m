classdef CCMParams_NoCsome < CCMParams
    % Object defining CCM parameters - encapsulates various dependent
    % calculations of rates and volumes. 
    
    properties (Dependent)
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
 
        Vmax    % uM/s RuBisCO max reaction rate/concentration (inside csome)
        Vba     % maximum rate of bicarbonate dehydration by CA (inside csome)
        Vca     % maximum rate of carbon dioxide hydration by CA (inside csome)
        GC  % grouped params = D/(Rc^2 kc) + 1/Rc - 1/Rb [1/cm]
        GH
    end
    
    methods
        function obj = CCMParams_NoCsome()
            obj@CCMParams(); 
        end
        
        function jc = CalcOptimalJc(obj, Hmax)
            p = obj;
            HcytopCell = @(jc) calcHcytoDiff_NoCsome(jc, p, Hmax);
            jc = fzero(HcytopCell, 1e-2);
        end
        
        function value = get.Vmax(obj)
            value = obj.VmaxCell;
        end
        function value = get.Vba(obj)
            value = obj.VbaCell;
        end
        function value = get.Vca(obj)
            value = obj.VcaCell;
        end
        function value = get.GC
            value = 'NaN';
            warning('trying to call GC in the No-Csome case')
        end
        function value = get.GH
            value = 'NaN';
            warning('trying to call GH in the No-Csome case')
        end 
        
        function value = get.xi(obj)
            value = obj.D * obj.Kba / (obj.VbaCell * obj.Rb^2);
        end
        function value = get.gamma(obj)
            value = obj.VcaCell / obj.VbaCell;
        end
        function value = get.kappa(obj)
            value = obj.Kba / obj.Kca;
        end
        function value = get.beta_c(obj)
            value = -(obj.alpha+ obj.kmC)*obj.Rb/obj.D;
        end
        function value = get.beta_c2(obj)
            value = -obj.alpha*obj.Rb/obj.D; 
        end
        function value = get.epsilon_c(obj)
            value = -obj.kmC*obj.Cout*obj.Rb/(obj.Kca*obj.D);
        end
        function value = get.epsilon_h(obj)
            value = -(obj.jc*obj.Hout + obj.kmH*obj.Hout)*obj.Rb/(obj.D*obj.Kba);
        end
        function value = get.beta_h(obj)
            value = -obj.kmH*obj.Rb/obj.D;
        end
    end
end
