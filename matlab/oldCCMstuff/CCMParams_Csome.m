classdef CCMParams_Csome < CCMParams
    % Object defining CCM parameters - encapsulates various dependent
    % calculations of rates and volumes. 
    
    properties (Dependent)
        % Non-dimensional params
        % See supplementary material and pdf document NonDimEqns2 for how these
        % Formulas contain the cytosol solutions.
        
        xi      % ratio of rate of diffusion across cell to rate of 
                %dehydration reaction of carbonic anhydrase (D/Rc^2)/(Vba/Kba)
        gamma   % ratio of forward and backward carbonic anhydrase rates (Vba/Vca)
        kappa   % ratio of 1/2 max forward and 1/2 max backward carbonic anhydrase rates
        beta_c  % dc/d\rho = beta_c*c - epsilon_c 
        beta_c2 % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        beta_h  % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        epsilon_c % dc/dp = beta_c*c - eps_c 
        epsilon_h % dh/d\rho = beta_h*h - epsilon_h - beta_c2*c
        GC  % grouped params = D/(Rc^2 kc) + 1/Rc - 1/Rb [1/cm]
        GH
%         G
        
        Vmax    % uM/s RuBisCO max reaction rate/concentration (inside csome)
        Vba     % maximum rate of bicarbonate dehydration by CA (inside csome)
        Vca     % maximum rate of carbon dioxide hydration by CA (inside csome)
    end
    
    methods
        function obj = CCMParams_Csome()
            obj@CCMParams(); 
        end
        
        function jc = CalcOptimalJc(obj, Hmax)
            p = obj;
            Hcytop = @(jc) calcHcytoDiff_Csome(jc, p, Hmax);
            jc = fzero(Hcytop, 1e-2); 
        end
        
        function value = get.Vmax(obj)
            value = obj.VmaxCsome;
        end
        function value = get.Vba(obj)
            value = obj.VbaCsome;
        end
        function value = get.Vca(obj)
            value = obj.VcaCsome;
        end
        
        function value = get.xi(obj)
            value = obj.D * obj.Kba / (obj.VbaCsome * obj.Rc^2);
        end
        function value = get.gamma(obj)
            value = obj.VcaCsome / obj.VbaCsome;
        end
        function value = get.kappa(obj)
            value = obj.Kba / obj.Kca;
        end
        function value = get.GC(obj)
            value = (obj.D/(obj.kcC*obj.Rc^2) + 1/obj.Rc - 1/obj.Rb);
        end
        function value = get.GH(obj)
            value = (obj.D/(obj.kcH*obj.Rc^2) + 1/obj.Rc - 1/obj.Rb);
        end
%         function value = get.G(obj)
%             value = (obj.D/(obj.k*obj.Rc^2) + 1/obj.Rc - 1/obj.Rb);
%         end
        function value = get.beta_c(obj)
            value = -(obj.alpha+ obj.kmC)/(obj.Rc*((obj.kmC+obj.alpha)*obj.GC + obj.D/obj.Rb^2));
        end
        function value = get.beta_c2(obj)
            value = -obj.alpha*((obj.alpha+ obj.kmC)*obj.G/((obj.alpha+obj.kmC)*obj.GC + obj.D/obj.Rb^2)-1)*obj.Kca/(obj.Kba*obj.Rc)/(obj.kmH*obj.G +obj.D/obj.Rb^2); 
        end
        function value = get.epsilon_c(obj)
            value = -obj.kmC*obj.Cout/(obj.Kca*obj.Rc*((obj.alpha+ obj.kmC)*obj.GC+obj.D/obj.Rb^2));
        end
        function value = get.epsilon_h(obj)
            value = -(obj.jc*obj.Hout + obj.kmH*obj.Hout + obj.alpha*obj.kmC*obj.Cout*obj.GC/((obj.alpha+obj.kmC)*obj.GC+obj.D/obj.Rb^2))...
                    /(obj.Kba*obj.Rc*(obj.kmH*obj.GH+obj.D/obj.Rb^2));
        end
        function value = get.beta_h(obj)
            value = -obj.kmH/(obj.Rc*(obj.kmH*obj.GH + obj.D/obj.Rb^2));
        end
    end
end
