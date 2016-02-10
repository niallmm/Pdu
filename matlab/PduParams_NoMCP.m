classdef PduParams_NoMCP < PduParams &  matlab.mixin.SetGet
    % Object defining CCM parameters - encapsulates various dependent
    % calculations of rates and volumes. 
    
    properties (Dependent)
        % Non-dimensional params
        % See supplementary material and pdf document NonDimEqns2 for how these
        % Formulas contain the cytosol solutions.
        xi      % ratio of rate of diffusion across cell to rate of 
                %dehydration reaction of carbonic anhydrase (D/Rc^2)/(VCDE/KPQ)
        gamma   % ratio of PduP/PduQ and PduCDE max rates (2*VPQ)/(VCDE)
        kappa   % ratio of 1/2 max PduCDE and 1/2 max PduP/PduQ concentrations (KCDE/KPQ)
        beta_a  % da/d\rho = beta_a*a + epsilon_a
        beta_p  % dp/d\rho = beta_p*p + epsilon_p
        epsilon_a % da/d\rho = beta_a*a + epsilon_a
        epsilon_p % dp/d\rho = beta_p*p + epsilon_p
        Xa  % grouped params = D/(Rc^2 kcA) + 1/Rc - 1/Rb [1/cm]
        Xp  % grouped params = D/(Rc^2 kcP) + 1/Rc - 1/Rb [1/cm]
        
        %Non-dim params for analytical solution assuming constant conc in
        %MCP
        U;
        V;
        W;
        Y;
        Z;
        E;
        F;
        Aprime;
        Bprime;
        Cprime;
        Dprime;
        Eprime;
        
        % Calculated appropriate to the volume in which the enzymes are
        % contained which depends on the situation (in cbsome or not).
        VCDE    % uM/s PduCDE max reaction rate/concentration
        VPQ     % maximum rate of aldehyde consumption by PduP/PduQ
    end
    
    methods
        function obj = PduParams_NoMCP()
            obj@PduParams(); 
        end
        
        function jc = CalcOptimalJc(obj, Hmax)
            p = obj;
            Hcytop = @(jc) calcHcytoDiff_Csome(jc, p, Hmax);
            jc = fzero(Hcytop, 1e-2); 
        end
        
        function value = get.VCDE(obj)
            value = obj.VCDECell;
        end
        function value = get.VPQ(obj)
            value = obj.VPQCell;
        end
        
        function value = get.xi(obj)
            value = obj.D * obj.KPQ / (obj.VCDE * obj.Rb^2);
        end
        function value = get.gamma(obj)
            value = 2*obj.VPQ / obj.VCDE;
        end
        function value = get.kappa(obj)
            value = obj.KCDE / obj.KPQ;
        end
        
        function value = get.Xa(obj)
            value = (obj.D/(obj.kcA*obj.Rc^2) + 1/obj.Rc - 1/obj.Rb);
        end
        function value = get.Xp(obj)
            value = (obj.D/(obj.kcP*obj.Rc^2) + 1/obj.Rc - 1/obj.Rb);
        end

        function value = get.beta_a(obj)
            value = -(obj.Rb*obj.kmA)/obj.D;
            
            %simplified form when jc=0 and kcA=>0
            %value = -obj.Rc*obj.kcA/obj.D;
        end
        function value = get.epsilon_a(obj)
            value = obj.Rb*obj.Aout*obj.kmA/(obj.D*obj.KPQ);
            
            %simplified form when jc=0 and kcA=>0
            %value = (obj.Aout*obj.Rc*obj.kcA)/(obj.D*obj.KPQ);
        end
        
        function value = get.beta_p(obj)
            value = -(obj.Rb*obj.kmP)/obj.D;
            
            %simplified form when jc=0 and kcA=>0
            %value = -obj.Rc*obj.kcP/obj.D;
        end
        
        function value = get.epsilon_p(obj)
            value = obj.Rb*obj.Pout*(obj.jc+obj.kmP)/(obj.D*obj.KCDE);
            
            %simplified form when jc=0 and kcA=>0
            %value = (obj.Pout*obj.Rc*obj.kcP)/(obj.D*obj.KCDE);
        end
        
        function value = get.U(obj)
            value =(2*obj.VPQ*obj.Rc^3*(obj.D/(obj.kmA*obj.Rb^2)+obj.Xa))/(3*obj.D*obj.KPQ);
        end
        
        function value = get.V(obj)
            value =obj.Aout/obj.KPQ;
        end
        
        function value = get.W(obj)
            value =1/2*obj.VCDE/obj.VPQ;
        end
        
        function value = get.Y(obj)
            value =(obj.VCDE*obj.Rc^3*(obj.D/(obj.kmP*obj.Rb^2)+obj.Xp))/(3*obj.D*obj.KCDE);
        end
        
        function value = get.Z(obj)
            value =(obj.Pout*(1+obj.jc/obj.kmP))/obj.KCDE;
        end
        
        function value = get.E(obj)
            value =1+obj.Y-obj.Z;
        end
        
        function value = get.F(obj)
            value =1+obj.U-obj.V;
        end
           
        function value = get.Aprime(obj)
            value = (obj.VCDE*obj.Rc)/(3*obj.KCDE*obj.kcP);
        end    
           
        function value = get.Bprime(obj)
            value = obj.Pout/obj.KCDE;
        end 
           
        function value = get.Cprime(obj)
            value = (2*obj.VPQ*obj.Rc)/(3*obj.KPQ*obj.kcA);
        end    
        
        function value = get.Dprime(obj)
            value = obj.Aout/obj.KPQ;
        end 
        
        function value = get.Eprime(obj)
            value = obj.VCDE/(2*obj.VPQ);
        end 
        
        
    end
end
