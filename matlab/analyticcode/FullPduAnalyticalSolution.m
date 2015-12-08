classdef FullPduAnalyticalSolution
    % Calculates the Analytic Solutions for the whole MCP
    properties
        pdu_params;     % params used to solve the model
        a_cyto_uM;      % uM concentration of total bicarbonate in cytoplasm.
        p_cyto_uM;      % uM concentration of CO2 in cytoplasm.
        a_cyto_mM;      % mM concentration of total bicarbonate in cytoplasm.
        p_cyto_mM;      % mM concentration of CO2 in cytoplasm.
        a_MCP_uM;     % uM concentration of total bicarbonate in carboxysome.
        p_MCP_uM;     % uM concentration of CO2 in carboxysome.
        a_MCP_mM;     % mM concentration of total bicarbonate in carboxysome.
        p_MCP_mM;     % mM concentration of CO2 in carboxysome.
        
        %concentrations based on assumption of constant a_MCP and p_MCP
        a_MCP_const_uM;
        p_MCP_const_uM;
        a_MCP_const_mM;
        p_MCP_const_mM;
        
        r;
        a_cyto_rad_uM;
        p_cyto_rad_uM;
        %==================================================================
        % intermediate values useful for calculating 
        B3;
        C3;
        U;
        V;
        W;
        Y;
        Z;
        a_nondim;
        p_nondim;
        
    end
    
    methods
        function obj = FullPduAnalyticalSolution(pdu_params,aMCP,pMCP)
            obj.pdu_params = pdu_params;
            
            % Start with MCP concs from numerical solution
            p = pdu_params;
           
            %obj.a_MCP_uM = aMCP;
            %obj.p_MCP_uM = pMCP;
            
            obj.B3=(aMCP-p.Aout)/(p.D/(p.kmA*p.Rb^2)+p.Xa); %to simplify expression
           
            obj.C3=(p.kmP*pMCP-p.Pout*(p.jc+p.kmP))/(p.D/p.Rb^2+p.kmP*p.Xp);
            
            % concentration in the cytosol at r = Rb
            obj.a_cyto_uM = obj.B3*(1./p.Rb-p.D/(p.kcA*p.Rc^2)-1/p.Rc) +aMCP;
            
            obj.p_cyto_uM = obj.C3*(1./p.Rb-p.D/(p.kcP*p.Rc^2)-1/p.Rc) +pMCP;
           
           % concentration across the cell
           obj.r = linspace(p.Rc, p.Rb, 100);
           
           obj.a_cyto_rad_uM = obj.B3*(1./obj.r-p.D/(p.kcA*p.Rc^2)-1/p.Rc) +aMCP;
            
           obj.p_cyto_rad_uM = obj.C3*(1./obj.r-p.D/(p.kcP*p.Rc^2)-1/p.Rc) +pMCP;

            
           % unit conversion to mM
           obj.a_cyto_mM = obj.a_cyto_uM * 1e-3;
           obj.p_cyto_mM = obj.p_cyto_uM * 1e-3;
           obj.a_MCP_mM = aMCP* 1e-3;
           obj.p_MCP_mM = pMCP * 1e-3;
           
           %Find a_MCP and p_MCP based on assumption of constant a_MCP and p_MCP
           obj.Y=(3*p.D*p.KCDE)/(p.VCDEMCP*p.Rc^3*(p.D/(p.kmP*p.Rb^2)+p.Xp));
           obj.Z=(p.Pout*(1+p.jc/p.kmP))/p.KCDE;
           
           obj.p_nondim=(sqrt(obj.Y^2*(obj.Z+1)^2-2*obj.Y*(obj.Z-1)+1)...
               +obj.Y*(obj.Z-1)-1)/(2*obj.Y);
           
           
           obj.U=(3*p.D*p.KPQ)/(p.VCDEMCP*p.Rc^3*(p.D/(p.kmA*p.Rb^2)+p.Xp));
           obj.V=p.Aout/p.KPQ;
           obj.W=1/2*p.VCDEMCP/p.VPQMCP*obj.p_nondim/(1+obj.p_nondim);
           
           obj.a_nondim=(obj.U*(obj.V-1)+obj.W-1 ...
                +sqrt(obj.U^2*(obj.V+1)^2+2*obj.U* ...
                (obj.V*(obj.W-1)+obj.W+1)+(obj.Z-1)^2))/(2*obj.U);
           
           obj.a_MCP_const_uM=obj.a_nondim*p.KPQ;
           obj.p_MCP_const_uM=obj.p_nondim*p.KCDE;
           
           obj.a_MCP_const_mM=obj.a_MCP_const_uM/1000;
           obj.p_MCP_const_mM=obj.p_MCP_const_uM/1000;
            
           
        end
        
    end
    
end

