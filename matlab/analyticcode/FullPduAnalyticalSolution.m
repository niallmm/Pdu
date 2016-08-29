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
        

        
        r;
        a_cyto_rad_uM;
        p_cyto_rad_uM;
        %==================================================================
        % intermediate values useful for calculating 
        B3;
        C3;

        
    end
    
    methods
        function obj = FullPduAnalyticalSolution(pdu_params,aMCP,pMCP)
            obj.pdu_params = pdu_params;
            
            % Start with MCP concs from numerical solution
            p = pdu_params;
           
           
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
           

            
           
        end
        
    end
    
end

