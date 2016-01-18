classdef ConstantMCPAnalyticalSolution
    % Calculates the Analytic Solutions assuming constant concentrations in
    % the MCP
    properties
        pdu_params;     % params used to solve the model
        
        %concentrations based on assumption of constant a_MCP and p_MCP
        a_full_uM;      %full analytical solution
        p_full_uM;      
        
        a_satsat_uM;        %both enzymes saturated
        a_unsatunsat_uM;    %both enzymes unsaturated
        a_unsatsat_uM;      %PduCDE unsaturated
        a_satunsat_uM;      %PduPQ unsaturated
        p_satsat_uM;
        p_unsatunsat_uM;
        p_unsatsat_uM;
        p_satunsat_uM;
        
        %==================================================================
        
        
    end
    
    methods
        function obj = ConstantMCPAnalyticalSolution(pdu_params)
            obj.pdu_params = pdu_params;
            
            p=pdu_params;
            
            obj.p_full_uM=(-p.E+sqrt(p.E^2+4*p.Z))/2*p.KCDE;
            obj.a_full_uM=(-(p.F-p.U*p.W*obj.p_full_uM/(1+obj.p_full_uM))+...
                sqrt((p.F-p.U*p.W*obj.p_full_uM/(1+obj.p_full_uM))^2+...
                4*(p.U*p.W*obj.p_full_uM/(1+obj.p_full_uM)+p.V)))/2*p.KPQ;
            
            obj.p_satsat_uM=(p.Z-p.Y)*p.KCDE;
            obj.p_unsatunsat_uM=p.Z/(p.Y+1)*p.KCDE;
            obj.p_unsatsat_uM=p.Z/(p.Y+1)*p.KCDE;
            obj.p_satunsat_uM=(p.Z-p.Y)*p.KCDE;
            
            obj.a_satsat_uM=p.V+p.U*(p.W-1)*p.KPQ;        
            obj.a_unsatunsat_uM=(p.V+p.U*p.W*obj.p_unsatunsat_uM)/(p.U+1)*p.KPQ;    
            obj.a_unsatsat_uM=(p.V+p.U*(p.W*obj.p_unsatsat_uM-1))*p.KPQ;      
        	obj.a_satunsat_uM=(p.V+p.U*p.W)/(p.U+1)*p.KPQ;      
            
           
        end
        
    end
    
end

