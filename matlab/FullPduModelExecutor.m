classdef FullPduModelExecutor < PduModelExecutor    
    methods
        function obj = FullPduModelExecutor(pdu_params)
            obj@PduModelExecutor(pdu_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.pdu_params;
            %need numerical solution for MCP interior
            initv = zeros(2, xnum); % initialize vectors for CO2 and HCO3- concentrations
            [r, h, c, fintime, t] = driverssnondim(xnum, p, initv);
            results = NumericalPduModelSolution(p, r, h, c, fintime, t);
            aMCP=results.a_mM(end,1)*10^3;
            pMCP=results.p_mM(end,1)*10^3;
            result = FullPduAnalyticalSolution(p,aMCP,pMCP);
        end
    end
    
end

