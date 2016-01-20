classdef FullPduModelExecutor < PduModelExecutor    
    methods
        function obj = FullPduModelExecutor(pdu_params)
            obj@PduModelExecutor(pdu_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.pdu_params;
            %need numerical solution for MCP interior
            initv = zeros(2, p.xnum); % initialize vectors for CO2 and HCO3- concentrations
            [r, d, a, fintime, t] = driverssnondim(p.xnum, p, initv);
            results = NumericalPduModelSolution(p, r, d, a, fintime, t);
            aMCP=results.a_mM(end,1)*10^3; %MCP concs at r=0
            pMCP=results.p_mM(end,1)*10^3;
            result = FullPduAnalyticalSolution(p,aMCP,pMCP); %?? are these the right MCP concs to pass? don't I want r=Rc?
            result.a_MCP_uM=results.a_MCP_uM; %these are vectors in r
            result.p_MCP_uM=results.p_MCP_uM;
        end
    end
    
end

