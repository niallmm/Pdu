classdef FullPduModelExecutor < PduModelExecutor    
    methods
        function obj = FullPduModelExecutor(pdu_params)
            obj@PduModelExecutor(pdu_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.pdu_params;
            initv = zeros(2, p.xnum); % initialize vectors for concentrations
            [r, d, a, fintime, t] = driverssnondim(p.xnum, p, initv);
            results = NumericalPduModelSolution(p, r, d, a, fintime, t);
            aMCP=results.a_mM(end,end)*10^3; %MCP concs at r=Rc
            pMCP=results.p_mM(end,end)*10^3;
            result = FullPduAnalyticalSolution(p,aMCP,pMCP);
            result.a_MCP_uM=results.a_MCP_uM; %these are vectors in r
            result.p_MCP_uM=results.p_MCP_uM;
        end
    end
    
end

