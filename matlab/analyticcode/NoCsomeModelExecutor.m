classdef NoCsomeModelExecutor < CCMModelExector    
    methods
        function obj = NoCsomeModelExecutor(ccm_params)
            obj@CCMModelExector(ccm_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
            result = NoCsomeAnalyticalSolution(obj.ccm_params);
        end
    end
    
end

