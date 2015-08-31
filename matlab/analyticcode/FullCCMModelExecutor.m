classdef FullCCMModelExecutor < CCMModelExector    
    methods
        function obj = FullCCMModelExecutor(ccm_params)
            obj@CCMModelExector(ccm_params); 
        end
        
        function result = RunAnalytical(obj)
            p = obj.ccm_params;
 %           [HcytoRb, CcytoRb, H, C] = CHconc_Csome(p.jc, p);
            result = FullCCMAnalyticalSolution(p);
        end
    end
    
end

