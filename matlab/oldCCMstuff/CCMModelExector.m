classdef CCMModelExector
    % class that executes the CCM model with given params
    % stores and outputs data appropriately.
    properties
        ccm_params;     % CCMParams instance for running the model.
    end
    
    methods
        function obj = CCMModelExector(ccm_params)
            obj.ccm_params = ccm_params;
        end
        
        % Runs a numerical simulation of the CCM system to find the
        % steady state behavior. Stores results in this class.
        % Note: numerical code is generic to the various cases we consider
        % (e.g. no carboxysome, no ccm, etc) so this method can be
        % implemented generically.
        function results = RunNumerical(obj, xnum)
            p = obj.ccm_params;  % shorthand
            initv = zeros(2, xnum); % initialize vectors for CO2 and HCO3- concentrations
            [r, h, c, fintime, t] = driverssnondim(xnum, p, initv);
            results = NumericalCCMModelSolution(p, r, h, c, fintime, t);
        end
    end
    
    methods (Abstract)
        % Run the analytical model
        RunAnalytical(obj)
    end
    
end

