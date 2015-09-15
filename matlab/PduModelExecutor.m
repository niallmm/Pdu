classdef PduModelExecutor
    % class that executes the Pdu model with given params
    % stores and outputs data appropriately.
    properties
        pdu_params;     % PduParams instance for running the model.
    end
    
    methods
        function obj = PduModelExecutor(pdu_params)
            obj.pdu_params = pdu_params;
        end
        
        % Runs a numerical simulation of the Pdu system to find the
        % steady state behavior. Stores results in this class.
        % Note: numerical code is generic to the various cases we consider
        % (e.g. no MCP, etc) so this method can be
        % implemented generically.
        function results = RunNumerical(obj, xnum)
            p = obj.pdu_params;  % shorthand
            initv = zeros(2, xnum); % initialize vectors for CO2 and HCO3- concentrations
            [r, h, c, fintime, t] = driverssnondim(xnum, p, initv);
            results = NumericalPduModelSolution(p, r, h, c, fintime, t);
        end
    end
    
    methods (Abstract)
        % Run the analytical model
        RunAnalytical(obj)
    end
    
    
end

