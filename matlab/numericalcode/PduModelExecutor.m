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
        function results = RunNumerical(obj)
            p = obj.pdu_params;  % shorthand
            initv = zeros(2, p.xnum); % initialize vectors for concentrations
            [r, h, c, fintime, t] = driverssnondim(p.xnum, p, initv);
            results = NumericalPduModelSolution(p, r, h, c, fintime, t);
        end
    end
    
    
    
end

