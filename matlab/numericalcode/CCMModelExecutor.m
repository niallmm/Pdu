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
        
        % Runs a numerical simulation of the CCM system to find the
        % steady state behavior. Stores results in this class.
        % Note: numerical code is generic to the various cases we consider
        % (e.g. no carboxysome, no ccm, etc) so this method can be
        % implemented generically.
        function results = RunNumerical(obj)
            xnum = 100; % number of points in discritization
            p = obj.pdu_params;  % shorthand
            initv = zeros(2, xnum); % initialize vectors for CO2 and HCO3- concentrations
            [r, h, c, fintime, t] = driverssnondim(xnum, p, initv);
            results = NumericalCCMModelSolution(p, r, h, c, fintime, t);
        end
    end
    
    
end

