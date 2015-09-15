classdef NumericalCCMModelSolution
    % Concentrations are in (TxN) matrices where T is the number of
    % timepoints and N is the number of points along the radius of the cell
    % considered in the discretization of the cell.
    properties
        ccm_params;     % params used to solve the model
        
        % Concentrations. Meaning depends on the model run. If you are
        % running a model with a carboxysome, then these are the
        % carboxysomal concentrations. If it is a whole cell model (i.e. no
        % carboxysome) then these are cytoplasmic. 
        h_nondim;       % nondimensional concentration of total bicarbonate over time and space.
        c_nondim;       % nondimensional concentration of co2 over time and space.
        h_mM;           % mM concentration of total bicarbonate over time and space.
        c_mM;           % mM concentration of co2 over time and space.
        
        fintime;        % final time of the numerical solution -- needs to be long enough to get to steady state
        t;              % vector of time values the numerical solver solved at -- this is only meaningful to check that we reached steady state
        r;              % radial points for concentration values.
    end
    
    methods
        function obj = NumericalCCMModelSolution(ccm_params, r, h_nondim, c_nondim, fintime, t)
            obj.ccm_params = ccm_params;
            obj.h_nondim = h_nondim;
            obj.c_nondim = c_nondim;
            obj.r = r;
            obj.t = t;
            obj.fintime = fintime;
            obj.h_mM = obj.DimensionalizeHTomM(h_nondim);
            obj.c_mM = obj.DimensionalizeCTomM(c_nondim);
        end
        
        % Converts C to mM from non-dimensional units
        function val = DimensionalizeCTomM(obj, c_nondim)
            p = obj.ccm_params;  % shorthand
            val = c_nondim * p.Kca * 1e-3;
        end
        
        % Converts H to mM from non-dimensional units
        function val = DimensionalizeHTomM(obj, h_nondim)
            p = obj.ccm_params;  % shorthand
            val = h_nondim * p.Kba * 1e-3;
        end
    end
    
end

