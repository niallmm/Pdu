classdef NumericalPduModelSolution
    % Concentrations are in (TxN) matrices where T is the number of
    % timepoints and N is the number of points along the radius of the cell
    % considered in the discretization of the cell.
    properties
        pdu_params;     % params used to solve the model
        
        % Concentrations. Meaning depends on the model run. If you are
        % running a model with a MCP, then these are the
        % MCP concentrations. If it is a whole cell model (i.e. no
        % MCP) then these are cytoplasmic. 
        p_nondim;       % nondimensional concentration of 1,2-PD over time and space.
        a_nondim;       % nondimensional concentration of propanal over time and space.
        p_mM;           % mM concentration of 1,2-PD over time and space.
        a_mM;           % mM concentration of propanal over time and space.
        a_MCP_uM;     % uM concentration of propanal in MCP.
        p_MCP_uM;     % uM concentration of 1,2-PD in MCP.
        a_MCP_mM;     % mM concentration of propanal in MCP.
        p_MCP_mM;     % mM concentration of 1,2-PD in MCP.

        fintime;        % final time of the numerical solution -- needs to be long enough to get to steady state
        t;              % vector of time values the numerical solver solved at -- this is only meaningful to check that we reached steady state
        r;              % radial points for concentration values inside MCP
        
    end
    
    methods
        function obj = NumericalPduModelSolution(pdu_params, r, p_nondim, a_nondim, fintime, t)
            obj.pdu_params = pdu_params;
            obj.p_nondim = p_nondim;
            obj.a_nondim = a_nondim;
            obj.r = r;
            obj.t = t;
            obj.fintime = fintime;
            obj.p_mM = obj.DimensionalizePTomM(p_nondim);
            obj.a_mM = obj.DimensionalizeATomM(a_nondim);
            obj.a_MCP_mM=obj.a_mM(end,1); %MCP concs at r=0
            obj.p_MCP_mM=obj.p_mM(end,1);
            obj.a_MCP_uM = obj.a_mM*1e3; %MCP conc vectors as f(r)
            obj.p_MCP_uM = obj.p_mM*1e3;
        end
        
        % Converts A to mM from non-dimensional units
        function val = DimensionalizeATomM(obj, a_nondim)
            p = obj.pdu_params;  % shorthand
            val = a_nondim * p.KPQ * 1e-3;
        end
        
        % Converts P to mM from non-dimensional units
        function val = DimensionalizePTomM(obj, p_nondim)
            p = obj.pdu_params;  % shorthand
            val = p_nondim * p.KCDE * 1e-3;
        end
    end
    
end

