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
        a_MCP_mM;     % mM concentration of propanal at center of MCP
        p_MCP_mM;     % mM concentration of 1,2-PD at center of MCP
        a_MCP_uM;     % uM concentration of propanal at center of MCP
        p_MCP_uM;     % uM concentration of 1,2-PD at center of MCP   
        a_cyto_uM;      % uM concentration of propanal across the cytosol
        p_cyto_uM;      % uM concentration of 1,2-PD across the cytosol 
        a_cyto_mM;      % mM concentration of propanal across the cytosol
        p_cyto_mM;      % mM concentration of 1,2-PD across the cytosol 
        
        fintime;        % final time of the numerical solution -- needs to be long enough to get to steady state
        t;              % vector of time values the numerical solver solved at -- this is only meaningful to check that we reached steady state
        r;              % radial points for concentration values inside carboxysome
        rb;             % radial points between carboxysome and cell membrane

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
            
            % concentrations at center of compartment
            obj.a_MCP_mM = obj.a_mM(end,1);
            obj.p_MCP_mM = obj.p_mM(end,1);
            obj.p_MCP_uM = obj.p_MCP_mM*1e3;
            obj.a_MCP_uM = obj.a_MCP_mM*1e3;
            
            % radial points between carboxysome and cell membrane
            obj.rb = linspace(obj.pdu_params.Rc, obj.pdu_params.Rb, 1e3);
            
            % concentrations at boundary of the cell
            
            B3=(obj.a_MCP_uM-obj.pdu_params.Aout)/(obj.pdu_params.D/(obj.pdu_params.kmA*obj.pdu_params.Rb^2)+obj.pdu_params.Xa); %to simplify expression
            C3=(obj.pdu_params.kmP*obj.p_MCP_uM-obj.pdu_params.Pout*(obj.pdu_params.jc+obj.pdu_params.kmP))/(obj.pdu_params.D/obj.pdu_params.Rb^2+obj.pdu_params.kmP*obj.pdu_params.Xp);
            
            obj.a_cyto_uM = B3*(obj.pdu_params.D/(obj.pdu_params.kcA*obj.pdu_params.Rc^2) + 1/obj.pdu_params.Rc - 1./obj.rb) + obj.a_MCP_uM;
            obj.p_cyto_uM = C3*(obj.pdu_params.D/(obj.pdu_params.kcP*obj.pdu_params.Rc^2) + 1/obj.pdu_params.Rc - 1./obj.rb) + obj.p_MCP_uM;
            obj.p_cyto_mM = obj.p_cyto_uM*1e-3;
            obj.a_cyto_mM = obj.a_cyto_uM*1e-3;
        end
        
        % Converts C to mM from non-dimensional units
        function val = DimensionalizeATomM(obj, a_nondim)
            p = obj.pdu_params;  % shorthand
            val = a_nondim * p.KPQ * 1e-3;
        end
        
        % Converts H to mM from non-dimensional units
        function val = DimensionalizePTomM(obj, p_nondim)
            p = obj.pdu_params;  % shorthand
            val = p_nondim * p.KCDE * 1e-3;
        end
    end
    
end

