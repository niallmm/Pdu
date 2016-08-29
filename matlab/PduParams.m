classdef PduParams< matlab.mixin.SetGet
    % Object defining Pdu parameters - encapsulates various dependent
    % calculations of rates and volumes. 
    
    % Mutable properties that may depend on the model context.
    % includes parameters for numerical integration -CMJ
    properties
        jc = 0;        % active uptake rate of 1,2-PD(cm/s)
        kcA = 1e-5;      % permeability of MCP to propanal (cm/s)
        kcP = 1e-5;      % permeability of MCP to 1,2-PD (cm/s)
        Rb = 5e-5;       % radius of cell (cm)
        Rc = 1e-5;       % radius of MCP (cm) %made larger to match real volume ratio
        D = 1e-5;        % diffusion constant (cm^2/s)

        kmA = 0.01;       % cm/s permeability of outer membrane to propanal
        kmP = 0.01;       % cm/s permeability of outer membrane to 1,2-PD
        
        alpha = 0;          % reaction rate of conversion of CO2 to HCO3- at the cell membrane (cm/s)
        Pout = 55*1e3;      % uM concentration of 1,2-PD outside
        Aout = 0;           % uM concentration of propanal outside

        kcatCDE = 300;          % rxns/s maximum reaction rate at single PduCDE active site
        NCDE = 1500;              % number of PduCDE active sites %updated based on MFS data
        KCDE= 0.5*1e3;        	% half max reaction rate of PduCDE, uM
        kcatPQ = 55;            % rxns/s maximum rate of aldehyde consumption by PduP/PduQ
        NPQ = 2500;              % number of PduP/PduQ active sites %updated based on MFS data
        KPQ = 15*1e3;           % uM half max reaction rate for PduP/PduQ
        
        x;
        dx;
        xnum=100;
        
    end
    
    % Values that cannot be edited by client code since they are physical
    % constants.
    properties (Constant)
        Na = 6.022e23;       % Avogadro's number is constant, of course.
        RT = 2.4788;           % (R = 8.314e-3 kJ/(K*mol))*(298.15 K)

    end
    
    properties (Dependent)
        Vcell   % volume of cell
        VMCP  % volume of MCP
        SAcell  % surface area of the cell
        
        % Dependent paramters for the case that PduCDE & PduP/Q are uniformly 
        % co-localized to the MCP.
        VCDEMCP    % uM/s PduCDE max reaction rate/concentration
        VPQMCP     % uM/s maximum rate of aldehyde consumption by PduP/PduQ
        
        % Dependent paramters for the case that PduCDE & PduP/Q are uniformly 
        % distributed through the cytoplasm.
        VCDECell    % uM/s PduCDE max reaction rate/concentration
        VPQCell     % uM/s maximum rate of aldehyde consumption by PduP/PduQ

    end
    
    properties (Abstract)
        % Non-dimensional params
        xi      % ratio of rate of diffusion across cell to rate of 
                %dehydration reaction of carbonic anhydrase (D/Rc^2)/(VCDE/KPQ)
        gamma   % ratio of PduP/PduQ and PduCDE max rates (2*VPQ)/(VCDE)
        kappa   % ratio of 1/2 max PduCDE and 1/2 max PduP/PduQ concentrations (KCDE/KPQ)
        beta_a  % da/d\rho = beta_a*a + epsilon_a
        beta_p  % dp/d\rho = beta_p*p + epsilon_p
        epsilon_a % da/d\rho = beta_a*a + epsilon_a
        epsilon_p % dp/d\rho = beta_p*p + epsilon_p
        Xa  % grouped params = D/(Rc^2 kcA) + 1/Rc - 1/Rb [1/cm]
        Xp  % grouped params = D/(Rc^2 kcP) + 1/Rc - 1/Rb [1/cm]
        
        % Calculated appropriate to the volume in which the enzymes are
        % contained which depends on the situation (in MCP or not).
        VCDE    % uM/s PduCDE max reaction rate/concentration
        VPQ     % maximum rate of aldehyde consumption by PduP/PduQ
    end
    
    methods
        function value = get.Vcell(obj)
            value = 4*pi*obj.Rb^3/3;
        end
        function value = get.VMCP(obj)
            value = 4*pi*obj.Rc^3/3;
        end
        
        function value = get.SAcell(obj)
            value = 4*pi*obj.Rb^2;
        end
       
        
        function value = get.VCDEMCP(obj)
            value = obj.kcatCDE * obj.NCDE*1e6/(obj.VMCP * obj.Na * 1e-3);
        end
        function value = get.VPQMCP(obj)
            value = obj.kcatPQ * obj.NPQ*1e6/(obj.VMCP * obj.Na * 1e-3);
        end
   
        function value = get.VCDECell(obj)
            value = obj.kcatCDE * obj.NCDE*1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        function value = get.VPQCell(obj)
            value = obj.kcatPQ * obj.NPQ*1e6/(obj.Vcell * obj.Na * 1e-3);
        end
        


    end
    
end
