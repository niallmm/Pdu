classdef NoCsomeAnalyticalSolution
    
    properties
        ccm_params;     % params used to solve the model
        h_cyto_uM;      % uM concentration of total bicarbonate.
        c_cyto_uM;      % uM concentration of total bicarbonate.
        h_cyto_mM;      % mM concentration of total bicarbonate.
        c_cyto_mM;      % mM concentration of co2.
        
        % =================================================================
        % Calculate CO2 and O2 fixation rates for whole cell
        VO;             % [uM/s] maximum rate of oxygen fixation calculated from the specificity of RuBisCO
        % intgrated over cell volume
        CratewO_um;        % [um/s] rate of CO2 fixation with oxygen accounted for
        OratewC_um;        % [um/s] rate of O2 fication with CO2 accounted for
        CratewO_pm;        % [pmole/s] rate of CO2 fixation with oxygen accounted for
        OratewC_pm;        % [pmole/s] rate of O2 fication with CO2 accounted for
        % =================================================================
        % Calculate CO2 and HCO3- flux rates at cell membrane
        % integrated over surface area of cell
        Hin_pm;            % [pmole/s] rate of active uptake of HCO3- jc*Hout
        Hleak_pm;          % [pmole/s] rate of HCO3- leakage out of cell kmH*(Hout-Hcyto)
        Cleak_pm;          % [pmole/s] rate of CO2 leakage out of cell kmC*(Cout-Ccyto)
        Hin_um;            % [um/s]
        Hleak_um;          % [um/s]
        Cleak_um;          % [um/s]
        
        error;             % the proportion of oxygen fixations to total fixation events
    end
    
    methods
        function obj = NoCsomeAnalyticalSolution(ccm_params)
            obj.ccm_params = ccm_params;
            
            % Uses analytic solution in cytosol to break down the fate of carbon in the
            % case that all the CO2 fixation maxhinery is in the cytosol.
            p = ccm_params;
            
            % CA equilibriates + Rub
            N2 = p.kmC*p.Cout+(p.jc+p.kmH)*p.Hout;
            M2 = p.kmC*(1+p.kmH*p.Vca*p.Kba/(p.kmC*p.Vba*p.Kca));
            
            CcytoRub = 0.5*(N2/M2 - p.Rb*p.Vmax/(3*M2) - p.Km)...
                +0.5*sqrt((p.Km-N2/M2 + p.Rb*p.Vmax/(3*M2)).^2 + 4*p.Km*N2/M2);
            
            HcytoRub= p.Vca*p.Kba*CcytoRub/(p.Vba*p.Kca);
            
            
            % CA saturated case
            CcytoCAsat0 = p.kmC*p.Cout/(p.alpha+p.kmC) + p.VbaCell*(p.Rb/(3*(p.alpha+p.kmC))+ p.Rb^2/(6*p.D));
            CcytoCAsatRb = p.kmC*p.Cout/(p.alpha+p.kmC) + p.VbaCell*(p.Rb/(3*(p.alpha+p.kmC)));
            HcytoCAsat0 = -p.VbaCell*(p.Rb^2/(6*p.D) + p.Rb/(3*p.kmH)) + (p.jc+p.kmH)*p.Hout/p.kmH ...
                + p.alpha*CcytoCAsatRb/p.kmH;
            
            diff= CcytoRub-CcytoCAsat0;
            
            % determine whether CA is saturated and choose apporpriate
            % analytic solution
            if  CcytoRub>CcytoCAsat0 || abs(diff)/(CcytoRub+CcytoCAsat0) <1e-3
                % if the predicted HCO3- concentration from the unsaturated case is larger 
                % than the saturated case then CA is saturated
                obj.c_cyto_uM = CcytoCAsat0;
                obj.h_cyto_uM = HcytoCAsat0;
                csat = 1;
            elseif CcytoRub<CcytoCAsat0 % if the predicted CO2 concentration
                % from the unsaturated case is smaller or equal, then CA is not
                % saturated
                obj.c_cyto_uM = CcytoRub;
                obj.h_cyto_uM = HcytoRub;
                CAunsat =1;
            end
            
            % convert to mM if we want it
            obj.h_cyto_mM = obj.h_cyto_uM * 1e-3;
            obj.c_cyto_mM = obj.c_cyto_uM * 1e-3;
            
            obj.VO = p.VmaxCell_pH8*p.KO/(p.Km_8*p.S_sat);

            
            obj.CratewO_pm = p.Vmax*obj.c_cyto_uM...
                ./(obj.c_cyto_uM+p.Km*(1+p.O/p.KO))*p.Vcell*1e3;
            obj.CratewO_um = p.Vmax*obj.c_cyto_uM...
                ./(obj.c_cyto_uM+p.Km*(1+p.O/p.KO))*p.Vcell*1e-3; % convert from uM*cm^3 to umoles
            obj.OratewC_pm = obj.VO*p.O...
                ./(p.O+p.KO*(1+obj.c_cyto_uM/p.Km))*p.Vcell*1e3;
            obj.OratewC_um = obj.VO*p.O...
                ./(p.O+p.KO*(1+obj.c_cyto_uM/p.Km))*p.Vcell*1e-3;
            
            obj.Hin_pm = p.jc*p.Hout*p.SAcell*1e3;
            obj.Hleak_pm = p.kmH*(p.Hout - obj.h_cyto_uM)*p.SAcell*1e3;
            obj.Cleak_pm = p.kmC*(p.Cout - obj.c_cyto_uM)*p.SAcell*1e3;
            
            obj.Hin_um = p.jc*p.Hout*p.SAcell*1e-3;
            obj.Hleak_um = p.kmH*(p.Hout - obj.h_cyto_uM)*p.SAcell*1e-3;
            obj.Cleak_um = p.kmC*(p.Cout - obj.c_cyto_uM)*p.SAcell*1e-3;
            
            obj.error = obj.OratewC_pm/(obj.CratewO_pm + obj.OratewC_pm);
        end
    end
    
end

