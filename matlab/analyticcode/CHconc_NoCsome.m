% Soluions w/ reactions in whole cell
% Concentrations returned are in uM units.
function [H, C] = CHconc_NoCsome(jc_in, ccm_params)

p = ccm_params;
jc = jc_in;
alpha = p.alpha;
kmH = p.kmH;
Rb = p.Rb;
D = p.D;
kmC = p.kmC;
kmH = p.kmH;
Hout = p.Hout; 
Cout = p.Cout;
Vmax = p.Vmax;
Km = p.Km;
Kca = p.Kca;
Kba = p.Kba;
Vba = p.Vba;
Vca = p.Vca;

% jc = jcsweep;
% CA equilibriates + Rub
N2 = kmC*Cout+(jc+kmH)*Hout;
M2 = kmC*(1+kmH*Vca*Kba/(kmC*Vba*Kca));

CcytoRub = 0.5*(N2/M2 - Rb*Vmax/(3*M2) - Km)...
          +0.5*sqrt((Km-N2/M2 + Rb*Vmax/(3*M2)).^2 + 4*Km*N2/M2);

HcytoRub= Vca*Kba*CcytoRub/(Vba*Kca);
      
CcytoRubsat = N2/M2 - Rb*Vmax/3;

CcytoRubunsat = N2/(M2+Rb*Vmax/(3*Km));

% CA saturated

CcytoCAsat0 = kmC*Cout/(alpha+kmC) + Vba*(Rb/(3*(alpha+kmC))+ Rb^2/(6*D));

CcytoCAsatRb = kmC*Cout/(alpha+kmC) + Vba*(Rb/(3*(alpha+kmC)));

HcytoCAsat0 = -Vba*(Rb^2/(6*D) + Rb/(3*kmH)) + (jc+kmH)*Hout/kmH ...
                + alpha*CcytoCAsatRb/kmH;

            
if  CcytoRub>=CcytoCAsat0
    C = CcytoCAsat0;
    H = HcytoCAsat0;
    %csat = 1
elseif CcytoRub<CcytoCAsat0
    C = CcytoRub;
    H = HcytoRub;
    %CAunsat =1
end

