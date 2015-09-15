% Calculate the CO2 and HCO3- concentrations in carboxysome and at cell
% membrane in the case that there is a carboxysome-based CCM.
function [HcytoRb, CcytoRb, H, C] = CHconc_Csome(jc_in, ccm_params)

% TODO: use the CCMParams to calculate the non-dimensional params rather
% than unpacking everything here and using old code...
p = ccm_params;
alpha = p.alpha;
Rb = p.Rb;
Rc = p.Rc;
Vmax = p.Vmax;
Vca = p.Vca;
Vba = p.Vba;
Kba = p.Kba;
Kca = p.Kca;
D = p.D;
kmC = p.kmC;
kmH = p.kmH;
Km = p.Km;

H = zeros(size(jc_in));
C = zeros(size(jc_in));
HcytoRb = zeros(size(jc_in));
CcytoRb = zeros(size(jc_in));
csat = zeros(size(jc_in));
for i = 1:length(jc_in)
    jc = jc_in(i);
% =========================================================================
% calculate CO2 and HCO3- concentrations in carboxysome and at cell
% membrane
% =========================================================================
G = p.G;

N = (jc + kmH)*p.Hout*((kmC+alpha)*G + D/Rb^2) ...
    + kmC*p.Cout*((kmH+alpha)*G +D/Rb^2);
M = (kmC + alpha)*(1+Vca*Kba/(Vba*Kca))*kmH*G + ...
    kmC*(1+kmH*Vca*Kba/(kmC*Vba*Kca))*D/Rb^2;
P = ((alpha + kmC)*G + D/Rb^2).*(kmH*G + D/Rb^2);


Ccsomep = 0.5*(N./M - Rc^3*Vmax*P./(3*M*D) - Km) ...
    + 0.5*sqrt((-N./M + Rc^3*Vmax*P./(3*M*D) + Km).^2 + 4*N*Km./M);

Hcsome = Ccsomep*Vca*Kba/(Vba*Kca);

% saturated CA forward reaction

CCAsat0 = Vba*(Rc^3)*(G+D/((alpha+kmC)*Rb^2))/(3*D) + ...
            Vba*(Rc^2)/(6*D) + kmC*p.Cout/(alpha+kmC);

HCAsat0 = -Vba*(Rc^2)/D -Vba*(Rc^3)*(G+D/(kmH*Rb^2))/(3*D)...
    +(jc+kmH)*p.Hout/kmH + alpha*kmC*p.Cout./(kmH*((alpha + kmC)*G+D/Rb^2)) ...
    +(alpha-(alpha*(alpha+kmC)*G./((alpha+kmC)*G+D/Rb^2))).*CCAsat0/kmH;

if  Ccsomep>=CCAsat0
    C(i) = CCAsat0;
    H(i) = HCAsat0;
    %csat(i) = 1;
elseif Ccsomep<CCAsat0
    C(i) = Ccsomep;
    H(i) = Hcsome;
    %CAunsat =1
end


CcytoRb(i) = (kmC*p.Cout - (alpha+kmC)*C(i))*G/...
    ((alpha+kmC)*G + D/Rb^2) + C(i);

HcytoRb(i) = ((jc+kmH)*p.Hout + alpha*CcytoRb(i) - kmH*H(i))*G/(kmH*G + D/Rb^2)+H(i);
end

