% find critical jc of 
function Hdiff = calcHcytoDiff_NoCsome(jc_in, ccm_params, Hmax)
    [H, unused_C] = CHconc_NoCsome(jc_in, ccm_params);
    Hdiff = H - Hmax;