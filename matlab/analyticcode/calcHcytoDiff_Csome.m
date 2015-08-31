% Calculate the difference between the cytoplasmic bicarbonate
% concentration and the maximum supplied. Used for finding optimal jc.
function Hdiff = calcHcytoDiff_Csome(jc, ccm_params, Hmax)
    ccm_params.h_cyto_exp = Hmax;
    ccm_params.jc = jc;
    exec = FullCCMModelExecutor(ccm_params);
    res = exec.RunAnalytical();
    Hdiff = res.hdiff;