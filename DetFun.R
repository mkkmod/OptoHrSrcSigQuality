library(signal)

filter_hr = function(sig, fsample, bp1, bp2)
{
    half_fsample = fsample / 2
    fltBP = butter(2, W = c(bp1 / half_fsample, bp2 / half_fsample),  type="pass")
    signal::filter(fltBP, sig, init.x = rep(sig[1], 4))
}

detFun = function(sig, fsample, bp1, bp2, lp)
{
    half_fsample = fsample / 2
    fltLP = butter(4, W = lp / half_fsample, type="low")
   
    bp_sig = filter_hr(sig, fsample, bp1, bp2)
    signal::filter(fltLP, bp_sig^2)
}

detFunDiff = function(sig, fsample, bp1, bp2, lp)
{
    half_fsample = fsample / 2
    fltLP = butter(4, W = lp / half_fsample, type="low")
    
    bp_sig = filter_hr(sig, fsample, bp1, bp2)
    #signal::filter(fltLP, abs(diff(bp_sig))) # słabsze
    signal::filter(fltLP, diff(bp_sig)^2)
}
