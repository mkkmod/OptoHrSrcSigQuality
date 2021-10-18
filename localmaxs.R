localmaxs = function(sig, wnd2, saveMaxId = c())
{
    wnd2 = as.integer(round(wnd2))
    wnd = 0:(wnd2 * 2 - 1)
    wndlen = length(wnd)
    siglen = length(sig)
    ci = 1
    
    pre = c()
    dbg_cis = c(1)
    pcmi = NA
    
    while(TRUE)
    {
        cwnd = ci + wnd
        swnd = cwnd[cwnd < siglen]
        mi = which.max(sig[swnd])
        cmi = mi + swnd[1]
        if (mi > 1 & mi < wndlen & cmi < siglen)
            pre = c(pre, cmi)
        else if(!is.na(pcmi) & abs(cmi - pcmi) <= 1)
            pre = c(pre, cmi)
        
        pcmi = cmi
        ci = ci + wndlen
        dbg_cis = c(dbg_cis, ci)
        
        if (ci > length(sig) - wnd2)
            break
    }
    
    maxs = c()
    
    symwnd = seq(-wnd2, wnd2)
    
    for(cmi in pre)
    {
        cwnd = cmi + symwnd
        swnd = cwnd[cwnd > 1 & cwnd < siglen]
        mi = which.max(sig[swnd])
        ncmi = mi + swnd[1]
        if((abs(ncmi - cmi) <= 1) | ((1 + length(maxs)) %in% saveMaxId ) )
            maxs = c(maxs, cmi)
        # else{
        #     cat("{rm}")
        # }
    }
    
    list(maxs = maxs, pre = pre, dbg_cis = dbg_cis)
}
