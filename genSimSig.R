
# generuje sygnał z impulsami w zadanych punktach
genSimSigAt = function(at, sigLen, peak_wd, fsample, sfq)
{
    simSig = sapply(
        1:sigLen, 
        function(k) peak_wd - min(min(abs(at - k)), peak_wd)
    )
    fltLP = butter(2, W = sfq * 3 / (fsample / 2), type = "low")
    dsim = filtfilt(fltLP, simSig) 
    dsim
    # simSig
}

# generuje pozycje impulsów
genSimPoints = function(sigLen, sbpm, fsample)
{
    bp1_sim = 3
    bp2_sim = 30
    sfq = sbpm / 60
    sp = fsample / sfq 
    stopifnot(sigLen > sp * 3)
    
    # intRng = min(sigLen, sp * 30)
    intRng = sigLen
    
    toIns = trunc(intRng / sp) + 1
    #if(sp / 2 + (toIns - 1) * sp > intRng)
    if(toIns * sp > intRng)
        toIns = toIns - 1
    #simSig = rep(0, intRng) 
    #insPts = sp / 2 + 0:(toIns - 1) * sp
    insPts = (intRng - (toIns - 1) * sp) / 2 + 0:(toIns - 1) * sp
    #simSig[insPts] = 1
 
    insPts   
}


# generator sygnału do korelowania z funkcją detekcyjną,
# zwraca sygnał o odpowiednio dobranej długości, trzeba tyle obciąć korelowany sygnał
genSimSig = function(sigLen, sbpm, fsample)
{
    insPts = genSimPoints(sigLen, sbpm, fsample)
    sfq = sbpm / 60
    peak_wd = fsample / 200 * 60 * 1/4
    dsim = genSimSigAt(insPts, sigLen, peak_wd, fsample, sfq)
    
    if(F)
    {
        plot(simSig, type = 'l')        
        plot(dsim, type = 'l')
    }
    
    dsim
}