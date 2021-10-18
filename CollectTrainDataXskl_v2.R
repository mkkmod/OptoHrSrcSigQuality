# skrypt do tworzenia zestawu danych do uczenia optym filtru opartego o sieci TCN 
# bezpośrednio odczytujący pliki .xskl

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())

library(bspec)
library(signal)
library(doParallel) 
library(reticulate)

source('DetFun.R')
source("sd_cor_to_rem_idx.R")
source('genSimSig.R')
source('localmaxs.R')

use_virtualenv("/home/mkrej/MyPy3Venv/", required = T)

py_run_string("import numpy as np")
py_run_string("import clr")
xsklLoaderLibPath = "./bin/Debug/XsklSimpLoaderLib"
py_run_string(sprintf('sys.path.append("%s")', xsklLoaderLibPath))
py_run_string('clr.AddReference("XsklSimpLoaderLib")')


py_run_string('from XsklSimpLoaderLib import XsklSimpLoader')
py_run_string('ld = XsklSimpLoader()')

hrAlgLibPath ="./bin/Debug/OptoOnnxHrAlg"
py_run_string(sprintf('sys.path.append("%s")', hrAlgLibPath))

py_run_string('clr.AddReference("OptoOnnxHrAlg")')
py_run_string('from OptoOnnxHrAlg import OptoHrAlgCompPlain')
py_run_string('hrAlg = OptoHrAlgCompPlain()')

recDir = path.expand('~/QnapSmb/badania2018/Received')
savDir = '../ExampData/_online_16s_v2'

# najpierw: sudo mount -t cifs -o username= ...
# p. HowToDocs/Opto.txt
if(!dir.exists(recDir))
    stop("brak katalogu z zapisami, może niezamontowany?? \n=> ", recDir) 

print(savDir)

# losowanie pliku
allSesFiles = dir(recDir, "iMonSession.xskl", recursive = T)
fname = file.path(recDir, sample(allSesFiles, 1))
dbgStartTm = NA


forcePitchBpm = NA

# debugowe zapisy
if(T)
{
    # bad500fname("25048-J-G") dodany
    fname = "/home/mkrej/QnapSmb/badania2018/Received/xskl_ses_2019_02_05__11-31-29_m564_sav_2019_02_05__11-51-47/iMonSession.xskl"
    dbgStartTm = 650
}
    
print(fname)
stopifnot(file.exists(fname))


py_run_string(sprintf('fname = "%s"', fname))
py_run_string("ld.Load(fname)")

# OptoSigNN/compTCN_xskl.py
py_run_string(
    paste(
        c(    
            "sigids = list(ld.SigIdList())",
            "peaksids = [i for i in sigids if i.startswith('peakch-')]",
            "allinvs = list(list(ld.GetSigStartAndInterval(sigid))[1] for sigid in peaksids)",
            "if not all(x==allinvs[0] for x in allinvs):",
            "    raise NotImplementedError('not all interv equal !! [', ', '.join(str(x) for x in allinvs), ']')",
            "interv = allinvs[0]",
            "start = min(list(list(ld.GetSigStartAndInterval(sigid))[0] for sigid in peaksids))",
            "siglst = []",
            "for sigid in peaksids:",
            "    ld.PrintSigInfo(sigid)",
            "    psig = ld.GetEqSampSigPadRep(sigid, start)",  # Repeat Padding !!
            "    sig = np.fromiter(psig, float)",
            "    siglst.append(sig)",
            "minlen = min(len(sig) for sig in siglst)",
            "siglst = list(list(sig[0:minlen] for sig in siglst))"
        ), collapse = "\n")
    )

siglst = py_eval("siglst")


if(F)
{
    plot(siglst[[1]][5000+1:15000], type="l")
    length(siglst)
    plot(siglst[[2]][5000+1:15000], type="l")
}


# num [1:8192] 1.53e-06 1.53e-06 1.53e-06 1.53e-06 1.53e-06 ...


fsample = 1 / py_eval("interv")

saveLen = 2^14

peak_wd = fsample * 60 / (250 * 2)

py$fsample = fsample
# czas pływającego normowania
normtm = 3.0
py$normlen = trunc(normtm * fsample)
py$normMaxOutAmp = 100

py$PitchToWindowMul = 0.4
py_run_string('hrAlg.SetPitchToWindowMul(PitchToWindowMul)')


normAndDet = function(selDet, startPitchBpm, constPitch, normOFF = FALSE)
{
    #py$selDet = c(rev(selDet[1:py$normlen]), selDet, selDet[(length(selDet)-fsample*6):length(selDet)])
    py$selDet = c(rev(selDet[1:py$normlen]), selDet)
    #plot(py$selDet, type = 'l', col='red4')
    
    if (!normOFF)
        selDetNorm = py_eval("np.fromiter(ld.MovingQuantilesNormalization(selDet, normlen, 0xF, normMaxOutAmp), float)")
    else
        selDetNorm= py$selDet
    
    #selDetNorm[1:py$normlen] = 0
    selDetNorm = c(selDetNorm, rep(0, fsample * 10))
    #selDetNorm = selDetNorm[(py$normlen + 1):length(selDetNorm)]
    py$selDetNorm = selDetNorm
    py$startPitchBpm = startPitchBpm
    py$constPitch = constPitch
    
    maxs = py_eval("np.fromiter(hrAlg.Comp(selDetNorm, fsample, startPitchBpm, constPitch), int)")
    if(F)
    {
        plot(selDetNorm, type = 'l', col='red4')
        lines(maxs, selDetNorm[maxs], type="p", col = 'blue4', lwd = 2)
        exp = data.frame(tm = 1:length(selDetNorm) / fsample, detNorm = selDetNorm)    
        write.table(exp, "../dbg_sig.csv", sep = "; ", row.names = F, quote = F)
    }
    maxs = maxs - py$normlen
    maxs = maxs[maxs > 0]
    
    selDetNorm = selDetNorm[(py$normlen + 1):length(selDetNorm)]
    selDetNorm = selDetNorm[1:length(selDet)]
    
    torem = sd_cor_to_rem_idx(maxs)
    rem_maxs = maxs[torem]
    if (length(torem))
    {
        maxs = maxs[-torem]
    }
    
    stopifnot(!any(rem_maxs %in% maxs))
    
    list(selDetNorm = selDetNorm, maxs = maxs, rem_maxs = rem_maxs)
}


bp1_sig = 3
bp2_sig = 10
lp_sig = 4

if(F)
    dbgStartTm = dbgStartTm + 10

# *****************************************
# check time range
{
    if(is.na(dbgStartTm))
    {
        # losowanie start time
        startTm = round(sample(((1+py$normlen):((length(siglst[[1]]) - 1.5 * saveLen)) / fsample), 1))
    }else
    {
        # debugowy start time
        startTm = dbgStartTm
    }

    cat(sprintf("startTm = %d\n", startTm))
    
    startId = startTm * fsample
    selRng = startId + 1:saveLen
    

    # new sel sig
    allMaxs = list()
    allDet = list()
    allHrSd = data.frame(id = 1:length(siglst))
    
    for(id in 1:length(siglst))
    {
        #selSig = allSig[, id]
        selSig = 1e9 * siglst[[id]][selRng]
        allHrSd$sd[id] = NA
        if (sd(selSig) < 1e-5)
        {
            cat(sprintf("fid: %d SKIPPING ZERO SIG\n", id))
            next
        }
        selDet = detFunDiff(selSig, fsample, bp1_sig, bp2_sig, lp_sig)
        
        if(is.na(forcePitchBpm))
            nad = normAndDet(selDet, startPitchBpm = 120, constPitch = FALSE)
        else
            nad = normAndDet(selDet, startPitchBpm = forcePitchBpm, constPitch = TRUE)
        
        selDetNorm = nad$selDetNorm
        maxs = nad$maxs
        rem_maxs = nad$rem_maxs
        
        allHrSd$sd[id] =  sd(diff(maxs))
        allMaxs[[id]] = maxs
        allDet[[id]] = selDetNorm
    }
    
    print(allHrSd[order(allHrSd$sd), ])
    best_fid = allHrSd[which.min(allHrSd$sd), ]$id
    cat(sprintf("best_fid = %g\n", best_fid))
    
    par(mfcol = c(4, 1), mar = c(2, 2, 2, 1))
    plot(NA, xlim = c(0, saveLen), ylim = c(-1, 10))
    tresh_sd = quantile(allHrSd$sd, c(.3), na.rm=T)
    
    selDiffsVec = c()
    selMedDiffs = c()
    
    for(id in 1:nrow(allHrSd))
    {
        if(is.na(allHrSd$sd[id]))
            next
        if (allHrSd$sd[id] > tresh_sd )
            next
        
        #selSig = allSig[, id]
        
        maxs = allMaxs[[id]]
        selDiffsVec = c(selDiffsVec, diff(maxs))
        selMedDiffs = c(selMedDiffs, median(diff(maxs)))
        selDetNorm = allDet[[id]]
        
        lines(selDetNorm, type = 'l', col='red4')
        lines(maxs, selDetNorm[maxs], type="p", col = 'black', lwd = 2)
        
    }
    
    pitchBpm = round(60 / (median(selMedDiffs) / fsample))
    pitchBpm2 = round(60 / (median(selDiffsVec) / fsample))
    cat(sprintf("pitchBpm = %g bpm, %g bpm \n", pitchBpm, pitchBpm2))
    
    plot(NA, xlim = c(0, saveLen), ylim = c(-6, 6))
    
    agr_det = 0
    for(id in 1:nrow(allHrSd))
    {
        if(is.na(allHrSd$sd[id]))
            next
        if (allHrSd$sd[id] > tresh_sd )
            next
        
        agr_det = agr_det + allDet[[id]]
        
        selSig = 1e9 * siglst[[id]][selRng]
        ftsig = filter_hr(selSig, fsample, 3, 30)
        ftsig = ftsig / sd(ftsig)
        lines(ftsig, type = 'l')
        
    }
    # filtr agr. fun
    if (F)
    {
        fltLP = butter(4, W = lp_sig / .5 / fsample, type = "low")
        agr_det = signal::filter(fltLP, agr_det)
    }
    nad = normAndDet(agr_det, startPitchBpm = pitchBpm, constPitch = TRUE, normOFF = T)
    selDetNorm = nad$selDetNorm
    maxs = nad$maxs
    rem_maxs = nad$rem_maxs
    
    plot(selDetNorm, type = 'l', col='red4')
    lines(maxs, selDetNorm[maxs], type="p", col = 'black', lwd = 2)
    lines(rem_maxs, selDetNorm[rem_maxs], type="p", col = 'red4')
    text(maxs+10, selDetNorm[maxs]+.5, labels = sprintf("(%d)", 1:length(maxs)), cex=.7)

    dsim = genSimSigAt(maxs, saveLen + fsample*2, peak_wd, fsample, sfq = pitchBpm / 60)
    dsim = dsim[1:saveLen]
    dsim = dsim / max(abs(dsim))
    dsim = pmax(0, dsim)
    stopifnot(!any(is.na(dsim)))

    plot(dsim, type = 'l', col = "red4", ylim = c(0,1))
    
    par(mfcol = c(1, 1))
    abline(v = maxs, col = "blue2", xpd = NA)
    
    # *******************************
    # if save ..
    # zapis inputu- i target-u dla TCF w Py
    # trzeba zwracać uwagę czy poprawny wybrany!
    # => i ew. zmienić fid !
    if(F)
    {
        tsig = data.frame(tsig = dsim)
        
        sc = trunc(range(selRng) / fsample)
        savPath = sprintf("%s_%d_%d_targ.Rdata", fname, sc[1], sc[2])
        savPath = sub(recDir, savDir, savPath)

        if (file.exists(savPath))
        {
            stop("*** plik target istnieje, możesz nadpisać ręcznie ***")
        }
        
        if(!dir.exists(dirname(savPath)))
            dir.create(dirname(savPath), recursive = T)
       
        save(tsig, file = savPath)
        cat(sprintf('TARGET: \n  "%s" \nSAVED!\n--\n', savPath))
        
        imgPath = sub("_targ.Rdata", "_img.png", savPath, ignore.case = T)
        dev.copy(png, imgPath, width = 1023, height = 768)
        dev.off()
        
        savPath = sprintf("%s_%d_%d_mx.Rdata", fname, sc[1], sc[2])
        savPath = sub(recDir, savDir, savPath)
        savName = sub(paste0(savDir, "/"), "", savPath)
        
        allSig = list()
        for(i in 1:length(siglst))
            allSig[[i]] = 1e9 * siglst[[i]][selRng]
        allSig = matrix(unlist(allSig), nrow=length(allSig[[1]]))
        dim = as.vector(dim(allSig))
        save(dim, allSig, file = savPath)

        cat(sprintf('INPUT: \n  "%s" \nSAVED!\n--\n', savPath))
        
        # upd targCompInfo.csv        
        info_path = file.path(savDir, "targCompInfo.csv")
        
        if(file.exists(info_path))
            targCompInfo = read.table(info_path, sep=",", header = T, stringsAsFactors = F)
        else
            targCompInfo = data.frame()
        infRow = data.frame(inpFile = savName, startSec= sc[1], fid = best_fid, bpm = pitchBpm)
        targCompInfo = rbind(targCompInfo, infRow)
        write.table(targCompInfo, info_path, sep=",", row.names = F)
        
        # sprawdzenie czy nie ma plików nieobecnych w csv    
        savedInp = dir(savDir, "*\\_mx.Rdata", recursive = T, ignore.case = T)
        f = savedInp[1]
        for(f in savedInp)
        {
            # f %in% targCompInfo$inpFile
            if (length(which(targCompInfo$inpFile == f)) != 1)
                print(sprintf("FILE '%s' not in info csv !!", f))
        }
    }

}


if(F)
    py_run_string("ld.Dispose()")






