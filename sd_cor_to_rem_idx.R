
#  jeżeli pominę dany punkt i zmienność tacho zmaleje to usuwam
sd_cor_to_rem_idx = function(maxs)
{
    if(!length(maxs))
        return(c())
    sd_th = sd(diff(maxs)) * .8
    torem = c()
    for(im in 1:length(maxs))
    {
        n_sd = sd(diff(maxs[-im]))
        if(n_sd < sd_th)
        {
            torem = c(torem, im)
            print(paste0("sd_cor_to_rem_idx(): rem max ", im, ", n_sd=", n_sd, " sd_th=", sd_th))
        }
    }
    torem
}