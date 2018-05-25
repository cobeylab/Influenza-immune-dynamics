# Load the panel pomp object  --------------------------------------
if(file.exists(pomp_filename)){
  load(pomp_filename)
}

# Set the parameters -------------------------------------------------
guess.shared <- shared_params
for( i in c(1:length(inferred_params))){
  guess.shared[which(names(guess.shared) == names(inferred_params[i]))] <- inferred_params[i]
}
guess.specific <- panelObject@pParams$specific  

start <- Sys.time()
mf <- mif2(
  panelObject,
  Nmif = n_mif,
  shared.start = unlist(guess.shared),
  specific.start = matrix(
    data =  guess.specific,
    nrow = nrow(guess.specific),
    ncol = ncol(guess.specific),
    dimnames = list(rownames(guess.specific),
                    colnames(guess.specific))                 
  ),
  rw.sd = rw_sd_vec,
  cooling.type = "geometric",
  cooling.fraction.50 = cooling_rate,
  Np = n_particles
)
end <- Sys.time()
filename <- chain_filename
save(mf, file = filename)
print(end - start)
## Evaluate the likelihood ----------------------------------------------------------------
if(evaluate_Lhood == TRUE){
  ll <- logmeanexp(replicate(n_reps_pfilter,logLik(pfilter(mf,Np=n_particles_pfilter))),se=TRUE)
  output <- (data.frame(as.list(coef(mf)$shared),loglik=ll[1],loglik_se=ll[2], n_mif = n_mif_updated, n_part = n_particles_pfilter, chain = chainId))
  write.table(output, file = output_filename, sep = ",",col.names = FALSE, append=TRUE)
}



