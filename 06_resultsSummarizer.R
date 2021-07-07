rm(list=ls())
load("~/nBox/sero_switch_th/out/01outVAR.RData")
load("~/nBox/sero_switch_th/out/04_aHistDecomps.RData")
# load("~/nBox/sero_switch_th/out/02assesment.RData")

summarizerVar <- function(varOut){
  IRFS <- varOut$IRFS
  IRFhorizons <- length(varOut$IRFS[1,,1])
  IRFhorizons <- IRFhorizons/c(IRFhorizons,IRFhorizons/10,IRFhorizons/20,IRFhorizons/30,IRFhorizons/50,IRFhorizons/100)
  IRF_out     <- IRFS[,IRFhorizons,] 
  IRF_out_m     <- apply(IRF_out,MARGIN=c(2,3),median)
  IRF_out_q1     <- apply(IRF_out,MARGIN=c(2,3),quantile,probs=0.16)
  IRF_out_q3     <- apply(IRF_out,MARGIN=c(2,3),quantile,probs=0.84)
  IRF_out_notes  <- c("Rows are impulse response horizons in the time scale specified,
                      rows refer to the response of the dependent variable to a 1 unit shock on that Series,
                      M is median, q1 is 16tilde,q3 is 84 tilde")
  
  out <- list(IRF_out_m,
              IRF_out_q1,
              IRF_out_q3,
              IRF_out_notes)
  return(out)
  }

DecompOut <- function(Decomps=decomps){
  #do historical decomposition approximation only on second half
  DecompsLength <- c(nrow(Decomps[[1]])/2,nrow(Decomps[[1]]))
  #summarize decomposition of each series
  DecompsOut <- list()
  for (i in 1:length(Decomps)){
    temp <- Decomps[[i]]
    summarize <- colMeans(temp)
    DecompsOut[[i]] <- summarize
  }
  
  Decomp_out_notes <- c("First Four elements are median historical decompositions using median IRF, next four 16Tilde,next four 84Tilde,
                      Elements within these lists represent for first four: The mean historical decomposition of innovations for other
                      time series on the dependent variable. Next four elements: Dependent variable - historical decomposition (counterfactual scenario)")
  

    out<- list(DecompsOut,
              Decomp_out_notes)
    
    return(out)
}

summarizerVar1 <- summarizerVar(var1) 
summarizerVar2 <- summarizerVar(var2) 
summarizerVar3 <- summarizerVar(var3) 
summarizerVar4 <- summarizerVar(var4) 
decomps <- DecompOut(decomps)
