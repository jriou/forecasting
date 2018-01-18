## ALL PREDICTIONS

args <- commandArgs(trailingOnly=TRUE)
library(dplyr);library(tidyr);library(rstan)
rstan_options(auto_write = TRUE) # multicore for cluster
options(mc.cores = parallel::detectCores())
load("datacluster_pred.Rdata")
load("island_specific_priors.Rdata")
load("region_specific_priors.Rdata")

# Functions
KL.div <- function (freqs1, freqs2) {
  # normalize
  freqs1 = freqs1/sum(freqs1) 
  freqs2 = freqs2/sum(freqs2)
  # replace 0 by small values
  freqs1[freqs1==0] = 1e-20
  freqs2[freqs2==0] = 1e-20
  # compute KL
  KL = sum(freqs1 * log(freqs1/freqs2))
  return(KL)
}

# Controls
dates <- data.frame(d=seq.Date(from=as.Date("2016-01-24"),to=as.Date("2017-02-19"),by=7))
Nchains <- 8
Nit <- 4000
Nwarmup <- 2000
Nthin <- 1
ResampleThreshold <- 1.05

for (sim in as.numeric(args[[1]]):as.numeric(args[[2]])) {
  datecut <- dates[sim,"d"]
  # compile stan model
  M_D1 = stan_model("T_D1.stan",model_name="M_D1")
  
  for(island in c("MARTINIQUE","GUADELOUPE","SAINTMARTIN")) {
    # extract data up to date
    D1_ISL <- filter(D1,ISLAND==island,DATE<=datecut)
    # extract future observations
    D1_FUT <- filter(D1,ISLAND==island,DATE>datecut)
    # limit to dates with at least 5 observations
    if(nrow(D1_ISL)>=5) {
      ## Models with island-specific priors --------------------------------------------------------------
      
      ## Non-informative prior (NIP)
      dataList_NIP <- list(
        # data
        W=nrow(D1_ISL), # number of records
        O_t=D1_ISL$NCASES, # number of reported cases
        sumO_t=D1_ISL$CUM_NCASES, # cumulative number of reported cases
        pop=unique(D1_ISL$POP), # island total population
        siL=5, # length of discretized serial interval distribution (in weeks)
        siW=dSI_ZIKV, # discretized serial interval distribution
        
        # priors
        P_R0_type=4, # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
        P_R0=c(1,0.2), # weakly-informative
        P_rho_type=1, # 1=beta; 2=uniform; 0=hyperprior
        P_rho=c(1,1), # weakly-informative
        
        # hyperpriors (only used for region-specific priors, type=0)
        P_mu_R0=c(1,1), # R0 mean (gamma distribution)
        P_sigma_R0=c(1,1), # R0 standard deviation (gamma distribution) 
        P_IL_mu_rho=c(0,1), # rho mean (inverse logit scale, normal distribution)
        P_IL_sigma_rho=c(1,1), # rho standard deviation (inverse logit scale, gamma distribution)
        
        # prediction
        pW=104 # number of weeks of prediction
      )
      
      S_D1_NIP <- sampling(object=M_D1, data=dataList_NIP, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                           control=list(adapt_delta=0.999))
      
      maxrhat <- max(summary(S_D1_NIP,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_NIP <- sampling(object=M_D1, data=dataList_NIP, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                             control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_NIP,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      
      ## Informative prior on R0 (IPR0)
      # island-specific
      dataList_IPR0 <- dataList_NIP
      dataList_IPR0$P_R0_type = 4 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPR0$P_R0 = P_R0[[island]]$estimate # informative
      
      S_D1_IPR0 <- sampling(object=M_D1, data=dataList_IPR0, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                            control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPR0,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPR0 <- sampling(object=M_D1, data=dataList_IPR0, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                              control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPR0,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      # region-specific
      dataList_IPR0_R <- dataList_IPR0
      dataList_IPR0_R$P_R0_type = 0 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPR0_R$P_mu_R0 = P_hyperpar_R0$mu$estimate # informative
      dataList_IPR0_R$P_sigma_R0 = P_hyperpar_R0$sigma$estimate # informative
      
      S_D1_IPR0_R <- sampling(object=M_D1, data=dataList_IPR0_R, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                              control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPR0_R,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPR0_R <- sampling(object=M_D1, data=dataList_IPR0_R, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                                control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPR0_R,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      
      ## Informative prior on rho (IPrho)
      # island-specific
      dataList_IPrho <- dataList_NIP
      dataList_IPrho$P_rho_type=1 # 1=beta; 2=uniform; 0=hyperprior
      dataList_IPrho$P_rho=P_rho[[island]]$estimate # informative
      
      S_D1_IPrho <- sampling(object=M_D1, data=dataList_IPrho, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                             control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPrho,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPrho <- sampling(object=M_D1, data=dataList_IPrho, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                               control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPrho,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      # region-specific
      dataList_IPrho_R <- dataList_IPrho
      dataList_IPrho_R$P_rho_type = 0 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPrho_R$P_IL_mu_rho = P_hyperpar_IL_rho$mu$estimate # informative
      dataList_IPrho_R$P_IL_sigma_rho = P_hyperpar_IL_rho$sigma$estimate # informative
      
      S_D1_IPrho_R <- sampling(object=M_D1, data=dataList_IPrho_R, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                               control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPrho_R,pars=c("rho","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPrho_R <- sampling(object=M_D1, data=dataList_IPrho_R, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                                 control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPrho_R,pars=c("rho","rho"))[[1]][,"Rhat"])
      }
      
      ## Informative prior on both (IPboth)
      # island-specific
      dataList_IPboth <- dataList_NIP
      dataList_IPboth$P_R0_type = 4 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPboth$P_R0 = P_R0[[island]]$estimate # informative
      dataList_IPboth$P_rho_type=1 # 1=beta; 2=uniform; 0=hyperprior
      dataList_IPboth$P_rho=P_rho[[island]]$estimate # informative
      
      S_D1_IPboth <- sampling(object=M_D1, data=dataList_IPboth, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                              control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPboth,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPboth <- sampling(object=M_D1, data=dataList_IPboth, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                                control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPboth,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      # region-specific
      dataList_IPboth_R <- dataList_IPboth
      dataList_IPboth_R$P_R0_type = 0 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPboth_R$P_mu_R0 = P_hyperpar_R0$mu$estimate # informative
      dataList_IPboth_R$P_sigma_R0 = P_hyperpar_R0$sigma$estimate # informative
      dataList_IPboth_R$P_rho_type = 0 # 1=exponential; 2=uniform; 3=normal; 4=gamma; 0=hyperprior
      dataList_IPboth_R$P_IL_mu_rho = P_hyperpar_IL_rho$mu$estimate # informative
      dataList_IPboth_R$P_IL_sigma_rho = P_hyperpar_IL_rho$sigma$estimate # informative
      
      S_D1_IPboth_R <- sampling(object=M_D1, data=dataList_IPboth_R, chains=Nchains , iter=Nit, warmup=Nwarmup, thin=Nthin,
                                control=list(adapt_delta=0.999))
      maxrhat <- max(summary(S_D1_IPboth_R,pars=c("R0","rho"))[[1]][,"Rhat"])
      nResample = 0
      while(maxrhat>ResampleThreshold & nResample<5) {
        nResample = nResample + 1
        S_D1_IPboth_R <- sampling(object=M_D1, data=dataList_IPboth_R, chains=Nchains, iter=Nit*(1+nResample/5), warmup=Nwarmup*(1+nResample/5), thin=Nthin,
                                  control=list(adapt_delta=0.999))
        maxrhat <- max(summary(S_D1_IPboth_R,pars=c("R0","rho"))[[1]][,"Rhat"])
      }
      
      
      ## Export relevant results
      S_ <- c("S_D1_NIP","S_D1_IPR0","S_D1_IPrho","S_D1_IPboth","S_D1_IPR0_R","S_D1_IPrho_R","S_D1_IPboth_R")
      NSAMPLE <- Nchains*(Nit-Nwarmup)/Nthin
      NOBS <- length(D1_ISL$NCASES)
      MINDATE <- min(D1_ISL$DATE)
      MAXDATE <- max(D1_ISL$DATE)
      REF <- D1_ISL$NCASES
      
      # posterior distribution of R0
      SUM_R0 <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="R0")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% tbl_df() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.)
      
      # posterior distribution of rho
      SUM_rho <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="rho")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% tbl_df() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.)
      
      # fitted values
      SUM_E_NCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) {
        FIT <- summary(x,pars="pred_lp")[[1]]
        FIT <- cbind(dplyr::select(D1_ISL,DATE,N_WEEK,NCASES),FIT)
      }) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7*NOBS),ISLAND=rep(island,7*NOBS),
                             PRIOR_R0=c(rep(0,NOBS),rep(1,NOBS),rep(0,NOBS),rep(1,NOBS),rep(2,NOBS),rep(0,NOBS),rep(2,NOBS)),
                             PRIOR_rho=c(rep(0,NOBS),rep(0,NOBS),rep(1,NOBS),rep(1,NOBS),rep(0,NOBS),rep(2,NOBS),rep(2,NOBS))),.) %>% tbl_df()
      
      # estimated number of total cases until date of analysis
      SUM_E_OBSERVEDCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="totlp")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% tbl_df() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.)
      
      # estimated number of total cases until date of analysis
      SUM_E_TOTALCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="totoverall")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% tbl_df() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.)
      
      # estimated attack rate until date of analysis
      SUM_E_AR <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="attackrate")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% tbl_df() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.)
      
      # forecasted evolution of the number of observed cases
      SUM_F_NCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) {
        PRED <- extract(x,pars="pO")[[1]] 
        MEAN <- apply(PRED,2,mean)
        SD <- apply(PRED,2,sd)
        MEDIAN <- apply(PRED,2,median)
        DELTA <- abs(sweep(PRED,2,colMeans(PRED)))
        RANK <- apply(DELTA,2,rank,ties.method="random")
        RANK <- rank(apply(RANK,1,max),ties.method="random")
        KEEP.95 <- which(RANK<=length(RANK)*0.95)
        KEEP.50 <- which(RANK<=length(RANK)*0.5)
        PRED2.95 <- PRED[KEEP.95,]
        PRED2.50 <- PRED[KEEP.50,]
        INTERVAL.95 <- t(apply(PRED2.95,2,function(y) c(min(y),max(y))))
        INTERVAL.50 <- t(apply(PRED2.50,2,function(y) c(min(y),max(y))))
        DATE <- MAXDATE+1:104*7
        OUT <- data.frame(mean=MEAN,se_mean=SD/sqrt(dim(PRED)[2]),sd=SD,`2.5%`=INTERVAL.95[,1],`25%`=INTERVAL.50[,1],`50%`=MEDIAN,`75%`=INTERVAL.50[,2],`97.5%`=INTERVAL.95[,2],n_eff=NA,Rhat=NA,DATE=DATE)
        names(OUT)[4:8] <- c("2.5%","25%","50%","75%","97.5%")
        return(OUT)
      }) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>%
        bind_cols(data.frame(DATECUT=rep(datecut,104*7),ISLAND=rep(island,104*7),
                             PRIOR_R0=c(rep(0,104),rep(1,104),rep(0,104),rep(1,104),rep(2,104),rep(0,104),rep(2,104)),
                             PRIOR_rho=c(rep(0,104),rep(0,104),rep(1,104),rep(1,104),rep(0,104),rep(2,104),rep(2,104))),.) %>%
        tbl_df()
      # bind cases observed after analysis
      if(nrow(D1_FUT)>1) {
        SUM_F_NCASES = left_join(SUM_F_NCASES,select(D1_FUT,DATE,NCASES)) 
      } else {
        SUM_F_NCASES$NCASES = 0
      }
      
      # RMSD between forecast and actual observed cases
      SUM_RMSD <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) {
        REF <- filter(D1,ISLAND==island,DATE>datecut)$NCASES
        L <- length(REF)
        if(L>4) {
          PRED <- extract(x,pars="pO")[[1]][,1:L]
          RMSD <- apply(PRED,1,function(y) sqrt(mean((y-REF)^2)))
          RMSD <- c(mean=mean(RMSD),se_mean=sd(RMSD)/sqrt(length(RMSD)),sd=sd(RMSD),
                    quantile(RMSD,probs=c(0.025,0.25,0.5,0.75,0.975)),n_eff=NA,Rhat=NA)
        } else {
          RMSD <- c(mean=NA,se_mean=NA,sd=NA,
                    rep(NA,5),n_eff=NA,Rhat=NA)
          names(RMSD) = c("mean","se_mean","sd","2.5%","25%","50%","75%","97.5%","n_eff","Rhat")
        }
      }) %>% 
        do.call("rbind",.) %>%
        as.data.frame()  %>%
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      
      # forecasted final number of observed cases
      SUM_F_OBSERVEDCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="ptotlp")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      # forecasted final number of total cases
      SUM_F_TOTALCASES <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="ptotoverall")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      # forecasted final attack rate
      SUM_F_AR <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) summary(x,pars="pattackrate")[[1]]) %>%
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      # forecasted date of peak incidence
      SUM_F_PEAK <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) {
        PRED <- cbind(extract(x,pars="pred_lp")[[1]],extract(x,pars="pO")[[1]])
        PEAK <- unlist(apply(PRED,1,function(y) which(y==max(y))))
        PEAK <- round(c(mean=mean(PEAK),se_mean=sd(PEAK)/length(PEAK),sd=sd(PEAK),quantile(PEAK,probs=c(0.025,0.25,0.5,0.75,0.975)),n_eff=NA,Rhat=NA))
      }) %>% 
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      # forecasted maximal incidence of clinical cases
      SUM_F_MAXINC <- lapply(list(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth,S_D1_IPR0_R,S_D1_IPrho_R,S_D1_IPboth_R),function(x) {
        PRED <- cbind(extract(x,pars="pred_lp")[[1]],extract(x,pars="pO")[[1]])
        PEAK <- unlist(apply(PRED,1,function(y) max(y)))
        PEAK <- round(c(mean=mean(PEAK),se_mean=sd(PEAK)/length(PEAK),sd=sd(PEAK),quantile(PEAK,probs=c(0.025,0.25,0.5,0.75,0.975)),n_eff=NA,Rhat=NA))
      }) %>% 
        do.call("rbind",.) %>%
        as.data.frame() %>% 
        bind_cols(data.frame(DATECUT=rep(datecut,7),ISLAND=rep(island,7),PRIOR_R0=c(0,1,0,1,2,0,2),PRIOR_rho=c(0,0,1,1,0,2,2)),.) %>% tbl_df()
      
      # KL distance
      SUM_KLD <- NULL
      
      # NIP ***
      R0_Q <- density(extract(S_D1_NIP,pars="R0")[[1]],from=0,to=10)
      if(dataList_NIP$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_NIP$P_R0[[1]],rate=dataList_NIP$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_NIP,pars="rho")[[1]],from=0,to=1)
      if(dataList_NIP$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_NIP$P_rho[[1]],shape2=dataList_NIP$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=0,PRIOR_rho=0,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPR0 ***
      R0_Q <- density(extract(S_D1_IPR0,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPR0$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPR0$P_R0[[1]],rate=dataList_IPR0$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPR0,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPR0$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPR0$P_rho[[1]],shape2=dataList_IPR0$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=1,PRIOR_rho=0,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPrho ***
      R0_Q <- density(extract(S_D1_IPrho,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPrho$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPrho$P_R0[[1]],rate=dataList_IPrho$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPrho,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPrho$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPrho$P_rho[[1]],shape2=dataList_IPrho$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=0,PRIOR_rho=1,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPboth ***
      R0_Q <- density(extract(S_D1_IPboth,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPboth$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPboth$P_R0[[1]],rate=dataList_IPboth$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPboth,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPboth$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPboth$P_rho[[1]],shape2=dataList_IPboth$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=1,PRIOR_rho=1,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPR0_R ***
      R0_Q <- density(extract(S_D1_IPR0_R,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPR0_R$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPR0_R$P_R0[[1]],rate=dataList_IPR0_R$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPR0_R,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPR0_R$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPR0_R$P_rho[[1]],shape2=dataList_IPR0_R$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=2,PRIOR_rho=0,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPrho_R ***
      R0_Q <- density(extract(S_D1_IPrho_R,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPrho_R$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPrho_R$P_R0[[1]],rate=dataList_IPrho_R$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPrho_R,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPrho_R$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPrho_R$P_rho[[1]],shape2=dataList_IPrho_R$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=0,PRIOR_rho=2,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      # IPboth_R ***
      R0_Q <- density(extract(S_D1_IPboth_R,pars="R0")[[1]],from=0,to=10)
      if(dataList_IPboth_R$P_R0_type==4) R0_P <- list(x=R0_Q$x,y=dgamma(R0_Q$x,shape=dataList_IPboth_R$P_R0[[1]],rate=dataList_IPboth_R$P_R0[[2]]))
      # plot(R0_Q$y ~ R0_Q$x,type="l",col="red");points(R0_P$y ~ R0_P$x,type="l")
      rho_Q <- density(extract(S_D1_IPboth_R,pars="rho")[[1]],from=0,to=1)
      if(dataList_IPboth_R$P_rho_type==1) rho_P <- list(x=rho_Q$x,y=dbeta(rho_Q$x,shape1=dataList_IPboth_R$P_rho[[1]],shape2=dataList_IPboth_R$P_rho[[2]]))
      # plot(rho_Q$y ~ rho_Q$x,type="l",col="red");points(rho_P$y ~ rho_P$x,type="l")
      SUM_KLD <- rbind(SUM_KLD,
                       data.frame(DATECUT=datecut,ISLAND=island,PRIOR_R0=2,PRIOR_rho=2,
                                  KLD_R0=KL.div(R0_P$y,R0_Q$y),
                                  KLD_rho=KL.div(rho_P$y,rho_Q$y)))
      
      
      save(S_D1_IPboth,S_D1_IPboth_R,S_D1_IPR0,S_D1_IPR0_R,S_D1_IPrho,S_D1_IPrho_R,S_D1_NIP,
           file=paste0("~/datashortcut/Sampling_",island,"_",as.character(datecut),"_15mai17.Rdata"))
      save(SUM_E_AR,SUM_E_NCASES,SUM_E_OBSERVEDCASES,SUM_E_TOTALCASES,SUM_F_AR,SUM_F_MAXINC,SUM_F_NCASES,
           SUM_F_OBSERVEDCASES,SUM_F_PEAK,SUM_F_TOTALCASES,SUM_KLD,SUM_R0,SUM_rho,SUM_RMSD,
           file=paste0("~/datashortcut/All_predictions_",island,"_",as.character(datecut),"_15mai17.Rdata"))
      rm(S_D1_NIP,S_D1_IPR0,S_D1_IPrho,S_D1_IPboth)
    } 
  }
}
