
require(dplyr);require(tidyr);require(ggplot2);library(cowplot)
require(rstan);require(coda);require(loo)
library(xtable)
library(MASS);select <- dplyr::select 
library(LaplacesDemon)
library(scales)
library(tibble)

# IMPORTATION & DATA-MANAGEMENT  ----------------------------------------------------------------------------------------------
key <- data.frame(ISLAND=c("TAHITI","MOOREA","SLV","TUAMOTU","MARQUISES","AUSTRALES","SAINTMARTIN","MARTINIQUE","GUADELOUPE"),
                  ISLAND2=c("Tahiti","Mo'orea","Sous-le-vent Isl.","Tuamotus","Marquesas Isl.","Austral Isl.","Saint-Martin","Martinique","Guadeloupe"),
                  ISLAND_ID=1:9,
                  ISLAND_ABB=c("TAH","MOO","SLV","TUA","MRQ","AUS","MAF","MTQ","GLP"),
                  POP=c(184000,16000,33000,16000,9000,7000,36000,385000,400000),
                  REGION=c("FP","FP","FP","FP","FP","FP","CA","CA","CA"),
                  REGION_ID=c(rep(1,6),rep(2,3)),stringsAsFactors=FALSE)
key$ISLAND2 <- factor(key$ISLAND2,levels=c("Austral Isl.","Mo'orea","Marquesas Isl.","Sous-le-vent Isl.","Tahiti","Tuamotus",
                                           "Guadeloupe","Martinique","Saint-Martin"))
zikachik <- tbl_df(as.data.frame(read.csv2("data/datazikachik6.csv",header=T,sep=",",stringsAsFactors=FALSE))) %>%
  filter(!ISLAND %in% c("PORTORICO","SAINTBARTH","GUYANE")) %>%
  left_join(key) %>%
  mutate(VIRUS_ID=ifelse(VIRUS=="CHIKV",0,1),
         DATE=as.Date(paste(YEAR,WEEK,"1",sep="-"),format="%Y-%U-%u")-1) %>%
  group_by(ISLAND,VIRUS) %>%
  mutate(N_WEEK=row_number()-1) %>%
  # filter(!(ISLAND=="SAINTMARTIN" & N_WEEK>25 & VIRUS=="CHIKV")) %>%
  dplyr::select(REGION,ISLAND,VIRUS,DATE,N_WEEK,NCASES,POP,REGION_ID,ISLAND2,ISLAND_ABB,ISLAND_ID,VIRUS_ID,WEEK,YEAR) %>%
  mutate(NCASES_M1=ifelse(is.na(lag(NCASES,1)),0,lag(NCASES,1)),
         NCASES_M2=ifelse(is.na(lag(NCASES,2)),0,lag(NCASES,2)),
         NCASES_M3=ifelse(is.na(lag(NCASES,3)),0,lag(NCASES,3)),
         NCASES_M4=ifelse(is.na(lag(NCASES,4)),0,lag(NCASES,4)),
         NCASES_M5=ifelse(is.na(lag(NCASES,5)),0,lag(NCASES,5)),
         CUM_NCASES=cumsum(NCASES_M1),
         INC=NCASES/POP,
         CUMINC=cumsum(NCASES)/POP,
         TOT_NCASES=sum(NCASES),
         AR=TOT_NCASES/POP) %>%
  ungroup()

ggplot(zikachik) +
  geom_area(aes(x=DATE,y=NCASES,fill=VIRUS)) +
  facet_wrap(~ISLAND2,scales="free_y",ncol=2,dir="v") +
  scale_fill_manual(values=c(colCHIKV,colZIKV)) +
  theme_bw()


# COMPUTE Ostar ----------------------------------------------------------------------------------------------
SI_ZIKV <- c(2.5,0.7)
SI_CHIKV <- c(1.6,0.6)

dSI_ZIKV = pgamma(1:5+.5,shape=SI_ZIKV[1]^2/SI_ZIKV[2]^2,rate=SI_ZIKV[1]/SI_ZIKV[2]^2)-
  pgamma(1:5-.5,shape=SI_ZIKV[1]^2/SI_ZIKV[2]^2,rate=SI_ZIKV[1]/SI_ZIKV[2]^2)
dSI_CHIKV = pgamma(1:5+.5,shape=SI_CHIKV[1]^2/SI_CHIKV[2]^2,rate=SI_CHIKV[1]/SI_CHIKV[2]^2)-
  pgamma(1:5-.5,shape=SI_CHIKV[1]^2/SI_CHIKV[2]^2,rate=SI_CHIKV[1]/SI_CHIKV[2]^2)

zikachik <- zikachik %>%
  mutate(SI_mean=ifelse(VIRUS=="ZIKV",SI_ZIKV[1],SI_CHIKV[1]),
         SI_sd=ifelse(VIRUS=="ZIKV",SI_ZIKV[2],SI_CHIKV[2]))

Ostar <- NULL
for (i in 1:dim(zikachik)[1]) {
  m = unlist(zikachik[i,"SI_mean"])
  s = unlist(zikachik[i,"SI_sd"])
  w = pgamma(1:5+.5,shape=m^2/s^2,rate=m/s^2)-pgamma(1:5-.5,shape=m^2/s^2,rate=m/s^2)
  w = w/sum(w)
  Ostar = c(Ostar,sum(zikachik[i,c("NCASES_M1","NCASES_M2","NCASES_M3","NCASES_M4","NCASES_M5")]*w))
}
zikachik$Ostar <- Ostar

D1 <- filter(zikachik,REGION_ID==2,VIRUS=="ZIKV",Ostar>0)
D2 <- filter(zikachik,REGION_ID==2,VIRUS=="CHIKV",Ostar>0)
D3 <- filter(zikachik,REGION_ID==1,Ostar>0)
save(zikachik,key,D1,D2,D3,dSI_ZIKV,file="datacluster_pred.Rdata")


# COMPUTE PRIORS ----------------------------------------------------------------------------------------------

library(dplyr);library(tidyr);library(rstan)
rstan_options(auto_write = TRUE) # multicore for cluster
options(mc.cores = parallel::detectCores())
load("datacluster_pred.Rdata")

Nchains <- 8
Nit <- 4000
Nwarmup <- 2000
Nthin <- 1


# Priors on R0 and rho for CHIKV in Martinique, Guadeloupe, Saint Martin & Guyane 

M_D2 <- stan_model(file="T_D2.stan",model_name="M_D2")

dataList_D2 <- list(
  W=length(D2$NCASES),
  K=length(unique(D2$ISLAND_ID)),
  island=as.numeric(factor(D2$ISLAND_ID,labels=paste0(1:length(unique(D2$ISLAND_ID))))),
  O_t=D2$NCASES,
  Ostar=D2$Ostar,
  sumO_t=D2$CUM_NCASES,
  pop=tapply(D2$POP,D2$ISLAND_ID,unique)
)

S_D2 <- sampling(M_D2,data=dataList_D2,chains=Nchains,iter=Nit,warmup=Nwarmup,thin=Nthin,control=list(adapt_delta=0.999))
intpars = c("mu_R0_CHIKV","sigma_R0_CHIKV","mu_rho_CHIKV","IL_sigma_rho_CHIKV","phi","R0_CHIKV","rho_CHIKV")
summary(S_D2,pars=intpars)[[1]]
stan_trace(S_D2,pars=intpars)
pairs(S_D2,pars=intpars[1:5])
pairs(S_D2,pars=intpars[c(2,6)])
pairs(S_D2,pars=intpars[c(4,7)])

save(S_D2,file="S_D2.Rdata")



# Priors on the transmission and reporting ratios of ZIKV compared to CHIKV in French Polynesia 

M_D3 <- stan_model(file="T_D3.stan",model_name="M_D3")

dataList_D3 <- list(
  W=length(D3$NCASES),
  K=length(unique(D3$ISLAND)),
  J=2,
  island=D3$ISLAND_ID,
  virus=D3$VIRUS_ID,
  O_t=D3$NCASES,
  Ostar=D3$Ostar,
  sumO_t=D3$CUM_NCASES,
  pop=tapply(D3$POP,D3$ISLAND_ID,unique)
)

S_D3 <- sampling(M_D3,data=dataList_D3,chains=Nchains,iter=Nit,warmup=Nwarmup,thin=Nthin,control=list(adapt_delta=0.999))
intpars = c("mu_R0_CHIKV","sigma_R0_CHIKV","beta_R0","mu_rho_CHIKV","IL_sigma_rho_CHIKV","beta_rho","phi","R0_CHIKV","rho_CHIKV")
summary(S_D3,pars=intpars)[[1]]
stan_trace(S_D3,pars=intpars[1:4])
pairs(S_D3,pars=intpars[1:5])
pairs(S_D3,pars=intpars[c(2,6)])
pairs(S_D3,pars=intpars[c(4,7)])

save(S_D3,file="S_D3.Rdata")


# Combine priors

load("S_D2.Rdata")
load("S_D3.Rdata")

# Island-specific priors
samples_beta_R0_D3 =  extract(S_D3,pars="beta_R0")[[1]]
samples_beta_rho_D3 =  extract(S_D3,pars="beta_rho")[[1]]

samples_R0_CHIKV_D2 = extract(S_D2,pars="R0_CHIKV")[[1]] %>%
  split(., rep(1:ncol(.), each = nrow(.)))
names(samples_R0_CHIKV_D2) = c("SAINTMARTIN","MARTINIQUE","GUADELOUPE")
samples_rho_CHIKV_D2 = extract(S_D2,pars="rho_CHIKV")[[1]] %>%
  split(., rep(1:ncol(.), each = nrow(.)))
names(samples_rho_CHIKV_D2) = c("SAINTMARTIN","MARTINIQUE","GUADELOUPE")

samples_R0_CHIKV_D23 = lapply(samples_R0_CHIKV_D2,function(x) x*samples_beta_R0_D3)
samples_rho_CHIKV_D2 = lapply(samples_rho_CHIKV_D2,function(x) x*samples_beta_rho_D3)

P_R0 = lapply(samples_R0_CHIKV_D23,fitdistr,dgamma,start=list(shape=1,rate=1))
P_rho = lapply(samples_rho_CHIKV_D2,fitdistr,dbeta,start=list(shape1=1,shape2=1))


# Region-specific hyperparameter priors
samples_mu_R0_CHIKV_D2 = extract(S_D2,pars="mu_R0_CHIKV")[[1]]
samples_sigma_R0_CHIKV_D2 = extract(S_D2,pars="sigma_R0_CHIKV")[[1]]
samples_mu_R0_CHIKV_D23 = samples_mu_R0_CHIKV_D2 * samples_beta_R0_D3

P_hyperpar_R0 = list(mu=fitdistr(samples_mu_R0_CHIKV_D23,dgamma,start=list(shape=1,rate=0.1)),
                     sigma=fitdistr(samples_sigma_R0_CHIKV_D2,dgamma,start=list(shape=1,rate=1)))


samples_IL_mu_rho_CHIKV_D2 = extract(S_D2,pars="IL_mu_rho_CHIKV")[[1]]
samples_IL_sigma_rho_CHIKV_D2 = extract(S_D2,pars="IL_sigma_rho_CHIKV")[[1]]
samples_IL_mu_rho_CHIKV_D23 = logit(invlogit(samples_IL_mu_rho_CHIKV_D2) * samples_beta_rho_D3)

P_hyperpar_IL_rho = list(mu=fitdistr(samples_IL_mu_rho_CHIKV_D23,dnorm,start=list(mean=0,sd=1)),
                         sigma=fitdistr(samples_IL_sigma_rho_CHIKV_D2,dgamma,start=list(shape=1,rate=1)))


#  Compute prior densities for computing KL-distance

getdens = function(x,from,to) {
  d = density(x,from=from,to=to)
  d$y[d$y==0] = 1e-20
  return(d)
} 

dens_R0_CHIKV_D23 = lapply(samples_R0_CHIKV_D23,getdens,from=0,to=10)
dens_rho_CHIKV_D2 = lapply(samples_R0_CHIKV_D23,getdens,from=0,to=1)

save(P_R0,P_rho,
     dens_R0_CHIKV_D23,dens_rho_CHIKV_D2,
     file="island_specific_priors.Rdata")

dens_mu_R0_CHIKV_D23 = getdens(samples_mu_R0_CHIKV_D23,from=0,to=10)
dens_sigma_R0_CHIKV_D2 = getdens(samples_sigma_R0_CHIKV_D2,from=0,to=10)
dens_IL_mu_rho_CHIKV_D23 = getdens(samples_IL_mu_rho_CHIKV_D23,from=-10,to=10)
dens_IL_sigma_rho_CHIKV_D2 = getdens(samples_IL_sigma_rho_CHIKV_D2,from=0,to=10)

save(P_hyperpar_R0,P_hyperpar_IL_rho,
     dens_mu_R0_CHIKV_D23,dens_sigma_R0_CHIKV_D2,
     dens_IL_mu_rho_CHIKV_D23,dens_IL_sigma_rho_CHIKV_D2,file="region_specific_priors.Rdata")


# text prior
lapply(P_R0,function(x) c(mean=x$estimate[1]/x$estimate[2],
                          q2.5=qgamma(0.025,x$estimate[1],x$estimate[2]),
                          q97.5=qgamma(0.975,x$estimate[1],x$estimate[2])))
lapply(P_rho,function(x) c(mean=x$estimate[1]/(x$estimate[1]+x$estimate[2]),
                          q2.5=qbeta(0.025,x$estimate[1],x$estimate[2]),
                          q97.5=qbeta(0.975,x$estimate[1],x$estimate[2])))




# IMPORT RESULTS (on cluster) ----------------------------------------------------------------------------------------------
library(dplyr)
dates <- data.frame(d=seq.Date(from=as.Date("2016-01-24"),to=as.Date("2017-02-19"),by=7))
TSUM_E_AR <- NULL
TSUM_E_NCASES <- NULL
TSUM_E_OBSERVEDCASES <- NULL
TSUM_E_TOTALCASES <- NULL
TSUM_F_AR <- NULL
TSUM_F_MAXINC <- NULL
TSUM_F_NCASES <- NULL
TSUM_F_OBSERVEDCASES <- NULL
TSUM_F_PEAK <- NULL
TSUM_F_TOTALCASES <- NULL
TSUM_KLD <- NULL
TSUM_R0 <- NULL
TSUM_rho <- NULL
TSUM_RMSD <- NULL

for(h in 1:dim(dates)[[1]]) {
  for(island in c("MARTINIQUE","SAINTMARTIN","GUADELOUPE")) {
    path <- paste0("All_predictions_",island,"_",dates[h,"d"],"_15mai17.Rdata")
    if(file.exists(path)) {
      l <- load(path)
      TSUM_E_AR <-rbind(TSUM_E_AR ,SUM_E_AR )
      TSUM_E_NCASES <-rbind( TSUM_E_NCASES,SUM_E_NCASES )
      TSUM_E_OBSERVEDCASES <-rbind( TSUM_E_OBSERVEDCASES,SUM_E_OBSERVEDCASES )
      TSUM_E_TOTALCASES <-rbind(TSUM_E_TOTALCASES , SUM_E_TOTALCASES)
      TSUM_F_AR <-rbind( TSUM_F_AR,SUM_F_AR)
      TSUM_F_MAXINC <-rbind( TSUM_F_MAXINC,SUM_F_MAXINC )
      TSUM_F_NCASES <-rbind(TSUM_F_NCASES ,SUM_F_NCASES )
      TSUM_F_OBSERVEDCASES <-rbind( TSUM_F_OBSERVEDCASES,SUM_F_OBSERVEDCASES )
      TSUM_F_PEAK <-rbind(TSUM_F_PEAK ,SUM_F_PEAK )
      TSUM_F_TOTALCASES <-rbind( TSUM_F_TOTALCASES,SUM_F_TOTALCASES )
      TSUM_KLD <-rbind( TSUM_KLD,SUM_KLD )
      TSUM_R0 <-rbind(TSUM_R0 ,SUM_R0)
      TSUM_rho <-rbind(TSUM_rho ,SUM_rho)
      if(h==1 | identical(names(TSUM_RMSD),names(SUM_RMSD))) TSUM_RMSD <- rbind(TSUM_RMSD ,SUM_RMSD)
    }
  }
  cat(".")
}
save(TSUM_E_AR,
     TSUM_E_NCASES,
     TSUM_E_OBSERVEDCASES,
     TSUM_E_TOTALCASES,
     TSUM_F_AR,
     TSUM_F_MAXINC,
     TSUM_F_NCASES,
     TSUM_F_OBSERVEDCASES,
     TSUM_F_PEAK,
     TSUM_F_TOTALCASES,
     TSUM_KLD,
     TSUM_R0,
     TSUM_rho,
     TSUM_RMSD,file="Outcome_Prediction_15mai17.Rdata")


# OUTPUT ----------------------------------------------------------------------------------------------------------------

require(dplyr);require(tidyr);require(ggplot2);library(cowplot);library(grid);library(scales);library(png)
require(rstan);require(coda);require(loo)
library(xtable)
library(MASS);select <- dplyr::select 
library(gridExtra)
library(LaplacesDemon)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
load("datacluster_pred.Rdata")
load("island_specific_priors.Rdata")
load("region_specific_priors.Rdata")
l = load("Outcome_Prediction_15mai17.Rdata")

dates <- data.frame(d=seq.Date(from=as.Date("2016-01-31"),to=as.Date("2016-10-02"),by=7))
island_labels <- c(MARTINIQUE="Martinique",SAINTMARTIN="Saint-Martin",GUADELOUPE="Guadeloupe")
priors_labels <-c(
  IPR0 = bquote("IP on R"["0"]),
  NIPR0 = bquote("N-IP on R"["0"]),
  IPrho = bquote("IP on \u03C1"),
  NIPrho = bquote("N-IP on \u03C1")
)

Sys.setlocale("LC_TIME", "English_United States")
windowsFonts(Arial=windowsFont("TT Arial")) 

theme_bw2 <- function(base_size = 12, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(          axis.title.y=element_text(size=8,angle=90,margin=margin(0,5,0,0)),
                    axis.title.x=element_text(size=8),
                    axis.ticks=element_line(size=.2),
                    axis.line = element_line(colour = "black"),
                    plot.title=element_text(size=9,vjust=1),
                    strip.background=element_rect(fill=NA,colour=NA,size=.5),
                    strip.text=element_text(size=9),
                    panel.background = element_rect(fill=NA,colour=NA),
                    panel.grid = element_blank(),
                    panel.border=element_rect(colour="black",fill=NA),
                    axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)
    )
}

rac = "~/Nextcloud/zikachik_shared/manuscript_prediction/figures/"
rac = "D:/julien/Nextcloud/zikachik_shared/manuscript_prediction/figures/"



# legend
colZIKV <- "steelblue"
colCHIKV <- "tomato"
colPRIOR = data.frame(APPROACH=factor(c("Non-informative\nprior distributions\t","Island-specific\nprior distributions\t","Region-specific\nprior distributions\t"),
                                      levels=c("Non-informative\nprior distributions\t","Island-specific\nprior distributions\t","Region-specific\nprior distributions\t")),
                      PRIOR_R0=c(0,1,2),
                      PRIOR_rho=c(0,1,2),
                      APPROACH2=factor(c("NI","L","R"),levels=c("NI","R","L")),
                      col=as.character(c("#009E73","#D55E00","purple3")),
                      stringsAsFactors = F)
colPRIOR2 = data.frame(APPROACH=factor(c("Non-informative\t","Regional\t","Local\t"),
                                      levels=c("Non-informative\t","Regional\t","Local\t")),
                      PRIOR_R0=c(0,1,2),
                      PRIOR_rho=c(0,1,2),
                      col=as.character(c("#009E73","purple3","#D55E00")),
                      stringsAsFactors = F)
Fig_leg = 
ggplot(colPRIOR2) +
  geom_line(aes(x=PRIOR_R0,y=PRIOR_rho,colour=APPROACH),size=1) +
  scale_colour_manual(values=colPRIOR2$col) +
  labs(colour="Choice of prior:\t") +
  theme_bw2() +
  theme(legend.title.align=0.5,
        legend.title=element_text(size=10,vjust=0,face="italic"),
        legend.text=element_text(size=10,lineheight=.9,vjust=0),
        legend.key.height = unit(.8,"cm"),
        legend.key.width = unit(1,"cm"))
Fig_legV = get_legend(Fig_leg)
Fig_legH = get_legend(Fig_leg+theme(legend.position="bottom"))

# timings
timing = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
                    ISLAND2=c("Guadeloupe","Martinique","Saint-Martin"),
                    POP=c(400000,385000,36000),
                    start=as.Date(c("2016-04-10","2016-01-10","2016-02-21")),
                    peak=c("2016-24","2016-12","2016-29"),
                    end=c("2016-42","2016-39","2017-4"),
                    top=c(3600,3600,360),
                    min=as.Date("2015-12-20"),
                    max=as.Date("2017-02-19"),
                    labS=c(3250,3250,325),
                    labP=c(3250,3250,325),
                    labE=c(3250,3250,325),stringsAsFactors = FALSE) %>%
  mutate(peak=as.Date(paste(peak,"1",sep="-"),format="%Y-%U-%u")-1,
         end=as.Date(paste(end,"1",sep="-"),format="%Y-%U-%u")-1)

fullD1 = filter(zikachik,VIRUS=="ZIKV",REGION=="CA") %>%
  mutate(ISLAND2=as.character(ISLAND2))

##################################################################################################################
# Figure 1 with phases 

modt =  theme(strip.text=element_text(size=10,margin = margin(0,0,0,0, "cm")),
              axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title.y=element_text(size=9),
              axis.title.x=element_text(size=9))

marg = .035
timing$threshold = c(200,200,20)
timing = left_join(fullD1,timing) %>%
  filter(DATE==start) %>%
  select(ISLAND,ISLAND2,NCASES.start=NCASES) %>%
  left_join(timing)
timing = left_join(fullD1,timing) %>%
  filter(DATE==peak) %>%
  select(ISLAND,ISLAND2,NCASES.peak=NCASES) %>%
  left_join(timing)
timing = left_join(fullD1,timing) %>%
  filter(DATE==end) %>%
  select(ISLAND,ISLAND2,NCASES.end=NCASES) %>%
  left_join(timing)

Fig1_A = ggplot() +
  geom_bar(data=fullD1,aes(x=DATE,y=NCASES),fill=colZIKV,stat="identity") +
  geom_point(data=timing,aes(x=start-20,y=top),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=NCASES.start+top*marg,yend=labS),alpha=.6) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=NCASES.peak+top*marg,yend=labP),alpha=.6) +
  geom_segment(data=timing,aes(x=end-1,xend=end-1,y=NCASES.end+top*marg,yend=labE),alpha=.6) +
  
  geom_point(data=timing,aes(x=start-1,y=NCASES.start+top*marg),shape=25,fill="black",size=1.2,alpha=.6) +
  geom_point(data=timing,aes(x=peak-1,y=NCASES.peak+top*marg),shape=25,fill="black",size=1.2,alpha=.6) +
  geom_point(data=timing,aes(x=end-1,y=NCASES.end+top*marg),shape=25,fill="black",size=1.2,alpha=.6) +
  
  geom_label(data=timing,aes(x=start-1,y=labS),label="S",size=3) +
  geom_label(data=timing,aes(x=peak-1,y=labP),label="P",size=3) +
  geom_label(data=timing,aes(x=end-1,y=labE),label="E",size=3) +
  
  geom_hline(data=timing,aes(yintercept=threshold),colour="red",linetype=2,size=.4) +
  
  facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
  theme_bw2() + 
  labs(x="",y="Incident cases (N)") +
  scale_x_date(breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
               labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
  scale_y_continuous(expand=c(0,0),labels=function(x) ifelse(x<1000,paste(x),paste0(x/1000,"K"))) +
  modt

Fig1_A

#ggsave(plot=Fig_datatiming, 
 #      filename=paste0(rac,"Fig1.pdf"),
  #     width=16,height=6,unit="cm")


tt1 = filter(zikachik,REGION=="CA",VIRUS=="CHIKV") %>%
  filter(DATE>="2013-12-20",DATE<="2015-02-01")
timing$chik = c(7000,7000,700)
timing$ISLAND2=c("Guadeloupe","Martinique","Saint-Martin")
Fig1_B1 = ggplot() +
  geom_bar(data=tt1,aes(x=DATE,y=NCASES),fill=colCHIKV,stat="identity") +
  geom_point(data=timing,aes(x=as.Date("2013-12-01"),y=chik),colour=NA) +
  geom_point(data=timing,aes(x=as.Date("2015-02-20"),y=chik),colour=NA) +
  facet_wrap(~ISLAND2,scales="free_y",ncol=1) +
  scale_y_continuous(expand=c(0,0),#breaks=c(0,10,20),#limits=c(0,22),
                     labels=function(x) ifelse(x<1000,paste(x),paste0(x/1000,"K"))) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2014-01-01","2015-01-01")),
               labels=c("2014","2015")) +
  theme_bw2() +
  labs(x="",y="Incident cases (N)") +
  modt + theme(strip.text=element_text(size=8),axis.text.y=element_text(size=6))
Fig1_B1
img <- rasterGrob(readPNG(paste0(rac,"CAR2.png")))
Fig1_B2 = ggplot(tt1) +
  geom_blank() +
  annotation_custom(img) +
  labs(title="French West Indies")+
  theme(plot.title=element_text(size=10,face="italic",hjust=0.1))


tt2 = filter(zikachik,REGION=="FP") %>%
  filter(DATE>="2013-10-01",DATE<="2015-04-01") 
tt2$ISLAND2 = gsub(" Isl.","",tt2$ISLAND2)
tt2_h = group_by(tt2,ISLAND2,POP) %>%
  summarise(mm=max(NCASES)*1.2)
tt2_h$mm = c(303,1515,1515,1515,9090,1515)
Fig1_C1 = ggplot() +
  geom_bar(data=tt2,aes(x=DATE,y=NCASES,fill=VIRUS),stat="identity") +
  geom_point(data=tt2_h,aes(x=as.Date("2013-10-01"),y=mm),colour=NA) +
  facet_wrap(~ISLAND2,scales="free_y",ncol=2) +
  scale_y_continuous(expand=c(0,0),#breaks=c(0,50,100)) +
                     labels=function(x) ifelse(x<1000,paste(x),paste0(x/1000,"K"))) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2014-01-01","2015-01-01")),
               labels=c("2014","2015")) +
  scale_fill_manual(values=c(colCHIKV,colZIKV),guide=FALSE) +
  theme_bw2() +
  labs(x="",y="Incident cases (N)") +
  modt + theme(strip.text=element_text(size=8),axis.text.y=element_text(size=6))
Fig1_C1

img2 <- rasterGrob(readPNG(paste0(rac,"FP3.png")))
Fig1_C2 = ggplot(tt1) +
  geom_blank() +
  annotation_custom(img2)+
  labs(title="French Polynesia") +
  theme(plot.title=element_text(size=10,face="italic",hjust=0.1))
Fig1_BC = plot_grid(Fig1_B1,Fig1_B2,Fig1_C1,Fig1_C2,rel_widths=c(1.3,2,2.5,2),nrow=1,labels=c("B","","C",""))

leg = ggplot(zikachik) +
  geom_bar(aes(x=N_WEEK,y=NCASES,fill=VIRUS),stat="identity") +
  scale_fill_manual(values=c(colCHIKV,colZIKV)) +
  labs(fill=NULL) +
  theme(legend.title=element_text(size=9,hjust=0),legend.text=element_text(size=9,hjust=.5))
Fig1_leg = get_legend(leg + theme(legend.position=c(0.4,.8),legend.direction="horizontal"))

Fig1 = plot_grid(Fig1_A,Fig1_BC,Fig1_leg,labels=c("A","",""),rel_heights=c(1,1,.13),ncol=1)
Fig1

ggsave(plot=Fig1,
     filename=paste0(rac,"Fig1.pdf"),
    width=18,height=13,unit="cm")




#########################################################################################
# Fig 2 : priors
# 
# # R0
# bp_R0 = qgamma(c(0.025,0.25,0.5,0.75,0.975),shape=1,rate=0.2)
# 
# sim_mu = rgamma(10000,shape=P_hyperpar_R0$mu$estimate[[1]],rate=P_hyperpar_R0$mu$estimate[[2]])
# sim_sigma = rgamma(10000,shape=P_hyperpar_R0$sigma$estimate[[1]],rate=P_hyperpar_R0$sigma$estimate[[2]])
# sim_tpt = as.vector(mapply(rnorm,n=1000,mean=sim_mu,sd=sim_sigma))
# bp_R0 = rbind(bp_R0,quantile(sim_tpt,probs=c(0.025,0.25,0.5,0.75,0.975)))
# 
# for(island in c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")) {
#   bp_R0 = rbind(bp_R0,qgamma(c(0.025,0.25,0.5,0.75,0.975),shape=P_R0[[island]]$estimate[[1]],rate=P_R0[[island]]$estimate[[2]]))
# }
# 
# bp_R0 = tbl_df(as.data.frame(bp_R0)) %>%
#   mutate(ISLAND=c("","","GUA","MAR","STM"),
#          APPROACH2=c("NI","R","L","L","L"),x=1,par="R0") %>%
#   left_join(colPRIOR)
# 
# 
# # rho
# bp_rho = qbeta(c(0.025,0.25,0.5,0.75,0.975),shape1=1,shape2=1)
# 
# sim_mu = rnorm(10000,mean=P_hyperpar_IL_rho$mu$estimate[[1]],sd=P_hyperpar_IL_rho$mu$estimate[[2]])
# sim_sigma = rgamma(10000,shape=P_hyperpar_IL_rho$sigma$estimate[[1]],rate=P_hyperpar_IL_rho$sigma$estimate[[2]])
# sim_tpt = invlogit(as.vector(mapply(rnorm,n=1000,mean=sim_mu,sd=sim_sigma)))
# bp_rho = rbind(bp_rho,quantile(sim_tpt,probs=c(0.025,0.25,0.5,0.75,0.975)))
# 
# for(island in c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")) {
#   bp_rho = rbind(bp_rho,qgamma(c(0.025,0.25,0.5,0.75,0.975),shape=P_rho[[island]]$estimate[[1]],rate=P_rho[[island]]$estimate[[2]]))
# }
# 
# bp_rho = tbl_df(as.data.frame(bp_rho)) %>%
#   mutate(ISLAND=c("","","GUA","MAR","STM"),
#          APPROACH2=c("NI","R","L","L","L"),x=1,par="rho") %>%
#   left_join(colPRIOR)
# 
# 
# mar = theme(plot.margin=unit(c(.8,0.05,0.3,.3),"cm"))
# expa =   scale_x_discrete(expand=c(.4,0))
# Fig2_A1 = 
#   filter(bp_R0,APPROACH2=="NI") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="#009E73",stat="identity") +
#   theme_bw2() +
#     theme(axis.ticks.x=element_blank()) +
#   mar + expa +
#   scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3,4,5),labels=c(" 0"," 1"," 2"," 3"," 4"," 5")) +
#   coord_cartesian(ylim=c(0,5)) +
#     labs(x=expression(R[0]),y=NULL)
# 
# Fig2_A2 = 
#   filter(bp_rho,APPROACH2=="NI") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="#009E73",stat="identity") +
#   theme_bw2() +
#   theme(axis.ticks.x=element_blank()) +
#   mar + expa +
#   scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,.2,.4,.6,.8,1)) +
#   labs(x=expression(rho),y=NULL)
# 
# Fig2_A = plot_grid(Fig2_A1,Fig2_A2,labels=c("A",""))
# 
# 
# Fig2_B1 = 
#   filter(bp_R0,APPROACH2=="R") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="purple3",stat="identity") +
#   theme_bw2() +
#   theme(axis.ticks.x=element_blank()) + mar + expa +
#   scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3,4,5),labels=c(" 0"," 1"," 2"," 3"," 4"," 5")) +
#   coord_cartesian(ylim=c(0,5)) +
#   labs(x=expression(R[0]),y=NULL)
# 
# Fig2_B2 = 
#   filter(bp_rho,APPROACH2=="R") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="purple3",stat="identity") +
#   theme_bw2() +
#   theme(axis.ticks.x=element_blank()) + mar + expa +
#   scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,.2,.4,.6,.8,1)) +
#   labs(x=expression(rho),y=NULL)
# 
# Fig2_B = plot_grid(Fig2_B1,Fig2_B2,labels=c("B",""))
# 
# 
# Fig2_C1 = 
#   filter(bp_R0,APPROACH2=="L") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="#D55E00",stat="identity") +
#   theme_bw2() +
# mar +  
#   scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3,4,5),labels=c(" 0"," 1"," 2"," 3"," 4"," 5")) +
#   scale_x_discrete(expand=c(0.1,0)) +
#   coord_cartesian(ylim=c(0,5)) +
#   labs(x=expression(R[0]),y=NULL)
# 
# Fig2_C2 = 
#   filter(bp_rho,APPROACH2=="L") %>%
#   ggplot(.,aes(x=ISLAND)) +
#   geom_boxplot(aes(middle=`50%`,lower=`25%`,upper=`75%`,ymin=`2.5%`,ymax=`97.5%`),fill="#D55E00",stat="identity") +
#   theme_bw2() +
# mar +  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=c(0,.2,.4,.6,.8,1)) +
#   scale_x_discrete(expand=c(.1,0)) +
#   labs(x=expression(rho),y=NULL)
# 
# Fig2_C = plot_grid(Fig2_C1,Fig2_C2,labels=c("C",""))
# 
# Fig2 = plot_grid(Fig2_A,Fig2_B,Fig2_C,rel_widths = c(1,1,2),nrow=1)
# Fig2
# ggsave(plot=Fig2,
#        filename=paste0(rac,"Fig2.pdf"),
#        width=18,height=6,unit="cm")
# 
# # 


### priors R_0

sim_mu = rgamma(10000,shape=P_hyperpar_R0$mu$estimate[[1]],rate=P_hyperpar_R0$mu$estimate[[2]])
sim_sigma = rgamma(10000,shape=P_hyperpar_R0$sigma$estimate[[1]],rate=P_hyperpar_R0$sigma$estimate[[2]])
sim_tpt = as.vector(mapply(rnorm,n=1000,mean=sim_mu,sd=sim_sigma))
library(coda)
HPDinterval(as.mcmc(sim_tpt), prob=0.95)
quantile(sim_tpt,probs=c(0.025,0.975))

xrange= seq(0,5,by=0.001)
probR0_NI = data.frame(x=xrange,
                       ISLAND=NA,
                       APPROACH2="NI",
                       prob=dgamma(xrange,shape=1,rate=0.2))
probR0_R = data.frame(x=xrange,
                      ISLAND=NA,
                      APPROACH2="R",
                      prob=density(sim_tpt,from=min(xrange),to=max(xrange),n=length(xrange))$y)
probR0_L = NULL
for(island in c("MARTINIQUE","SAINTMARTIN","GUADELOUPE")) {
  probR0_L = rbind(probR0_L,
                   data.frame(x=xrange,
                              ISLAND=island,
                              APPROACH2="L",
                              prob=dgamma(xrange,shape=P_R0[[island]]$estimate[[1]],rate=P_R0[[island]]$estimate[[2]])))
  library(coda)
  rr=rgamma(10000,shape=P_R0[[island]]$estimate[[1]],rate=P_R0[[island]]$estimate[[2]])
  cat("\n-----\n",island,"\n")
  print(HPDinterval(as.mcmc(rr), prob=0.95))
  cat("mean=");print(mean(rr));cat(" sd=");print(sd(rr))
  cat("\n")
  print(quantile(rr,probs=c(0.025,0.975)))
}

probR0 = rbind(probR0_NI,probR0_R,probR0_L) %>%
  left_join(key) %>%
  left_join(colPRIOR) %>%
  mutate(APPROACH=factor(APPROACH,levels=c("Non-informative\nprior distributions\t","Region-specific\nprior distributions\t","Island-specific\nprior distributions\t")))
invlogit = function(x) exp(x)/(1+exp(x))
sim_mu = rnorm(10000,mean=P_hyperpar_IL_rho$mu$estimate[[1]],sd=P_hyperpar_IL_rho$mu$estimate[[2]])
sim_sigma = rgamma(10000,shape=P_hyperpar_IL_rho$sigma$estimate[[1]],rate=P_hyperpar_IL_rho$sigma$estimate[[2]])
sim_tpt = invlogit(as.vector(mapply(rnorm,n=1000,mean=sim_mu,sd=sim_sigma)))
HPDinterval(as.mcmc(sim_tpt), prob=0.95)
quantile(sim_tpt,probs=c(0.025,0.975))

xrange= seq(0,1,by=0.001)
mar = 35
probrho_NI = data.frame(x=xrange,
                        ISLAND=NA,
                        APPROACH2="NI",
                        prob=dbeta(xrange,shape1=1,shape2=1))
probrho_R = data.frame(x=xrange,
                       ISLAND=NA,
                       APPROACH2="R",
                       prob=density(sim_tpt,from=min(xrange),to=max(xrange),n=length(xrange))$y)
probrho_L = NULL
for(island in c("MARTINIQUE","SAINTMARTIN","GUADELOUPE")) {
  probrho_L = rbind(probrho_L,
                    data.frame(x=xrange,
                               ISLAND=island,
                               APPROACH2="L",
                               prob=dgamma(xrange,shape=P_rho[[island]]$estimate[[1]],rate=P_rho[[island]]$estimate[[2]])))
  library(coda)
  rr=rgamma(10000,shape=P_rho[[island]]$estimate[[1]],rate=P_rho[[island]]$estimate[[2]])
  cat("\n-----\n",island,"\n")
  print(HPDinterval(as.mcmc(rr), prob=0.95))
  cat("mean=");print(mean(rr));cat(" sd=");print(sd(rr))
  cat("\n")
  print(quantile(rr,probs=c(0.025,0.975)))
  cat("\nmean")
  print(mean(rr))
}

probrho = rbind(probrho_NI,probrho_R,probrho_L) %>%
  left_join(key) %>%
  left_join(colPRIOR) %>%
  mutate(APPROACH=factor(APPROACH,levels=c("Non-informative\nprior distributions\t","Region-specific\nprior distributions\t","Island-specific\nprior distributions\t")))


Fig2C_NI = ggplot(filter(probR0,APPROACH2=="NI")) +
  geom_ribbon(aes(x=x,ymax=prob,ymin=0),fill="#009E73") +
  geom_line(aes(x=x,y=prob)) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,3.5),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5.5),expand=c(0,0),labels=function(x) paste0(" ",x)) +
  annotate("text",x=0.2,y=5.1,label="Non-informative",hjust=0,size=3.5,fontface="italic") +
  labs(y="PDF",x=expression(R[0][Z])) 

Fig2C_R = ggplot(filter(probR0,APPROACH2=="R")) +
  geom_ribbon(aes(x=x,ymax=prob,ymin=0),fill="purple3") +
  geom_line(aes(x=x,y=prob)) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,3.5),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5.5),expand=c(0,0),labels=function(x) paste0(" ",x)) +
  annotate("text",x=0.2,y=5.1,label="Regional",hjust=0,size=3.5,fontface="italic") +
  labs(y=NULL,x=expression(R[0][Z])) 

tmp1 = filter(probR0,APPROACH2=="L") %>%
  mutate(shift=ifelse(ISLAND2=="Guadeloupe",-2,
                      ifelse(ISLAND2=="Martinique",0,2)))
tmp2 = filter(tmp1,ISLAND2=="Guadeloupe")
Fig2C_L = 
  ggplot() +
  geom_ribbon(data=tmp1,aes(x=x,ymax=prob,ymin=0,group=ISLAND2),fill="#D55E00") +
  geom_line(data=tmp1,aes(x=x,y=prob,group=ISLAND2)) +
  geom_ribbon(data=tmp2,aes(x=x,ymax=prob,ymin=0,group=ISLAND2),fill="#D55E00") +
  geom_line(data=tmp2,aes(x=x,y=prob,group=ISLAND2)) +
  geom_line(data=tmp1,aes(x=x,y=prob,group=ISLAND2),linetype=2,size=.3) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,3.5),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5.5),expand=c(0,0),labels=function(x) paste0(" ",x)) +
  geom_segment(aes(x=1.75,xend=1.3,y=4.75,yend=3.85)) +
  annotate("text",x=1.80,y=4.8,label="Martinique",hjust=0,size=3) +
  geom_segment(aes(x=1.9,xend=1.55,y=3.75,yend=3.35)) +
  annotate("text",x=1.95,y=3.8,label="Guadeloupe",hjust=0,size=3) +
  geom_segment(aes(x=1.35,xend=1.95,y=2.07,yend=2.75)) +
  annotate("text",x=2,y=2.8,label="St-Martin",hjust=0,size=3) +
  annotate("text",x=0.2,y=5.1,label="Local",hjust=0,size=3.5,fontface="italic") +
  labs(x=expression(R[0][Z]),y=NULL)


### priors for rho




Fig2D_NI = ggplot(filter(probrho,APPROACH2=="NI")) +
  geom_ribbon(aes(x=x,ymax=prob,ymin=0),fill="#009E73") +
  geom_line(aes(x=x,y=prob)) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,.25,.5,.75,1),labels=c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous(limits=c(0,mar),expand=c(0,0))+
  annotate("text",x=0.06,y=mar*.92,label="Non-informative",hjust=0,size=3.5,fontface="italic") +
  labs(y="PDF",x=expression (rho[Z])) 

Fig2D_R = ggplot(filter(probrho,APPROACH2=="R")) +
  geom_ribbon(aes(x=x,ymax=prob,ymin=0),fill="purple3") +
  geom_line(aes(x=x,y=prob)) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,.25,.5,.75,1),labels=c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous(limits=c(0,mar),expand=c(0,0))+
  annotate("text",x=0.06,y=mar*.92,label="Regional",hjust=0,size=3.5,fontface="italic") +
  labs(y=NULL,x=expression(rho[Z])) 

tmp1 = filter(probrho,APPROACH2=="L") %>%
  mutate(shift=ifelse(ISLAND2=="Guadeloupe",-10,
                      ifelse(ISLAND2=="Martinique",0,10)))
tmp2 = filter(tmp1,ISLAND2=="Guadeloupe")
tmp3 = data.frame(shift=c(-10,0,10),lab=c("Guadeloupe","Martinique","Saint-Martin"))
Fig2D_L = 
  ggplot() +
  geom_ribbon(data=tmp1,aes(x=x,ymax=prob,ymin=0,group=ISLAND2),fill="#D55E00") +
  geom_line(data=tmp1,aes(x=x,y=prob,group=ISLAND2)) +
  geom_ribbon(data=tmp2,aes(x=x,ymax=prob,ymin=0,group=ISLAND2),fill="#D55E00") +
  geom_line(data=tmp2,aes(x=x,y=prob,group=ISLAND2)) +
  geom_line(data=tmp1,aes(x=x,y=prob,group=ISLAND2),linetype=2,size=.3) +
  theme_bw2() +
  scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,.25,.5,.75,1),labels=c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous(limits=c(0,mar),expand=c(0,0)) +
  geom_segment(aes(x=0.2,xend=0.35,y=24,yend=24.8)) +
  annotate("text",x=0.36,y=25,label="Guadeloupe",hjust=0,size=3) +
  geom_segment(aes(x=0.25,xend=0.40,y=14,yend=16.8)) +
  geom_segment(aes(x=0.14,xend=0.25,y=4.6,yend=14)) +
  annotate("text",x=0.42,y=17,label="St-Martin",hjust=0,size=3) +
  geom_segment(aes(x=0.32,xend=0.5,y=7,yend=10)) +
  annotate("text",x=0.51,y=10.2,label="Martinique",hjust=0,size=3) +
  annotate("text",x=0.06,y=mar*.92,label="Local",hjust=0,size=3.5,fontface="italic") +
  labs(x=expression(rho[Z]),y=NULL)

Fig2C = plot_grid(Fig2C_NI+modt,
                  Fig2C_R+modt,
                  Fig2C_L+modt,nrow=1)
Fig2D = plot_grid(Fig2D_NI + modt,
                  Fig2D_R + modt,
                  Fig2D_L + modt,nrow=1, rel_heights=c(.88,.88,1))

# Join priors

FigPRIOR <- plot_grid(Fig2D,Fig2C,nrow=2,labels=c("A","B"))
FigPRIOR
ggsave(plot=FigPRIOR,
       filename=paste0(rac,"Fig2.pdf"),
       width=18,height=9,unit="cm")
# convert -density 1200 FigPRIOR.pdf -background white -alpha remove FigPRIOR2.pdf




#########################################################################################
# Fig 3 : Prediction


plot_point_example <- function(prior_r0,prior_rho,datechoice,...) {
  ex_fit <- filter(TSUM_E_NCASES,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
    left_join(datechoice) %>%
    filter(DATECUT==choice,DATE<="2017-02-19")
  ex_pred <- filter(TSUM_F_NCASES,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
    left_join(datechoice) %>%
    filter(DATECUT==choice,DATE<="2017-02-19") %>%
    mutate(`97.5%`=ifelse(`97.5%`>topline,topline,`97.5%`),
           `75%`=ifelse(`75%`>topline,topline,`75%`),
           mean=ifelse(mean>topline,topline,mean)) %>%
    left_join(colPRIOR)
  ex_real <- left_join(fullD1,datechoice,by="ISLAND") %>%
    filter(DATE>choice,DATE<="2017-02-19")
  ex_acc = filter(TSUM_RMSD,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
    left_join(datechoice) %>%
    filter(DATECUT==choice) %>%
    mutate(ACC=paste0("ACC=",round(mean)))
  ex_sharp = filter(TSUM_F_NCASES,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
    left_join(datechoice) %>%
    filter(DATECUT==choice,DATE<="2017-02-19") %>%
    mutate(SHARP50=(`75%`-`25%`),SHARP95=(`97.5%`-`2.5%`)) %>%
    group_by(ISLAND,DATECUT,PRIOR_R0,PRIOR_rho) %>%
    summarise(SHARP50=mean(SHARP50),SHARP95=mean(SHARP95)) %>%
    left_join(datechoice) %>%
    mutate(SHA=paste0("SHA=",round(SHARP95)))
  g <- ggplot() +
    geom_ribbon(data=ex_pred,aes(x=DATE,ymin=`2.5%`,ymax=`97.5%`),alpha=.1) +
    geom_ribbon(data=ex_pred,aes(x=DATE,ymin=`25%`,ymax=`75%`),alpha=.2) +
    geom_point(data=ex_fit,aes(x=DATE,y=NCASES),shape=21,fill="black",size=1) +
    geom_point(data=ex_real,aes(x=DATE,y=NCASES),shape=21,size=.9) +
    geom_line(data=ex_pred,aes(x=DATE,y=mean),colour=unique(as.character(ex_pred$col)),size=0.8) +
    geom_vline(data=datechoice,aes(xintercept=as.numeric(choice)+7),linetype=2) +
    geom_point(data=datechoice,aes(x=as.Date("2015-12-20"),y=topline),col="white") +
    facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
    scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
                 labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
    scale_colour_discrete(guide=FALSE) +
    scale_y_continuous(expand=c(0,0),label=function(x) paste0(x/1000,"K")) +
    theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA),
                        strip.background=element_blank(),
                        axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
    modt +
    labs(x="Date",y="Incident cases (N)")
  return(g)
}

dc1 = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
                 choice=14+as.Date(c("2016-04-10","2016-01-10","2016-02-21")),
                 topline=c(12500,12500,1250)) %>%
  left_join(timing,by="ISLAND")
Fig3_A = plot_point_example(prior_r0=0,prior_rho=0,datechoice=dc1)
Fig3_B = plot_point_example(prior_r0=2,prior_rho=2,datechoice=dc1)
Fig3_C = plot_point_example(prior_r0=1,prior_rho=1,datechoice=dc1)

Fig3 = plot_grid(Fig3_A,Fig3_B,Fig3_C,ncol=1,labels=c("A","B","C"))
Fig3 = plot_grid(Fig3,Fig_legH,ncol=1,rel_heights = c(1,.07))
Fig3
ggsave(plot=Fig3,
       filename=paste0(rac,"Fig3.pdf"),
       width=18,height=16,unit="cm")


## Animation
for(i in 5:nrow(dates)) {
  dc1 = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
                   choice=as.Date(rep(dates[i,],3)),
                   topline=c(12500,12500,1250)) %>%
    left_join(timing,by="ISLAND")
  tmp_A = plot_point_example(prior_r0=0,prior_rho=0,datechoice=dc1)
  tmp_B = plot_point_example(prior_r0=2,prior_rho=2,datechoice=dc1)
  tmp_C = plot_point_example(prior_r0=1,prior_rho=1,datechoice=dc1)
  
  tmp = plot_grid(tmp_A,tmp_B,tmp_C,ncol=1)
  tmp = plot_grid(tmp,Fig_legH,ncol=1,rel_heights = c(1,.07),labels=i)
  ggsave(filename=paste0("animation/frame_",ifelse(i<10,"0",""),i,".png"),
         plot=tmp,device="png",width=16,height=14,unit="cm")
  cat(".")
}

# convert 




#########################################################################################
# Fig 4: accuracy and sharpness 

# accuracy

tmpacc = filter(TSUM_RMSD,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(colPRIOR) %>%
  # left_join(key) %>%
  left_join(timing) %>%
  filter(DATECUT<=end)
maxscale1 = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
                 topline=c(10,10,1)) %>%
  left_join(timing,by="ISLAND") %>%
  mutate(min=as.Date("2015-12-20",max=as.Date("2017-02-19")))
Fig4A = ggplot() +
  geom_line(data=tmpacc,aes(x=DATECUT,y=mean/1000,colour=APPROACH),size=1.1) +
  
  geom_point(data=maxscale1,aes(x=min,y=topline*1.15),colour="white") +
  geom_point(data=maxscale1,aes(x=max,y=topline*1.15),colour="white") +
  
  geom_segment(data=maxscale1,aes(x=start-1,xend=start-1,y=0,yend=topline),alpha=.4) +
  geom_segment(data=maxscale1,aes(x=peak-1,xend=peak-1,y=0,yend=topline),alpha=.4) +
  geom_segment(data=maxscale1,aes(x=end-1,xend=end-1,y=0,yend=topline),alpha=.4) +
  
  geom_label(data=maxscale1,aes(x=start-1,y=topline),label="S",size=2.5) +
  geom_label(data=maxscale1,aes(x=peak-1,y=topline),label="P",size=2.5) +
  geom_label(data=maxscale1,aes(x=end-1,y=topline),label="E",size=2.5) +
  
  facet_wrap(~ ISLAND2,scales="free_y") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0),label=unit_format("K")) +
  #scale_y_continuous(expand=c(0,0),limits=c(10,10^(log10(maxscale1)*1.1)),breaks=c(10,100,1000,10000,100000),
                # labels=c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5))) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
               labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
  labs(x="Date of assessment",y="Accuracy") 
Fig4A

# sharpness

tmpsharp = filter(TSUM_F_NCASES,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  filter(!is.na(NCASES)) %>%
  mutate(SHARP50=(`75%`-`25%`)/1000,SHARP95=(`97.5%`-`2.5%`)/1000) %>%
  group_by(ISLAND,DATECUT,PRIOR_R0,PRIOR_rho) %>%
  summarise(SHARP50=mean(SHARP50),SHARP95=mean(SHARP95)) %>%
  left_join(colPRIOR) %>%
  left_join(timing) %>%
  filter(DATECUT<=end)
maxscale2 = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
                       topline=c(32,32,3.2)) %>%
  left_join(timing,by="ISLAND") %>%
  mutate(min=as.Date("2015-12-20",max=as.Date("2017-02-19")))
maxscale2$scale_sharp=c(36,36,3.6)
Fig4B = ggplot() +
  geom_line(data=tmpsharp,aes(x=DATECUT,y=SHARP95,colour=APPROACH),size=1.1) +
  geom_point(data=maxscale2,aes(x=min,y=scale_sharp*1.15),colour="white") +
  geom_point(data=maxscale2,aes(x=max),y=8.9,colour="white") +
  
  geom_segment(data=maxscale2,aes(x=start-1,xend=start-1,y=0,yend=scale_sharp),alpha=.4) +
  geom_segment(data=maxscale2,aes(x=peak-1,xend=peak-1,y=0,yend=scale_sharp),alpha=.4) +
  geom_segment(data=maxscale2,aes(x=end-1,xend=end-1,y=0,yend=scale_sharp),alpha=.4) +
  
  geom_label(data=maxscale2,aes(x=start-1,y=scale_sharp),label="S",size=2.5) +
  geom_label(data=maxscale2,aes(x=peak-1,y=scale_sharp),label="P",size=2.5) +
  geom_label(data=maxscale2,aes(x=end-1,y=scale_sharp),label="E",size=2.5) +
  
  facet_wrap(~ ISLAND2,scales="free_y") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0),label=unit_format("K")) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
               labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
  labs(x="Date of assessment",y="Sharpness")
Fig4B


gA <- ggplotGrob(Fig4A+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.2, "cm")) + modt)
gB <- ggplotGrob(Fig4B+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.3, "cm")) + modt)
maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth

Fig4 = plot_grid(gA,gB,ncol=1,labels=LETTERS[1:6])
Fig4 = plot_grid(Fig4,Fig_legH,ncol=1,rel_heights = c(1,.1))
Fig4
ggsave(plot=Fig4,
       filename=paste0(rac,"Fig4.pdf"),
       width=18,height=12,unit="cm")




#########################################################################################
# Fig 5: posterior distributions

tmpR0 = filter(TSUM_R0,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(colPRIOR) %>%
  left_join(timing) %>%
  # left_join(key) %>%
  mutate(min=`2.5%`,max=`97.5%`)
Fig5A = ggplot() +
  geom_line(data=tmpR0,aes(x=DATECUT,y=mean,colour=APPROACH),size=1.1) +
  geom_ribbon(data=tmpR0,aes(x=DATECUT,ymin=min,ymax=max,fill=APPROACH),alpha=.1) +
  geom_line(data=tmpR0,aes(x=DATECUT,y=min,colour=APPROACH),size=.2,linetype=3) +
  geom_line(data=tmpR0,aes(x=DATECUT,y=max,colour=APPROACH),size=.2,linetype=3) +
  geom_point(data=timing,aes(x=min,y=9.2),colour="white") +
  geom_point(data=timing,aes(x=max,y=.5),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=0,yend=8.5),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=0,yend=8.5),alpha=.4) +
  geom_segment(data=timing,aes(x=end-1,xend=end-1,y=0,yend=8.5),alpha=.4) +
  
  geom_label(data=timing,aes(x=start-1,y=8.5),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=8.5),label="P",size=2.5) +
  geom_label(data=timing,aes(x=end-1,y=8.5),label="E",size=2.5) +
  
  facet_wrap(~ ISLAND2) +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
               labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
  labs(x="Date of assessment",y=expression(R[0][Z]))
Fig5A

tmprho = filter(TSUM_rho,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(colPRIOR) %>%
  left_join(timing) %>%
  mutate(min=`2.5%`,max=`97.5%`)
Fig5B = ggplot() +
  geom_line(data=tmprho,aes(x=DATECUT,y=mean,colour=APPROACH),size=1.1) +
  geom_ribbon(data=tmprho,aes(x=DATECUT,ymin=min,ymax=max,fill=APPROACH),alpha=.1) +
  geom_line(data=tmprho,aes(x=DATECUT,y=min,colour=APPROACH),size=.2,linetype=3) +
  geom_line(data=tmprho,aes(x=DATECUT,y=max,colour=APPROACH),size=.2,linetype=3) +
  geom_point(data=timing,aes(x=min,y=1.2),colour="white") +
  geom_point(data=timing,aes(x=max,y=1.2),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=0,yend=1.1),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=0,yend=1.1),alpha=.4) +
  geom_segment(data=timing,aes(x=end-1,xend=end-1,y=0,yend=1.1),alpha=.4) +
  
  geom_label(data=timing,aes(x=start-1,y=1.1),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=1.1),label="P",size=2.5) +
  geom_label(data=timing,aes(x=end-1,y=1.1),label="E",size=2.5) +
  
  facet_wrap(~ ISLAND2) +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0),limits=c(0,1.2),label=scales::percent,
                     breaks=c(0,0.25,0.50,0.75,1)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
               labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) +
  labs(x="Date of assessment",y=expression(rho[Z]))
Fig5B

gA <- ggplotGrob(Fig5A+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.2, "cm")) + modt)
gB <- ggplotGrob(Fig5B+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.3, "cm"))+ modt)

maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])

gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth

Fig5 = plot_grid(gA,gB,ncol=1,labels=c("A","B"))
Fig5 = plot_grid(Fig5,Fig_legH,ncol=1,rel_heights = c(1,.07))
Fig5
ggsave(plot=Fig5,
       filename=paste0(rac,"Fig5.pdf"),
       width=18,height=12,unit="cm")





#########################################################################################
# Fig 6: indicators" 



# total cases

tmptot = filter(TSUM_F_OBSERVEDCASES,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  mutate(ISLAND=as.character(ISLAND)) %>%
  left_join(colPRIOR) %>%
  # left_join(key) %>%
  left_join(timing) %>%
  filter(DATECUT<=peak)
timing$tmptot = c(205,205,20.5)
tmp.totobs2 = group_by(fullD1,ISLAND) %>%
  summarise(TT=sum(NCASES)/1000) 
timing$tmp.totobs2 = tmp.totobs2$TT
Fig6E = ggplot() +
  geom_line(data=tmptot,aes(x=DATECUT,y=mean/1000,colour=APPROACH),size=1.1) +
  
  geom_hline(data=timing,aes(yintercept = tmp.totobs2),linetype=2) +
  
  geom_point(data=timing,aes(x=as.Date("2016-02-01"),y=tmptot*1.15),colour="white") +
  geom_point(data=timing,aes(x=start-7,y=0),colour="white") +
  geom_point(data=timing,aes(x=peak+(peak-start)*.1,y=0),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=0,yend=tmptot),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=0,yend=tmptot),alpha=.4) +
  
  geom_label(data=timing,aes(x=start-1,y=tmptot),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=tmptot),label="P",size=2.5) +
  
  facet_wrap(~ ISLAND2,scales="free") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(label=function(x) paste0(x,"K"),expand=c(0,0)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-02-01","2016-03-01","2016-04-01","2016-05-01","2016-06-01","2016-07-01","2016-08-01")),
               labels=c("Jan.\n2016","Feb.\n2016","March\n","April\n","May\n","June\n","July\n","Aug.\n")) +
  labs(x="Date of assessment",y=(("Total epidemic size\n(N)")))
Fig6E
# 
# filter(tmptot,DATECUT<=peak) %>%
#   left_join(timing) %>%
#   select(tmptot,ISLAND,DATECUT,APPROACH,mean,tmp.totobs2) %>%
#   mutate(mean2=mean/1000/tmp.totobs2) %>%
#   group_by(APPROACH) %>%
#   summarise(mean=mean(mean2),min=min(mean2),max=max(mean2))
# 
# filter(tmptot,DATECUT<=peak) %>%
#   filter(ISLAND=="MARTINIQUE") %>%
#   group_by(APPROACH) %>%
#   View()
# 
# filter(tmptot,DATECUT<=peak) %>%
#   filter(ISLAND=="SAINTMARTIN",PRIOR_R0==0) %>%
#   group_by(APPROACH) %>%
#   View()


# maximal incidence

timing = D1 %>%
  group_by(ISLAND) %>%
  summarise(MAXINC=max(NCASES)) %>%
  right_join(timing) 

timing$maxincmarg = c(108,108,10.8)
timing = left_join(timing,key)
tmpmaxinc = filter(TSUM_F_MAXINC,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  mutate(ISLAND=as.factor(ISLAND)) %>%
  left_join(colPRIOR) %>%
  left_join(timing) %>%
  filter(DATECUT<=peak) %>%
  mutate(mean=mean/1000,
         `97.5%`=`97.5%`/1000,
         `2.5%`=`2.5%`/1000)
Fig6C = ggplot() +
  geom_line(data=tmpmaxinc,aes(x=DATECUT,y=`97.5%`,colour=APPROACH),size=1.1) +
  
  geom_hline(data=timing,aes(yintercept = MAXINC/1000),linetype=2) +
  
  geom_point(data=timing,aes(x=as.Date("2016-02-01"),y=maxincmarg*1.15),colour="white") +
  geom_point(data=timing,aes(x=start-7,y=0),colour="white") +
  geom_point(data=timing,aes(x=peak+(peak-start)*.1,y=0),colour="white") +

  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=0,yend=maxincmarg),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=0,yend=maxincmarg),alpha=.4) +

  geom_label(data=timing,aes(x=start-1,y=maxincmarg),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=maxincmarg),label="P",size=2.5) +

  facet_wrap(~ ISLAND2,scales="free") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(label=function(x) paste0(x,"K"),expand=c(0,0)) +
  # scale_y_log10(labels=trans_format('log10',math_format(10^.x)),breaks=c(100,1000,10000,100000,1000000)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-02-01","2016-03-01","2016-04-01","2016-05-01","2016-06-01","2016-07-01","2016-08-01")),
               labels=c("Jan.\n2016","Feb.\n2016","March\n","April\n","May\n","June\n","July\n","Aug.\n")) +
  labs(x="Date of assessment",y=(("Maximal weekly\nincidence (N)")))
Fig6C

tmpmaxinc %>%
  filter(DATECUT<=peak) %>%
  select(DATECUT,ISLAND,APPROACH2,MAXINC,`97.5%`) %>%
  mutate(mean2=`97.5%`/MAXINC*1000) %>%
  group_by(APPROACH2) %>%
  summarise(mean=mean(mean2),min=min(mean2),max=max(mean2))

filter(TSUM_F_MAXINC,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(colPRIOR) %>%
  left_join(timing) %>%
  left_join(key) %>%
  mutate(mean=mean/1000,
         `97.5%`=`97.5%`/1000,
         `2.5%`=`2.5%`/1000) %>%
  filter(DATECUT>peak) %>%
  select(DATECUT,ISLAND,APPROACH2,MAXINC,`97.5%`) %>%
  mutate(mean2=`97.5%`/MAXINC*1000) %>%
  group_by(APPROACH2) %>%
  summarise(min=min(mean2),mean=mean(mean2),max=max(mean2))


# date of peak incidence

timing = left_join(fullD1,timing) %>%
  mutate(peak.nweek=ifelse(DATE==peak,1,0)) %>% 
  filter(peak.nweek==1) %>%
  select(ISLAND,peak.week=N_WEEK) %>%
  right_join(timing)

tmppeak = filter(TSUM_F_PEAK,PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(colPRIOR) %>%
  # mutate(ISLAND=as.character(ISLAND)) %>%
  left_join(timing) %>%
  filter(DATECUT<=peak) %>%
  mutate(PEAK.GAP=mean-peak.week,PEAK.GAP.MIN=`2.5%`-peak.week,PEAK.GAP.MAX=`97.5%`-peak.week)
maxscale2 = 5
Fig6D = ggplot() +
  geom_line(data=tmppeak,aes(x=DATECUT,y=PEAK.GAP*7*12/365.25,colour=APPROACH),size=1.1) +

  geom_hline(data=timing,aes(yintercept = 0),linetype=2) +
  
  geom_point(data=timing,aes(x=as.Date("2016-02-01"),y=maxscale2*1.15),colour="white") +
  geom_point(data=timing,aes(x=start-7,y=0),colour="white") +
  geom_point(data=timing,aes(x=peak+(peak-start)*.1,y=0),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=-3.8,yend=maxscale2),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=-3.8,yend=maxscale2),alpha=.4) +

  geom_label(data=timing,aes(x=start-1,y=maxscale2),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=maxscale2),label="P",size=2.5) +

  facet_wrap(~ ISLAND2,scales="free") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0),limits=c(-3.9,(maxscale2)*1.3),breaks=c(-2,0,2,4,6)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-02-01","2016-03-01","2016-04-01","2016-05-01","2016-06-01","2016-07-01","2016-08-01")),
               labels=c("Jan.\n2016","Feb.\n2016","March\n","April\n","May\n","June\n","July\n","Aug.\n")) +
  labs(x="Date of assessment",y="Peak date difference\n(months)")
Fig6D

tmppeak %>%
  filter(DATECUT<=peak) %>%
  group_by(APPROACH2) %>%
  summarise(mean=mean(PEAK.GAP*7*12/365.25),min=min(PEAK.GAP*7*12/365.25),max=max(PEAK.GAP*7*12/365.25))

filter(tmppeak,DATECUT<peak,DATECUT>=start,!(ISLAND=="SAINTMARTIN" & DATECUT<"2017-06-01")) %>%
  left_join(timing) %>%
  group_by(APPROACH2) %>%
  summarise(mean=mean(PEAK.GAP),min=min(PEAK.GAP),max=max(PEAK.GAP))



# date of end

tmpend = filter(bind_rows(TSUM_E_NCASES,TSUM_F_NCASES),PRIOR_R0==0 & PRIOR_rho==0 | PRIOR_R0==1 & PRIOR_rho==1 | PRIOR_R0==2 & PRIOR_rho==2) %>%
  left_join(timing) %>%
  filter(DATECUT<=peak) %>%
  mutate(THRESHOLD=ifelse(ISLAND=="SAINTMARTIN",20,200)) %>%
  arrange(DATECUT,ISLAND,PRIOR_R0,PRIOR_rho,DATE) %>%
  group_by(DATECUT,ISLAND,PRIOR_R0,PRIOR_rho) %>%
  mutate(ISTHRESHOLD=mean<THRESHOLD,
         ISTHRESHOLD1=lag(ISTHRESHOLD,1),
         ISTHRESHOLD2=lag(ISTHRESHOLD,2),
         END=ifelse(ISTHRESHOLD+ISTHRESHOLD1+ISTHRESHOLD2==3,1,0),
         END1=lag(END,1),
         END=ifelse(END==1 & END1==0,1,0)
         ) %>% 
  filter(END==1) %>%
  mutate(MAXDATE=ifelse(DATE==max(DATE),1,0)) %>%
  filter(MAXDATE==1) %>%
  left_join(colPRIOR) %>%
  mutate(ENDDATE=(DATE-end)/7,
         PHEA=(DATE-start)/7)
maxscale3 = 40
# Fig6F = ggplot() +
#   geom_line(data=tmpend,aes(x=DATECUT,y=ENDDATE,colour=APPROACH),size=1.1) +
#   
#   geom_hline(data=timing,aes(yintercept = 0),linetype=2) +
#   
#   geom_point(data=timing,aes(x=min,y=maxscale3),colour="white") +
#   geom_point(data=timing,aes(x=peak+(peak-start)*.1,y=0),colour="white") +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=-40,yend=maxscale3),alpha=.4) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=-40,yend=maxscale3),alpha=.4) +
#   
#   geom_label(data=timing,aes(x=start-1,y=maxscale3),label="S",size=2.5) +
#   geom_label(data=timing,aes(x=peak-1,y=maxscale3),label="P",size=2.5) +
#   
#   
#   facet_wrap(~ ISLAND2,scales="free") +
#   theme_bw2() +
#   scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
#   scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
#   scale_y_continuous(expand=c(0,0),limits=c(-43,maxscale3*1.2),breaks=c(-40,-20,0,20,40),
#                      labels=c("-40","-20","0","+20","+40")) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-02-01","2016-03-01","2016-04-01","2016-05-01","2016-06-01","2016-07-01","2016-08-01")),
#                labels=c("Jan.\n2016","Feb.\n","March\n","April\n","May\n","June\n","July\n","Aug.\n")) +
#   labs(x="Date",y=expression(paste("End date (",Delta,"week)")))
# Fig6F
timing$dur = as.numeric((timing$end - timing$start)/7)
maxscale3 = 20
Fig6F = ggplot() +
  geom_line(data=tmpend,aes(x=DATECUT,y=as.numeric(PHEA)*7/365.25*12,colour=APPROACH),size=1.1) +

  geom_hline(data=timing,aes(yintercept = dur*7/365.25*12),linetype=2) +
  
  geom_point(data=timing,aes(x=as.Date("2016-02-01"),y=maxscale3),colour="white") +
  geom_point(data=timing,aes(x=start-7,y=0),colour="white") +
  geom_point(data=timing,aes(x=peak+(peak-start)*.1,y=0),colour="white") +
  
  geom_segment(data=timing,aes(x=start-1,xend=start-1,y=0,yend=maxscale3),alpha=.4) +
  geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=0,yend=maxscale3),alpha=.4) +

  geom_label(data=timing,aes(x=start-1,y=maxscale3),label="S",size=2.5) +
  geom_label(data=timing,aes(x=peak-1,y=maxscale3),label="P",size=2.5) +

  
  facet_wrap(~ ISLAND2,scales="free") +
  theme_bw2() +
  scale_colour_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_fill_manual(values=as.character(colPRIOR$col),guide=F) +
  scale_y_continuous(expand=c(0,0),limits=c(0,maxscale3*1.15),
                     breaks=c(0,6,12,18)) +
  scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-02-01","2016-03-01","2016-04-01","2016-05-01","2016-06-01","2016-07-01","2016-08-01")),
               labels=c("Jan.\n2016","Feb.\n2016","March\n","April\n","May\n","June\n","July\n","Aug.\n")) +
  labs(x="Date of assessment",y=(("Epidemic duration\n(months)")))
Fig6F


filter(tmpend,DATECUT<=peak) %>%
  group_by(APPROACH2) %>%
  left_join(timing) %>%
  mutate(dur = dur*7/365.25*12,
         PHEA=PHEA*7/365.25*12) %>%
  summarise(mean=mean(PHEA/dur),min=min(PHEA/dur),max=max(PHEA/dur))




gE <- ggplotGrob(Fig6E+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(.2, "cm")) + modt)
gC <- ggplotGrob(Fig6C+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.2, "cm")) + modt)
gD <- ggplotGrob(Fig6D+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(0.5, "cm")) + modt)
gF <- ggplotGrob(Fig6F+ theme(plot.margin=unit(rep(.2,4),"cm"),panel.spacing = unit(.4, "cm")) + modt)
maxWidth = unit.pmax(gC$widths[2:3], gD$widths[2:3],
                     gE$widths[2:3], gF$widths[2:3])

# Set the widths
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth


Fig6 = plot_grid(gE,gC,gD,gF,ncol=1,labels=LETTERS[1:6])
Fig6 = plot_grid(Fig6,Fig_legH,ncol=1,rel_heights = c(1,.07))
Fig6

ggsave(plot=Fig6,
       filename=paste0(rac,"Fig6.pdf"),
       width=18,height=22,unit="cm")

 # 
# 
# 
# 
# 
# 
# 
#  # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Piecewise regression with 2 knots
# 
# pw = filter(zikachik,ISLAND=="MARTINIQUE",VIRUS=="ZIKV")
# plot(NCASES ~ N_WEEK,data=pw,type="l")
# lin.mod = lm(NCASES ~ N_WEEK,data=pw)
# segm.mod = segmented(lin.mod,seg.Z=~N_WEEK)
# summary(segm.mod)
# points(segm.mod,col="red")
# points(NCASES ~ N_WEEK,data=pw,type="l")
# points(lowess(x=pw$N_WEEK,y=pw$NCASES,f=.3),type="l",col="blue")
# 
# getbreakingpoints <- function(nc,threshold=200) {
#   require(segmented)
#   dt = data.frame(nc=nc,nw=1:length(nc))
#   lin = lm(nc ~ nw,data=dt)
#   seg = segmented(lin,seg.Z=~nw)
#   peak = seg$psi[,2:3]
#   end = dt$nw[which(dt$nc<threshold)]
#   end = min(end[end>peak[1]])
#   names(end) = "Est."
#   return(list(peak=peak,end=end))
# }
# getbreakingpoints(pw$NCASES,pw$N_WEEK)
# 
# 
# 
# 
# 
# ##################################################################################################################
# # Figure with evolution of accuracy, sharpness, total observed
# 
# # Total observed -------------------------------------------------------------------------
# 
# tmp.totobs2 = group_by(fullD1,ISLAND) %>%
#   summarise(TT=sum(NCASES)/1000)
# tmp.totobs= mutate(TSUM_TOTALOBSCASES,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                    METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                    METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                    METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                    METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                    ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# 
# marg.totobs = .07
# Fig_totobs = 
#   ggplot() +
#   geom_line(data=tmp.totobs,aes(x=DATECUT,y=mean/1000,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=top*marg.totobs*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=top*marg.totobs*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=top*marg.totobs*.9,yend=0),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=top*marg.totobs*.9),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=top*marg.totobs*.9),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=top*marg.totobs*.9),label="E",size=3.5) +
#   
#   geom_point(data=timing,aes(x=min,y=top*marg.totobs),colour="white") +
#   geom_point(data=timing,aes(x=max,y=0),colour="white") +
#   
#   geom_hline(data=tmp.totobs2,aes(yintercept=TT),linetype=2) +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="N (x1,000)",colour="Approach") +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# # Accuracy -------------------------------------------------------------------------
# 
# tmp.rmsd = mutate(TSUM_RMSD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                   METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                   METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                   ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# 
# marg.rmsd = 3.3
# Fig_rmsd = 
#   ggplot() +
#   geom_line(data=tmp.rmsd,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=top*marg.rmsd*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=top*marg.rmsd*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=top*marg.rmsd*.9,yend=0),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=top*marg.rmsd*.9),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=top*marg.rmsd*.9),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=top*marg.rmsd*.9),label="E",size=3.5) +
#   
#   geom_point(data=timing,aes(x=min,y=top*marg.rmsd),colour="white") +
#   geom_point(data=timing,aes(x=max,y=0),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="RMSD",colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# # Sharpness -------------------------------------------------------------------------
# 
# tmp.sharp = mutate(TSUM_PRED,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                    METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                    METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                    METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                    METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                    ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
#   left_join(fullD1) %>%
#   filter(!is.na(NCASES)) %>%
#   select(DATECUT,ISLAND,METHOD,PRIOR_R0,PRIOR_rho,`2.5%`,`25%`,`75%`,`97.5%`,DATE,NCASES) %>%
#   mutate(SHARP50=(`75%`-`25%`)/1000,SHARP95=(`97.5%`-`2.5%`)/1000) %>%
#   group_by(ISLAND,DATECUT,METHOD) %>%
#   summarise(SHARP50=mean(SHARP50),SHARP95=mean(SHARP95))
# 
# marg.sharp = 0.013
# Fig_sharp = 
#   ggplot() +
#   geom_line(data=tmp.sharp,aes(x=DATECUT,y=SHARP95,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=top*marg.sharp*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=top*marg.sharp*.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=top*marg.sharp*.9,yend=0),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=top*marg.sharp*.9),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=top*marg.sharp*.9),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=top*marg.sharp*.9),label="E",size=3.5) +
#   
#   geom_point(data=timing,aes(x=min,y=top*marg.sharp),colour="white") +
#   geom_point(data=timing,aes(x=max,y=0),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="95%PI width (x1,000)",colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# # GROUP -------------------------------------------------------------------------
# 
# Fig4_A <- ggplotGrob(Fig_totobs + theme(legend.position="none")) # convert to gtable
# Fig4_B <- ggplotGrob(Fig_rmsd + theme(legend.position="none")) # convert to gtable
# Fig4_C <- ggplotGrob(Fig_sharp + theme(legend.position="none")) # convert to gtable
# 
# 
# Fig4_leg = mutate(TSUM_RMSD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                   METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                   METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                   ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
#   ggplot() +
#   geom_line(aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
#   scale_colour_manual(values=as.character(colPRIOR$col),
#                       labels=c("N-IP"=expression(paste("neither ",beta," nor ",rho)),
#                                "IP on R0"=expression(paste("only ",beta)),
#                                "IP on rho"=expression(paste("only ",rho)),
#                                "IP on both"=expression(paste("both ",beta," and ",rho)))) +
#   theme(legend.title=element_text(size=9),legend.text=element_text(size=8,lineheight=.6),legend.text.align=0) + 
#   labs(colour="Approach involves\ninformative priors on:")
# Fig4_leg = get_legend(Fig4_leg)
# 
# max.widths <- unit.pmax(Fig4_A$widths[1:4], Fig4_B$widths[1:4],Fig4_C$widths[1:4]) # calculate maximum widths
# Fig4_A$widths[1:4] <- max.widths
# Fig4_B$widths[1:4] <- max.widths
# Fig4_C$widths[1:4] <- max.widths
# 
# Fig4 = plot_grid(Fig4_A,Fig4_B,Fig4_C,
#                  ncol=1,labels=LETTERS[1:3],label_size=10)
# Fig4 = plot_grid(Fig4,Fig4_leg,rel_widths = c(3,.7))
# Fig4
# 
# ggsave(plot=Fig4,
#        filename="~/owncloud/zikachik_shared/manuscript_prediction/figures/Fig4.pdf",
#        width=18,height=16,unit="cm")
# 
# 
# 
# ##################################################################################################################
# # Figure with prediction of epidemic dynamics
# 
# TSUM_DYN = NULL
# for(i in unique(TSUM_F_NCASES$DATECUT)) {
#   for(j in 0:1) {
#     for(k in 0:1) {
#       tmp1 = filter(fullD1,DATE<=i) %>%
#         mutate(DATECUT=as.Date(i,origin="1970-01-01"),PRIOR_R0=j,PRIOR_rho=k,TYPE="ACTUAL") %>%
#         select(DATECUT,DATE,ISLAND,PRIOR_R0,PRIOR_rho,mean=NCASES,TYPE)
#       tmp2 = filter(TSUM_F_NCASES,DATECUT==i,PRIOR_R0==j,PRIOR_rho==k,mean>0) %>%
#         mutate(TYPE="PRED") %>%
#         select(DATECUT,DATE,ISLAND,PRIOR_R0,PRIOR_rho,mean,TYPE)
#       TSUM_DYN = bind_rows(TSUM_DYN,tmp1,tmp2)
#     }
#   }
# }
# 
# TSUM_DYN = arrange(TSUM_DYN,DATECUT,ISLAND,PRIOR_R0,PRIOR_rho) %>%
#   group_by(DATECUT,ISLAND,PRIOR_R0,PRIOR_rho) %>%
#   mutate(THRESHOLD=ifelse(ISLAND=="SAINTMARTIN",20,200),
#          ISTHRESHOLD=mean>THRESHOLD,
#          IST_M1=lag(ISTHRESHOLD,1),
#          IST_M2=lag(ISTHRESHOLD,2),
#          ITS_M012=ISTHRESHOLD+IST_M1+IST_M2,
#          ITS_M012=ifelse(is.na(ITS_M012),0,ITS_M012),
#          CUMITS=cumsum(ITS_M012),
#          START=as.Date(ifelse(ITS_M012==3,DATE,NA),origin="1970-01-01"),
#          END=as.Date(ifelse(CUMITS==max(CUMITS) & ITS_M012==1,DATE+7,NA),origin="1970-01-01"),
#          PEAK=as.Date(ifelse(mean==max(mean),DATE,NA),origin="1970-01-01")) %>%
#   summarise(START=min(START,na.rm=T),END=max(END,na.rm=T),PEAK=max(PEAK,na.rm=T))
# 
# TSUM_DYN = left_join(TSUM_DYN,timing) %>% 
#   mutate(START_REL=ifelse(start>=DATECUT,(START-start)/7,NA),
#          END_REL=ifelse(end>=DATECUT,(END-end)/7,NA),
#          PEAK_REL=ifelse(peak>=DATECUT,(PEAK-peak)/7,NA)) %>%
#   ungroup() %>%
#   mutate(METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#          METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#          METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#          METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#          METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#          ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# 
# # Epidemic start -------------------------------------------------------------------------
# 
# marg.dyn1 = 29
# Fig_epstart = ggplot() +
#   geom_line(data=TSUM_DYN,aes(x=DATECUT,y=START_REL,colour=METHOD),size=.8) +
#   
#   geom_point(data=timing,aes(x=min,y=marg.dyn1),colour="white") +
#   geom_point(data=timing,aes(x=max,y=-marg.dyn1),colour="white") +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=marg.dyn1*.9,yend=-marg.dyn*.8),alpha=.6) +
#   geom_label(data=timing,aes(x=start-1,y=marg.dyn1*.8),label="S",size=3.5) +
#   
#   geom_hline(yintercept=0,linetype=2,alpha=.6) +
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="Weeks",colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0),breaks=c(-20,-10,0,10,20),labels=c("-20","-10","0","+10","+20")) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# # Epidemic peak -------------------------------------------------------------------------
# marg.dyn2 = 32
# Fig_eppeak = ggplot() +
#   geom_line(data=TSUM_DYN,aes(x=DATECUT,y=PEAK_REL,colour=METHOD),size=.8) +
#   
#   geom_point(data=timing,aes(x=min,y=marg.dyn2),colour="white") +
#   geom_point(data=timing,aes(x=max,y=-marg.dyn2),colour="white") +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=marg.dyn2*.8,yend=-marg.dyn2*.8),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=marg.dyn2*.8,yend=-marg.dyn2*.8),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=marg.dyn2*.8),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=marg.dyn2*.8),label="P",size=3.5) +
#   
#   geom_hline(yintercept=0,linetype=2,alpha=.6) +
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="Weeks",colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0),breaks=c(-30,-20,-10,0,10,20,30),labels=c("-30","-20","-10","0","+10","+20","+30")) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# # Epidemic end -------------------------------------------------------------------------
# marg.dyn3 = 72
# Fig_epend = ggplot() +
#   geom_line(data=TSUM_DYN,aes(x=DATECUT,y=END_REL,colour=METHOD),size=.8) +
#   
#   geom_point(data=timing,aes(x=min,y=marg.dyn3),colour="white") +
#   geom_point(data=timing,aes(x=max,y=-marg.dyn3),colour="white") +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=marg.dyn3*.8,yend=-marg.dyn3*.8),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=marg.dyn3*.8,yend=-marg.dyn3*.8),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=marg.dyn3*.8,yend=-marg.dyn3*.8),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=marg.dyn3*.8),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=marg.dyn3*.8),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=marg.dyn3*.8),label="E",size=3.5) +
#   
#   geom_hline(yintercept=0,linetype=2,alpha=.6) +
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y="Weeks",colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0),breaks=c(-60,-40,-20,0,20,40,60),labels=c("-60","-40","-20","0","+20","+40","+60")) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# # GROUP -------------------------------------------------------------------------
# 
# Fig5_A <- ggplotGrob(Fig_epstart + theme(legend.position="none")) # convert to gtable
# Fig5_B <- ggplotGrob(Fig_eppeak + theme(legend.position="none")) # convert to gtable
# Fig5_C <- ggplotGrob(Fig_epend + theme(legend.position="none")) # convert to gtable
# 
# max.widths <- unit.pmax(Fig5_A$widths[1:4], Fig5_B$widths[1:4],Fig5_C$widths[1:4]) # calculate maximum widths
# Fig5_A$widths[1:4] <- max.widths
# Fig5_B$widths[1:4] <- max.widths
# Fig5_C$widths[1:4] <- max.widths
# 
# Fig5 = plot_grid(Fig5_A,Fig5_B,Fig5_C,
#                  ncol=1,labels=LETTERS[1:3],label_size=10)
# Fig5 = plot_grid(Fig5,Fig4_leg,rel_widths = c(3,.7))
# Fig5
# 
# ggsave(plot=Fig5,
#        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig5.pdf",
#        width=18,height=15,unit="cm")
# 
# 
# 
# ##################################################################################################################
# # Figure posteriors and KL-distance
# 
# # Evolution of R0 ----------------------------------------------------------------
# 
# tmp.r0 = mutate(TSUM_R0,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                 METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                 METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                 METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                 METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                 ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# tmpend.r0 = filter(tmp.r0,DATECUT==max(DATECUT)) %>%
#   group_by(ISLAND) %>%
#   mutate(squaremax=max(`97.5%`),squaremin=min(`2.5%`),squaremeanmax=max(mean),squaremeanmin=min(mean))
# 
# Fig7_A = 
#   ggplot() +
#   
#   geom_line(data=tmp.r0,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=3.6,yend=1),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=3.6,yend=1),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=3.6,yend=1),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=3.6),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=3.6),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=3.6),label="E",size=3.5) +
#   
#   geom_linerange(data=tmpend.r0,aes(x=DATECUT+40,ymax=`75%`,ymin=`25%`,colour=METHOD),position=position_dodge(40),size=1) +
#   geom_linerange(data=tmpend.r0,aes(x=DATECUT+40,ymax=`97.5%`,ymin=`2.5%`,colour=METHOD),position=position_dodge(40)) +
#   geom_rect(data=tmpend.r0,aes(xmin=DATECUT+40-30,xmax=DATECUT+40+30,ymin=squaremin-.07,ymax=squaremax+.07),fill=NA,colour="black",size=.1) +
#   
#   geom_point(data=timing,aes(x=min,y=1),colour="white") +
#   geom_point(data=timing,aes(x=max,y=1),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y=expression(R[0]),colour="Approach") +
#   
#   scale_y_continuous(limits=c(.9,3.9),expand=c(0,0),breaks=c(1,1.5,2,2.5,3,3.5)) +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# # KL div of R0 ----------------------------------------------------------------
# 
# tmp.klr0 = mutate(TSUM_KLD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                   METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                   METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                   ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# 
# Fig7_B =  ggplot() +
#   geom_line(data=tmp.klr0,aes(x=DATECUT,y=KLD_R0,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=36,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=36,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=36,yend=0),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=36),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=36),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=36),label="E",size=3.5) +
#   
#   geom_point(data=timing,aes(x=min,y=1),colour="white") +
#   geom_point(data=timing,aes(x=max,y=1),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y=expression(paste("K-L divergence of ",R[0])),colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0),limits=c(0,40)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# # Evolution of rho ----------------------------------------------------------------
# 
# tmp.rho = mutate(TSUM_rho,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
#                  METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
#                  METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
#                  METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
#                  METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
#                  ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN")))
# tmpend.rho = filter(tmp.rho,DATECUT==max(DATECUT)) %>%
#   group_by(ISLAND) %>%
#   mutate(squaremax=max(`97.5%`),squaremin=min(`2.5%`),squaremeanmax=max(mean),squaremeanmin=min(mean))
# 
# Fig7_C = ggplot() +
#   
#   geom_line(data=tmp.rho,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=.9,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=.9,yend=0),alpha=.6) +
#   
#   geom_rect(data=tmpend.rho,aes(xmin=DATECUT+40-30,xmax=DATECUT+40+30,ymin=squaremin-.02,ymax=squaremax+.02),fill=NA,colour="black",size=.1) +
#   
#   geom_label(data=timing,aes(x=start-1,y=.9),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=.9),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=.9),label="E",size=3.5) +
#   
#   geom_linerange(data=tmpend.rho,aes(x=DATECUT+40,ymax=`75%`,ymin=`25%`,colour=METHOD),position=position_dodge(40),size=1) +
#   geom_linerange(data=tmpend.rho,aes(x=DATECUT+40,ymax=`97.5%`,ymin=`2.5%`,colour=METHOD),position=position_dodge(40)) +
#   
#   geom_point(data=timing,aes(x=min,y=1),colour="white") +
#   geom_point(data=timing,aes(x=max,y=1),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y=expression(rho),colour="Approach") +
#   
#   scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,0.2,0.4,0.6,0.8,1),labels=scales::percent) +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# # KL div of rho ----------------------------------------------------------------
# 
# Fig7_D = ggplot() +
#   geom_line(data=tmp.klr0,aes(x=DATECUT,y=KLD_rho,colour=METHOD),size=.8) +
#   
#   geom_segment(data=timing,aes(x=start-1,xend=start-1,y=36,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=peak-1,xend=peak-1,y=36,yend=0),alpha=.6) +
#   geom_segment(data=timing,aes(x=end-1,xend=end-1,y=36,yend=0),alpha=.6) +
#   
#   geom_label(data=timing,aes(x=start-1,y=36),label="S",size=3.5) +
#   geom_label(data=timing,aes(x=peak-1,y=36),label="P",size=3.5) +
#   geom_label(data=timing,aes(x=end-1,y=36),label="E",size=3.5) +
#   
#   geom_point(data=timing,aes(x=min,y=1),colour="white") +
#   geom_point(data=timing,aes(x=max,y=1),colour="white") +
#   
#   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
#   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
#   labs(x="Cut-off date",y=expression(paste("K-L divergence of ",rho)),colour="Approach") +
#   scale_colour_manual(values=as.character(colPRIOR$col)) +
#   scale_y_continuous(expand=c(0,0),limits=c(0,40)) +
#   scale_x_date(expand=c(0,0),breaks=as.Date(c("2016-01-01","2016-04-01","2016-07-01","2016-10-01","2017-01-01")),
#                labels=c("Jan.\n2016","April\n","July\n","Oct.\n","Jan.\n2017")) 
# 
# 
# 
# 
# 
# Fig7_A <- ggplotGrob(Fig7_A + theme(legend.position="none")) # convert to gtable
# Fig7_B <- ggplotGrob(Fig7_B + theme(legend.position="none")) # convert to gtable
# Fig7_C <- ggplotGrob(Fig7_C + theme(legend.position="none")) # convert to gtable
# Fig7_D <- ggplotGrob(Fig7_D + theme(legend.position="none")) # convert to gtable
# 
# max.widths <- unit.pmax(Fig7_A$widths[1:3], Fig7_B$widths[1:3], Fig7_C$widths[1:3], Fig7_D$widths[1:3]) # calculate maximum widths
# Fig7_A$widths[1:3] <- max.widths
# Fig7_B$widths[1:3] <- max.widths
# Fig7_C$widths[1:3] <- max.widths
# Fig7_D$widths[1:3] <- max.widths
# 
# Fig7 = plot_grid(Fig7_A,Fig7_B,Fig7_C,Fig7_D,
#                  ncol=1,labels=LETTERS[1:4],label_size = 10)
# Fig7 = plot_grid(Fig7,Fig4_leg,rel_widths = c(3,.6))
# Fig7
# 
# ggsave(plot=Fig7,
#        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig7.pdf",
#        width=18,height=19,unit="cm")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # 
# # 
# # tmp = group_by(D1,ISLAND) %>%
# #   summarise(TT=sum(NCASES)/1000) 
# # TOTOBS= mutate(TSUM_TOTALOBSCASES,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #                METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #                METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #                METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #                METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #                ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18")
# # lims = data.frame(ISLAND=c("MARTINIQUE","MARTINIQUE","GUADELOUPE","GUADELOUPE","SAINTMARTIN","SAINTMARTIN"),
# #                   lim=c(0,220,0,220,0,22))
# # phases$phase.y.NOBS = c(200,200,20)
# # 
# # Fig3_C = ggplot() +
# #   geom_line(data=TOTOBS,aes(x=DATECUT,y=mean/1000,colour=METHOD),size=.8) +
# #   geom_segment(data=tmp,aes(y=TT,yend=TT),x=as.numeric(as.Date("2016-02-07")),xend=as.numeric(as.Date("2016-09-18")),linetype=2) +
# #   geom_point(data=lims,aes(x=as.Date("2016-01-31"),y=lim),colour="white",alpha=0) + 
# #   scale_y_continuous(expand=c(0,0)) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   # geom_segment(data=phases,aes(y=phase.y.NOBS),yend=0,x=as.numeric(as.Date("2016-02-07")),xend=as.numeric(as.Date("2016-02-07")),colour="grey") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.NOBS*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.NOBS*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.NOBS*1.05),label="P3",size=2) +
# #   labs(x="Cut-off date",y="Pred. total reported\ncases (x1000)",colour="Approach") +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # TRD = mutate(TSUM_TREND,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18")
# # phases$phase.y.TREND = c(1.1,1.1,1.1)
# # Fig3_D = ggplot() +
# #   geom_line(data=TRD,aes(x=DATECUT,y=1-TREND_5WEEKS,colour=METHOD),size=.8) +
# #   scale_y_continuous(labels= scales::percent,limits=c(0,1.21),expand=c(0,0),breaks=c(0,.25,.5,.75,1)) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y="Pred. downward trend\nin the next 5 weeks",colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.TREND*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.TREND*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.TREND*1.05),label="P3",size=2)  +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # require(grid) # for unit.pmax(), unit.list()
# # Fig3_A <- ggplotGrob(Fig3_A + theme(legend.position="none")) # convert to gtable
# # Fig3_B <- ggplotGrob(Fig3_B + theme(legend.position="none")) # convert to gtable
# # Fig3_C <- ggplotGrob(Fig3_C + theme(legend.position="none")) # convert to gtable
# # Fig3_D <- ggplotGrob(Fig3_D + theme(legend.position="none")) # convert to gtable
# # 
# # max.widths <- unit.pmax(Fig3_A$widths[1:3], Fig3_B$widths[1:3],Fig3_C$widths[1:3],Fig3_D$widths[1:3]) # calculate maximum widths
# # Fig3_A$widths[1:3] <- max.widths
# # Fig3_B$widths[1:3] <- max.widths
# # Fig3_C$widths[1:3] <- max.widths
# # Fig3_D$widths[1:3] <- max.widths
# # 
# # Fig3 = plot_grid(Fig2_A,Fig3_A,Fig3_B,Fig3_C,Fig3_D,
# #                  ncol=1,labels=LETTERS[1:6],label_size=10)
# # Fig3 = plot_grid(Fig3,Fig2_leg,rel_widths = c(3,.5))
# # Fig3
# # 
# # ggsave(plot=Fig3, 
# #        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig3.pdf",
# #        width=20,height=18,unit="cm")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # phases = rbind(
# #   data.frame(ISLAND="GUADELOUPE",phase.y=3100,phase1.x=as.Date("2016-03-27"),
# #              phase1.xend=as.Date("2016-05-15"),phase2.xend=as.Date("2016-07-10")),
# #   data.frame(ISLAND="MARTINIQUE",phase.y=3100,phase1.x=as.Date("2016-01-10"),
# #              phase1.xend=as.Date("2016-02-14"),phase2.xend=as.Date("2016-05-15")),
# #   data.frame(ISLAND="SAINTMARTIN",phase.y=300,phase1.x=as.Date("2016-04-03"),
# #              phase1.xend=as.Date("2016-05-22"),phase2.xend=as.Date("2016-07-24")))
# # phases$phase2.x=phases$phase1.xend+7
# # phases$phase3.x=phases$phase2.xend+7
# # phases$phase3.xend=as.Date("2016-10-02")
# # phases$phase1.n=phases$phase1.x + round(as.numeric(phases$phase1.xend-phases$phase1.x)/2)
# # phases$phase2.n=phases$phase2.x + round(as.numeric(phases$phase2.xend-phases$phase2.x)/2)
# # phases$phase3.n=phases$phase3.x + round(as.numeric(phases$phase3.xend-phases$phase3.x)/2)
# # phases$xmin = as.Date("2015-12-31")
# # 
# # Fig2_A =   ggplot() +
# #   geom_bar(data=D1,aes(x=DATE,y=NCASES),stat="identity") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y,yend=phase.y)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y,yend=phase.y)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y,yend=phase.y)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y*1.08),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y*1.08),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y*1.08),label="P3",size=2) +
# #   geom_point(data=phases,aes(x=xmin,y=phase.y*1.18),colour="white") +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + 
# #   labs(x="",y="N") +
# #   scale_x_date(expand=c(0,0)) +
# #   scale_y_continuous(expand=c(0,0))
# # 
# # 
# # 
# # plot_point_example <- function(prior_r0,prior_rho,datechoice,...) {
# #   ex_fit <- filter(TSUM_FIT,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
# #     left_join(datechoice) %>%
# #     filter(DATECUT==choice,DATE<="2016-10-02")
# #   ex_pred <- filter(TSUM_PRED,PRIOR_R0==prior_r0,PRIOR_rho==prior_rho) %>%
# #     left_join(datechoice) %>%
# #     filter(DATECUT==choice,DATE<="2016-10-02") %>%
# #     mutate(`97.5%`=ifelse(`97.5%`>max,max,`97.5%`),
# #            `75%`=ifelse(`75%`>max,max,`75%`)) %>%
# #     left_join(colPRIOR)
# #   ex_real <- left_join(D1,datechoice) %>%
# #     filter(DATE>choice,DATE<="2016-10-02") 
# #   g <- ggplot() +
# #     geom_ribbon(data=ex_pred,aes(x=DATE,ymin=`2.5%`,ymax=`97.5%`),alpha=.1) +
# #     geom_ribbon(data=ex_pred,aes(x=DATE,ymin=`25%`,ymax=`75%`),alpha=.2) +
# #     geom_point(data=ex_fit,aes(x=DATE,y=NCASES),shape=21,fill="black",size=.9) +
# #     geom_point(data=ex_real,aes(x=DATE,y=NCASES),shape=21,size=.9) +
# #     geom_line(data=ex_pred,aes(x=DATE,y=mean),colour=unique(as.character(ex_pred$col)),size=0.7) +
# #     geom_vline(data=datechoice,aes(xintercept=as.numeric(choice)+7),linetype=2) +
# #     geom_point(data=datechoice,aes(x=xmin,y=max),col="white") +
# #     facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #     scale_x_date(expand=c(0,0))  +
# #     scale_colour_discrete(guide=FALSE) +
# #     scale_y_continuous(expand=c(0,0)) +
# #     theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA),
# #                         strip.background=element_blank(),
# #                         axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
# #     labs(x="",y="N") 
# #   return(g)
# # }
# # 
# # plot_acc_example = function(datechoice,...) {
# #   if(length(datechoice)==1) datechoice = rep(datechoice,3)
# #   stopifnot(length(datechoice)==3)
# #   datechoice = cbind(phases,dc=datechoice)
# #   ex_acc <- left_join(TSUM_RMSD,datechoice) %>%
# #     mutate(dc=as.numeric(as.Date(as.character(dc))-DATECUT)) %>% 
# #     filter(dc>0) %>% filter(dc==min(dc))
# #   
# #   
# # }
# # 
# # 
# # dc1 = data.frame(ISLAND=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"),
# #                  choice=phases$phase1.xend,
# #                  max=c(14500,14500,1300),
# #                  xmin=as.Date("2015-12-31"))
# # Fig2_B = plot_point_example(prior_r0=0,prior_rho=0,datechoice=dc1)
# # Fig2_C = plot_point_example(prior_r0=1,prior_rho=0,datechoice=dc1)
# # Fig2_D = plot_point_example(prior_r0=0,prior_rho=1,datechoice=dc1)
# # Fig2_E = plot_point_example(prior_r0=1,prior_rho=1,datechoice=dc1)
# # 
# # 
# # 
# # Fig2 = plot_grid(Fig2_A,Fig2_B,Fig2_C,Fig2_D,Fig2_E,labels=LETTERS[1:6],ncol=1,label_size = 10)
# # Fig2_leg = mutate(TSUM_RMSD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #                   METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #                   METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #                   METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #                   ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   ggplot() +
# #   geom_line(aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
# #   scale_colour_manual(values=as.character(colPRIOR$col),
# #                       labels=c("N-IP"="Non-informative\npriors",
# #                                "IP on R0"=expression(paste("IP on ",beta)),
# #                                "IP on rho"=expression(paste("IP on ",rho)),
# #                                "IP on both"=expression(paste("IP on ",beta," and ",rho)))) +
# #   theme(legend.title=element_text(size=9),legend.text=element_text(size=8,lineheight=.8),legend.text.align=0) + labs(colour="Approach")
# # 
# # Fig2_leg = get_legend(Fig2_leg)
# # Fig2 = plot_grid(Fig2,Fig2_leg,rel_widths = c(3,.6))
# # Fig2
# # ggsave(plot=Fig2, 
# #        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig2.pdf",
# #        width=18,height=21,unit="cm")
# # 
# # # Fig 3: RMSD, PEAK, Observed cases
# # 
# # phases$phase.y.RMSD = c(13000,13000,1300)
# # phases$phase1.x.RMSD = as.Date("2016-01-31")
# # phases$phase1.n.RMSD=phases$phase1.x.RMSD + round(as.numeric(phases$phase1.xend-phases$phase1.x)/2)
# # tmp = mutate(TSUM_RMSD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18")
# # 
# # Fig3_A = ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
# #   scale_y_log10() +
# #   # geom_segment(data=phases,aes(y=phase.y.RMSD),yend=0,x=as.numeric(as.Date("2016-02-07")),xend=as.numeric(as.Date("2016-02-07")),colour="grey") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.RMSD,yend=phase.y.RMSD)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.RMSD,yend=phase.y.RMSD)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.RMSD,yend=phase.y.RMSD)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.RMSD*1.2),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.RMSD*1.2),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.RMSD*1.2),label="P3",size=2) +
# #   geom_point(data=phases,aes(x=phase1.x,y=phase.y*1.3),colour="white") +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y="Accuracy\n(RMSD, log scale)",colour="Approach") +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # tmp = mutate(TSUM_PRED,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18") %>%
# #   left_join(D1) %>%
# #   filter(!is.na(NCASES)) %>%
# #   select(DATECUT,ISLAND,METHOD,PRIOR_R0,PRIOR_rho,`2.5%`,`25%`,`75%`,`97.5%`,DATE,NCASES) %>%
# #   mutate(ACC50=ifelse(NCASES<=`75%` & NCASES>=`25%`,1,0),
# #          ACC95=ifelse(NCASES<=`97.5%` & NCASES>=`2.5%`,1,0),
# #          SHARP50=(`75%`-`25%`),
# #          SHARP95=(`97.5%`-`2.5%`)) %>%
# #   group_by(ISLAND,DATECUT,METHOD) %>%
# #   summarise(ACC50=mean(ACC50),
# #             ACC95=mean(ACC95),
# #             SHARP50=mean(SHARP50),
# #             SHARP95=mean(SHARP95))
# # phases$phase.y.SHARP = c(80000,50000,6000)
# # Fig3_B = ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=SHARP95,colour=METHOD),size=.8) +
# #   scale_y_log10(breaks=c(100,1000,10000),labels=c("100","1000","10000")) +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.SHARP,yend=phase.y.SHARP)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.SHARP,yend=phase.y.SHARP)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.SHARP,yend=phase.y.SHARP)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.SHARP*1.2),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.SHARP*1.2),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.SHARP*1.2),label="P3",size=2) +
# #   geom_point(data=phases,aes(x=phase1.x,y=phase.y*1.3),colour="white") +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y="Sharpness\n(95%CI width, log scale)",colour="Approach") +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # 
# # tmp = group_by(D1,ISLAND) %>%
# #   summarise(TT=sum(NCASES)/1000)
# # TOTOBS= mutate(TSUM_TOTALOBSCASES,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #                METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #                METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #                METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #                METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #                ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18")
# # lims = data.frame(ISLAND=c("MARTINIQUE","MARTINIQUE","GUADELOUPE","GUADELOUPE","SAINTMARTIN","SAINTMARTIN"),
# #                   lim=c(0,220,0,220,0,22))
# # phases$phase.y.NOBS = c(200,200,20)
# # 
# # Fig3_C = ggplot() +
# #   geom_line(data=TOTOBS,aes(x=DATECUT,y=mean/1000,colour=METHOD),size=.8) +
# #   geom_segment(data=tmp,aes(y=TT,yend=TT),x=as.numeric(as.Date("2016-02-07")),xend=as.numeric(as.Date("2016-09-18")),linetype=2) +
# #   geom_point(data=lims,aes(x=as.Date("2016-01-31"),y=lim),colour="white",alpha=0) + 
# #   scale_y_continuous(expand=c(0,0)) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   # geom_segment(data=phases,aes(y=phase.y.NOBS),yend=0,x=as.numeric(as.Date("2016-02-07")),xend=as.numeric(as.Date("2016-02-07")),colour="grey") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.NOBS,yend=phase.y.NOBS)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.NOBS*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.NOBS*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.NOBS*1.05),label="P3",size=2) +
# #   labs(x="Cut-off date",y="Pred. total reported\ncases (x1000)",colour="Approach") +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # TRD = mutate(TSUM_TREND,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<="2016-09-18")
# # phases$phase.y.TREND = c(1.1,1.1,1.1)
# # Fig3_D = ggplot() +
# #   geom_line(data=TRD,aes(x=DATECUT,y=1-TREND_5WEEKS,colour=METHOD),size=.8) +
# #   scale_y_continuous(labels= scales::percent,limits=c(0,1.21),expand=c(0,0),breaks=c(0,.25,.5,.75,1)) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels),scales="free_y") +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y="Pred. downward trend\nin the next 5 weeks",colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.TREND,yend=phase.y.TREND)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.TREND*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.TREND*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.TREND*1.05),label="P3",size=2)  +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # require(grid) # for unit.pmax(), unit.list()
# # Fig3_A <- ggplotGrob(Fig3_A + theme(legend.position="none")) # convert to gtable
# # Fig3_B <- ggplotGrob(Fig3_B + theme(legend.position="none")) # convert to gtable
# # Fig3_C <- ggplotGrob(Fig3_C + theme(legend.position="none")) # convert to gtable
# # Fig3_D <- ggplotGrob(Fig3_D + theme(legend.position="none")) # convert to gtable
# # 
# # max.widths <- unit.pmax(Fig3_A$widths[1:3], Fig3_B$widths[1:3],Fig3_C$widths[1:3],Fig3_D$widths[1:3]) # calculate maximum widths
# # Fig3_A$widths[1:3] <- max.widths
# # Fig3_B$widths[1:3] <- max.widths
# # Fig3_C$widths[1:3] <- max.widths
# # Fig3_D$widths[1:3] <- max.widths
# # 
# # Fig3 = plot_grid(Fig2_A,Fig3_A,Fig3_B,Fig3_C,Fig3_D,
# #                  ncol=1,labels=LETTERS[1:6],label_size=10)
# # Fig3 = plot_grid(Fig3,Fig2_leg,rel_widths = c(3,.5))
# # Fig3
# # 
# # ggsave(plot=Fig3, 
# #        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig3.pdf",
# #        width=20,height=18,unit="cm")
# # 
# # 
# # # Theoretical considerations
# # 
# # tmp = mutate(TSUM_R0,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) 
# # tmpend = filter(tmp,DATECUT=="2016-10-02") %>%
# #   group_by(ISLAND) %>%
# #   mutate(squaremax=max(`97.5%`),squaremin=min(`2.5%`),squaremeanmax=max(mean),squaremeanmin=min(mean))
# # phases$phase.y.R0 = rep(3.7,3)
# # 
# # Fig4_A = ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
# #   geom_linerange(data=tmpend,aes(x=DATECUT+30,ymax=`75%`,ymin=`25%`,colour=METHOD),position=position_dodge(30),size=1) +
# #   geom_linerange(data=tmpend,aes(x=DATECUT+30,ymax=`97.5%`,ymin=`2.5%`,colour=METHOD),position=position_dodge(30)) +
# #   scale_y_continuous(limits=c(.8,4.1),expand=c(0,0),breaks=c(1,1.5,2,2.5,3,3.5)) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y=expression(R[0]),colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.R0,yend=phase.y.R0)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.R0,yend=phase.y.R0)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.R0,yend=phase.y.R0)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.R0*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.R0*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.R0*1.05),label="P3",size=2) +
# #   geom_rect(data=tmpend,aes(xmin=DATECUT+30-20,xmax=DATECUT+30+20,ymin=squaremin-.1,ymax=squaremax+.1),fill=NA,colour="black",size=.1) +
# #   # geom_segment(data=tmpend,aes(x=DATECUT,xend=DATECUT+10,y=squaremeanmax,yend=squaremax+.1),size=.1) +
# #   # geom_segment(data=tmpend,aes(x=DATECUT,xend=DATECUT+10,y=squaremeanmin,yend=squaremin-.1),size=.1) +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # 
# # tmp = mutate(TSUM_KLD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<"2016-10-02") 
# # phases$phase.y.R0KLD <- rep(35,3)
# # Fig4_B = ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=KLD_R0,colour=METHOD),size=.8) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y="K-L divergence",colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.R0KLD*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.R0KLD*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.R0KLD*1.05),label="P3",size=2) +
# #   geom_rect(data=tmpend,aes(xmin=DATECUT+30-20,xmax=DATECUT+30+20,ymin=squaremin-.1,ymax=squaremax+.1),fill=NA,colour="white",size=.1) +
# #   scale_colour_manual(values=as.character(colPRIOR$col)) +
# #   scale_y_continuous(expand=c(0,0),limits=c(0,40))
# # 
# # 
# # 
# # require(grid) # for unit.pmax(), unit.list()
# # Fig4_A <- ggplotGrob(Fig4_A + theme(legend.position="none")) # convert to gtable
# # Fig4_B <- ggplotGrob(Fig4_B + theme(legend.position="none")) # convert to gtable
# # max.widths <- unit.pmax(Fig4_A$widths[1:3], Fig4_B$widths[1:3]) # calculate maximum widths
# # Fig4_A$widths[1:3] <- max.widths
# # Fig4_B$widths[1:3] <- max.widths
# # 
# # Fig4 = plot_grid(Fig4_A,Fig4_B,
# #                  ncol=1,labels=LETTERS[1:2],label_size = 10)
# # Fig4 = plot_grid(Fig4,Fig2_leg,rel_widths = c(3,.5))
# # Fig4
# # 
# # ggsave(plot=Fig4, 
# #        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig4.pdf",
# #        width=22,height=10,unit="cm")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # tmp = mutate(TSUM_rho,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) 
# # tmpend = filter(tmp,DATECUT=="2016-10-02") %>%
# #   group_by(ISLAND) %>%
# #   mutate(squaremax=max(`97.5%`),squaremin=min(`2.5%`),squaremeanmax=max(mean),squaremeanmin=min(mean))
# # phases$phase.y.rho = rep(.9,3)
# # 
# # Fig5_A = 
# #   ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=mean,colour=METHOD),size=.8) +
# #   geom_linerange(data=tmpend,aes(x=DATECUT+30,ymax=`75%`,ymin=`25%`,colour=METHOD),position=position_dodge(30),size=1) +
# #   geom_linerange(data=tmpend,aes(x=DATECUT+30,ymax=`97.5%`,ymin=`2.5%`,colour=METHOD),position=position_dodge(30)) +
# #   scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,.2,.4,.6,.8,1),labels=scales::percent) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   labs(x="Cut-off date",y=expression(rho),colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.rho,yend=phase.y.rho)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.rho,yend=phase.y.rho)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.rho,yend=phase.y.rho)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.rho*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.rho*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.rho*1.05),label="P3",size=2) +
# #   geom_rect(data=tmpend,aes(xmin=DATECUT+30-20,xmax=DATECUT+30+20,ymin=squaremin-.02,ymax=squaremax+.02),fill=NA,colour="black",size=.1) +
# #   # geom_segment(data=tmpend,aes(x=DATECUT,xend=DATECUT+10,y=squaremeanmax,yend=squaremax+.02),size=.1) +
# #   # geom_segment(data=tmpend,aes(x=DATECUT,xend=DATECUT+10,y=squaremeanmin,yend=squaremin-.02),size=.1) +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # 
# # tmp = mutate(TSUM_KLD,METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==0,"N-IP",NA),
# #              METHOD=ifelse(PRIOR_R0==0 & PRIOR_rho==1,"IP on rho",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==0,"IP on R0",METHOD),
# #              METHOD=ifelse(PRIOR_R0==1 & PRIOR_rho==1,"IP on both",METHOD),
# #              METHOD=factor(METHOD,levels=c("N-IP","IP on R0","IP on rho","IP on both")),
# #              ISLAND=factor(ISLAND,levels=c("GUADELOUPE","MARTINIQUE","SAINTMARTIN"))) %>%
# #   filter(DATECUT<"2016-10-02") 
# # phases$phase.y.R0KLD <- rep(40,3)
# # Fig5_B = ggplot() +
# #   geom_line(data=tmp,aes(x=DATECUT,y=KLD_rho,colour=METHOD),size=.8) +
# #   facet_wrap(~ISLAND,labeller=labeller(ISLAND=island_labels)) +
# #   theme_bw2() + theme(panel.border=element_rect(colour="black",fill=NA)) +
# #   scale_y_continuous(limits=c(0,45),expand=c(0,0)) +
# #   labs(x="Cut-off date",y="K-L divergence",colour="Approach") +
# #   geom_segment(data=phases,aes(x=phase1.x,xend=phase1.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_segment(data=phases,aes(x=phase2.x,xend=phase2.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_segment(data=phases,aes(x=phase3.x,xend=phase3.xend,y=phase.y.R0KLD,yend=phase.y.R0KLD)) +
# #   geom_text(data=phases,aes(x=phase1.n,y=phase.y.R0KLD*1.05),label="P1",size=2) +
# #   geom_text(data=phases,aes(x=phase2.n,y=phase.y.R0KLD*1.05),label="P2",size=2) +
# #   geom_text(data=phases,aes(x=phase3.n,y=phase.y.R0KLD*1.05),label="P3",size=2) +
# #   geom_rect(data=tmpend,aes(xmin=DATECUT+30-20,xmax=DATECUT+30+20,ymin=squaremin-.1,ymax=squaremax+.1),fill=NA,colour="white",size=.1) +
# #   scale_colour_manual(values=as.character(colPRIOR$col))
# # 
# # 
# # 
# # require(grid) # for unit.pmax(), unit.list()
# # Fig5_A <- ggplotGrob(Fig5_A + theme(legend.position="none")) # convert to gtable
# # Fig5_B <- ggplotGrob(Fig5_B + theme(legend.position="none")) # convert to gtable
# # max.widths <- unit.pmax(Fig5_A$widths[1:3], Fig5_B$widths[1:3]) # calculate maximum widths
# # Fig5_A$widths[1:3] <- max.widths
# # Fig5_B$widths[1:3] <- max.widths
# # 
# # Fig5 = plot_grid(Fig5_A,Fig5_B,
# #                  ncol=1,labels=LETTERS[1:2],label_size = 10)
# # Fig5 = plot_grid(Fig5,Fig2_leg,rel_widths = c(3,.5))
# # Fig5
# # 
# # ggsave(plot=Fig5, 
# #        filename="D:/julien/owncloud/zikachik_shared/manuscript_prediction/figures/Fig5.pdf",
# #        width=22,height=10,unit="cm")

