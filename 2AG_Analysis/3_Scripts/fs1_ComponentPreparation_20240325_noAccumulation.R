#INPUTS-----------
mainwd <- "path_to_folder/" # change to path to folder!
mainData <- "2_CreatedData/A_PreprocessedData.csv"
fitDir <- "2_CreatedData/B_ControlFits_20240325_NoSignalingPool/"
fitData <- "FOCEi_result_NN2A_RE.Rdata"
output_evalPlot <- "2_CreatedData/B_ControlFits_20240325_NoSignalPool_eval.png"

# fixed information
fixedPars <- c(m_DAG = 1.12923469444637*10^-12,
               N_N2A = 4*10^4,
               N_HEK = 3.5*10^4,
               
               n_vesicle = 65, 
               N_2AG_vesicle = 2048, 
               N_avogadro = 6.022*10^23,
               
               m_pool_ss = 3.53809807438018 * 10^-12, # value is total (all N2A cell)
               m_pool_inhibited = 7.37572418275198 * 10^-14, # value is total (all N2A cell)
               
               m_signal_ss = 8.80461295845166 * 10^-15, # value is total (all N2A cell)
               
               Vd = 2.65*10^-15,
               Smax = 2.2,
               Ki = 4.17*10^-5)

#LIBRARIES----------
library(dplyr)
library(rxode2)
library(reshape2)
library(ggplot2)
library(cowplot)
set.seed(436)

#FUNCTIONS---------
appendFixedPars <- function(target_df, source_vec, par_names=NULL){
  for(i in c(1:length(par_names))){
    target_df[,par_names[i]] <- source_vec[which(names(source_vec)==par_names[i])]
  }
  
  return(target_df)
}

#MODELS-------------
# postATP_model
postATP_model <- rxode2({
  # random effect (where relevant)
  N_N2A_eff = N_N2A * exp(eta.N_N2A)
  N_HEK_eff = N_HEK 
  
  # state initiation
  m_extracellular(0) = 0
  m_intracellular(0) = 0
  
  # ODEs
  d/dt(m_extracellular) = k_release * N_N2A_eff * m_pool_ss_N2A - (k_absorption * N_HEK_eff + k_removal) * m_extracellular
  d/dt(m_intracellular) = k_absorption * m_extracellular - k_degradation * m_intracellular
  
  # post-calculation
  C = m_intracellular / Vd
  S = Smax/N_HEK * C / (C + Ki) * N_HEK_eff
})

#READS------------
# normalized main data
mainData_normalized <- read.csv(paste0(mainwd, mainData), header=T)
mainData_normalized$dv <- mainData_normalized$dv - min(subset(mainData_normalized, compound=='DMSO')$dv)

# collect fitted parameters
load(paste0(mainwd, fitDir, fitData))
fitPars <- fitRes$parFixed[,1] %>% as.numeric() %>% exp()
names(fitPars) <- gsub("t[.]", "", rownames(fitRes$parFixed))
fitPars <- fitPars[1:3]

# collect random effect matrix
omega_matrix <- fitRes$omega
# rm("fitRes")

#DERIVATION----------
# standardize all value to values per-N2A cells
fixedPars['m_pool_ss_N2A'] <- fixedPars['m_pool_ss'] / fixedPars['N_N2A']
fixedPars['m_pool_inhibited_N2A'] <- fixedPars['m_pool_inhibited'] / fixedPars['N_N2A']
fixedPars['m_dose_N2A'] <- fixedPars['n_vesicle'] * fixedPars['N_2AG_vesicle'] / fixedPars['N_avogadro']

# direct derivation
fixedPars['k_release'] <- fixedPars['m_dose_N2A'] / (30*60) / fixedPars['m_pool_ss_N2A'] # 1/s

#PARAMETER ISOLATION--------------
# preparation
core_parameters <- c(fixedPars[c("N_HEK", "N_N2A", "m_pool_ss_N2A", "k_release", "Smax", "Ki", "Vd")],
                     fixedPars[grepl("k_", names(fixedPars))],
                     fitPars[grepl("k_", names(fitPars))])
parameter_multiplier <- replicate(length(core_parameters), 1)

# apply multiplier
core_parameters <- core_parameters * parameter_multiplier

#PREPARATION-------------
evMat_postATP <- data.frame(time = seq(0, 0.5 * 60 * 60 , 5)) 

#SIMULATIONS--------------
# Post-ATP simulation
simRes_postATP <- postATP_model$solve(core_parameters, evMat_postATP, omega=omega_matrix, nSub=500) %>%
  data.frame()

#PLOT PREPARATION------------
# Post-ATP
simRes_postATP_summary <- dplyr::select(simRes_postATP, sim.id, time, S) %>%
  dcast(time ~ sim.id, value.var='S')
simRes_postATP_summary <- apply(simRes_postATP_summary[,2:ncol(simRes_postATP_summary)], 1,
                                function(x) quantile(x, c(0.5, 0.05, 0.95))) %>% t() %>%
  data.frame() %>% rename("median"=1, "lower"=2, "upper"=3) %>%
  mutate(time = subset(simRes_postATP, sim.id==1)$time)

# obspred 
obspred_data <- data.frame(fitRes) %>% filter(TIME > 0)

#PLOTTING----------
pltPostATP <- ggplot(simRes_postATP_summary, aes(x=time, y=median))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), fill='lightblue', alpha=0.45)+
  geom_point(data=subset(mainData_normalized, compound=="DMSO"), aes(x=time, y=dv), shape=1)+
  geom_line()+theme_bw()+
  scale_x_continuous(breaks=seq(0, max(simRes_postATP_summary$time), 600))+
  xlab("Time (seconds)")+ylab("Signal (a.u.)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face='bold'))

plt2 <- ggplot(obspred_data, aes(x=PRED, y=DV))+
  geom_point(shape=1)+theme_bw()+
  geom_abline(intercept=0)+
  xlab("Prediction")+ylab("Observation")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face='bold'))

plt_combination <- plot_grid(pltPostATP, plt2, ncol=2, labels=c("a)", "b)"), scale=0.95)

#OUTPUT--------------
ggsave(filename=paste0(mainwd, output_evalPlot), plot=plt_combination,
       device='png', width=60*4.5, height=25*4.5, units='mm', dpi=350)