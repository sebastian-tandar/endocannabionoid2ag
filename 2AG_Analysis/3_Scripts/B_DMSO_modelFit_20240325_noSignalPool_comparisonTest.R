#INPUTS-----------
mainwd <- "path_to_folder/" # change to path to folder!
mainData <- "2_CreatedData/A_PreprocessedData.csv"
outputwd <- "2_CreatedData/B_ControlFits_20240325_NoSignalingPool/"

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

#LIBRARIES-----------
library(dplyr)
library(ggplot2)
library(nlmixr2)
library(pso)
library(rxode2)
library(cowplot)
set.seed(436)

#FUNCTIONS-------------
appendFixedPars <- function(target_df, source_vec, par_names=NULL){
  for(i in c(1:length(par_names))){
    target_df[,par_names[i]] <- source_vec[which(names(source_vec)==par_names[i])]
  }
  
  return(target_df)
}
objFn_PSO <- function(pars, par_names, dataset){
  # parameter name adjustment
  names(pars) <- par_names
  
  # run simulation
  sim_res <- signalModel_rxode$solve(pars, dataset) %>% data.frame()
  
  # calculate error
  err <- sum((sim_res$S - dataset$dv)^2)
  
  # return value
  if(err < grandErr){
    grandErr <<- err
    print(err)
  }
  return(err)
}

#MODELS-------------
signalModel_nlmixr <- function(){
  ini({
    # fixed effect parameters
    t.k_removal = psoPrefit_parameters['t.k_removal'] %>% log()
    t.k_absorption = psoPrefit_parameters['t.k_absorption'] %>% log()
    t.k_degradation = psoPrefit_parameters['t.k_degradation'] %>% log()
    
    # random effect parameter
    eta.N_N2A ~ 0.1
    
    # residual effect
    add.err = 0.1
  })
  model({
    # parameter pre-calculation
    k_removal = exp(t.k_removal)
    k_absorption = exp(t.k_absorption)
    k_degradation = exp(t.k_degradation)
    
    # random effect (where relevant)
    N_N2A_eff = N_N2A * exp(eta.N_N2A)
    N_HEK_eff = N_HEK 
    # N_HEK_eff = N_HEK
    
    # state initiation
    m_extracellular(0) = 0
    m_intracellular(0) = 0
    
    # ODEs
    d/dt(m_extracellular) = k_release * N_N2A_eff * m_pool_ss_N2A - (k_absorption * N_HEK_eff + k_removal) * m_extracellular
    d/dt(m_intracellular) = k_absorption * m_extracellular - k_degradation * m_intracellular
    
    # post-calculation
    C = m_intracellular / Vd
    S = Smax/N_HEK * C / (C + Ki) * N_HEK_eff
    S ~ add(add.err)
  })
}
signalModel_rxode <- rxode2({
  # parameter pre-calculation
  k_removal = exp(t.k_removal)
  k_absorption = exp(t.k_absorption)
  k_degradation = exp(t.k_degradation)
  
  # random effect (where relevant)
  N_N2A_eff = N_N2A
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
mainData <- read.csv(paste0(mainwd, mainData), header=T) %>% 
  filter(compound=="DMSO") %>% 
  mutate(dv = dv - min(dv)) %>%  # normalize bottom to zero
  dplyr::select(-setID, -compound)

#PRE-CALCULATION-------------
# standardize all value to values per-N2A cells
fixedPars['m_pool_ss_N2A'] <- fixedPars['m_pool_ss'] / fixedPars['N_N2A']
fixedPars['m_pool_inhibited_N2A'] <- fixedPars['m_pool_inhibited'] / fixedPars['N_N2A']
fixedPars['m_dose_N2A'] <- fixedPars['n_vesicle'] * fixedPars['N_2AG_vesicle'] / fixedPars['N_avogadro']

# direct derivation
fixedPars['k_release'] <- fixedPars['m_dose_N2A'] / (30*60) / fixedPars['m_pool_ss_N2A'] # 1/s

#PREPARATION------------
# append fixed parameters to main data
mainData <- appendFixedPars(mainData, fixedPars, 
                            par_names = c("N_HEK", "N_N2A", 'k_release', 
                                          'm_pool_ss_N2A',
                                          "Vd", "Smax", "Ki"))

# parameter search boundaries
initPars <- c(t.k_removal = 5*10^-12, t.k_absorption = 3*10^-4, t.k_degradation = 1*10^-2) %>% log()
lowerBound <- c(t.k_removal = 10^-20, t.k_absorption = 10^-20, t.k_degradation = 10^-20) %>% log()
upperBound <- c(t.k_removal = 1, t.k_absorption = 1, t.k_degradation = 1) %>% log()

# #MAIN-------------------
# A | PSO Prefit
grandErr <- 10^10
psoPrefit <- psoptim(initPars, objFn_PSO,
                     par_names = names(initPars), dataset=mainData,
                     lower=lowerBound, upper=upperBound,
                     control=list(maxit=100))

# grab PSO result
psoPrefit_parameters <- psoPrefit$par[1:length(initPars)] %>% unlist() %>% exp()
names(psoPrefit_parameters) <- names(initPars)

# B | Main Fit
fitRes <- nlmixr2(signalModel_nlmixr, mainData, est='focei', control=list(outerOpt='bobyqa'))

# grab FOCEi result
fitRes_pars <- fitRes$parFixed[,1] %>% as.numeric() %>% exp()
names(fitRes_pars) <- rownames(fitRes$parFixed)
fitRes_pars <- fitRes_pars[names(fitRes_pars) %in% names(initPars)]

#OUTPUT-------------
dir.create(paste0(mainwd, outputwd))
save(psoPrefit, file=paste0(mainwd, outputwd, "PSO_result_NN2A_RE.Rdata"))
save(fitRes, file=paste0(mainwd, outputwd, "FOCEi_result_NN2A_RE.Rdata"))

# #SIMULATION--------------
# prepare event matrix
evMat <- data.frame(time = seq(0, 1800, 1)) %>%
  appendFixedPars(fixedPars,
                  par_names = c("N_HEK", "N_N2A", 'm_dose_N2A', 'k_release',
                                "m_pool_ss_N2A", 
                                "Vd", "Smax", "Ki"))

# run simulation for PSO result
psoPrefit_parameters <- exp(initPars)
simRes_pso <- signalModel_rxode$solve(log(psoPrefit_parameters), evMat) %>% data.frame()

# typical ID simulation
simRes_typical <- signalModel_rxode$solve(log(fitRes_pars), evMat) %>% data.frame()

#PLOTTING------------
pltPSO <- ggplot(data.frame(simRes_pso), aes(x=time, y=S))+
  geom_line()+theme_bw()+
  geom_point(data=mainData, aes(x=time, y=dv), shape=1)+
  geom_line(data=simRes_typical, aes(x=time, y=S), col='blue')
# pltPSO
pltNLMIXR <- ggplot(data.frame(fitRes), aes(x=TIME, y=DV))+
  geom_point(shape=1)+theme_bw()+
  geom_line(aes(x=TIME, y=IPRED, group=ID))+
  geom_line(data=simRes_typical, aes(x=time, y=S), col='blue', linewidth=1.2)

plot_combined <- plot_grid(pltPSO, pltNLMIXR, ncol=1, labels=c("A)", "B)"), align='hv', scale=0.95)
plot_combined