#LIBRARIES-------------
library(dplyr)
library(rxode2)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)

#INPUTS------------
mainwd <- "C:/Users/sebas/Documents/PROJECTS/2024_2GAFluorescence/3_Scripts/Signaling2AG/"
# mainwd <- "/home/project2AG/Signaling2AG/"
simComponents <- "SimulationPackage.Rdata"
fitData <- "FOCEi_result_NN2A_RE.Rdata"

#MODELS-----------
preATP_model <- rxode2({
  # state initiation
  m_pool(0) = m_pool_ss_N2A
  m_signal(0) = m_signal_ss_N2A * signal_fraction
  
  # ODEs
  d/dt(m_pool) = k_production * m_DAG_N2A - (k_metabolism + k_fill) * m_pool + k_sink * m_signal
  d/dt(m_signal) = k_fill * m_pool - k_sink * m_signal
})

postATP_model <- rxode2({
  # random effect (where relevant)
  N_N2A_eff = N_N2A * exp(eta.N_N2A)
  N_HEK_eff = N_HEK
  
  # state initiation
  m_pool(0) = m_pool_ss_N2A
  m_extracellular(0) = m_signal_ss_N2A * N_N2A_eff * signal_fraction
  m_intracellular(0) = 0
  
  # ODEs
  d/dt(m_pool) = k_production * m_DAG_N2A - (k_metabolism + k_fill) * m_pool
  d/dt(m_extracellular) = k_fill * N_N2A_eff * m_pool - (k_absorption * N_HEK_eff + k_removal) * m_extracellular
  d/dt(m_intracellular) = k_absorption * m_extracellular - k_degradation * m_intracellular
  
  # post-calculation
  C = m_intracellular / Vd
  S = Smax/N_HEK * C / (C + Ki) * N_HEK_eff
})


#INPUT PREPARATION----------
load(paste0(mainwd, fitData))
load(paste0(mainwd, simComponents))
setID_selection <- simPackage[['data']]$setID %>% unique()
pars <- data.frame(names = names(simPackage[['core_parameters']]),
                   value = simPackage[['core_parameters']]) %>% distinct() 

#FUNCTIONS----------------
runSimulation_preATP <- function(multiplier, signal_fraction, simulation_tmax = 180, simulation_dt = 0.1, as_fraction=T){
  #parse values
  applied_pars <- pars$value * multiplier
  names(applied_pars) <- pars$names
  
  # back-derivation
  back_derived <- c(m_pool_ss_N2A = unname(applied_pars['k_production'] * applied_pars['m_DAG_N2A'] / applied_pars['k_metabolism']))
  back_derived['m_signal_ss_N2A'] <- unname(applied_pars['k_fill'] * back_derived['m_pool_ss_N2A'] / applied_pars['k_sink'])
  back_derived['f_burst'] <- unname(back_derived['m_signal_ss_N2A'] / back_derived['m_pool_ss_N2A'])
  
  # event matrix
  evMat_preATP <- data.frame(time = seq(0, simulation_tmax * 60 , simulation_dt*60))
  
  # simulation
  simRes_preATP <- preATP_model$solve(c(signal_fraction = signal_fraction, applied_pars, back_derived), evMat_preATP) %>% 
    data.frame() %>% dplyr::select(time, m_pool, m_signal) %>% mutate(time = time / 60)
  
  if(as_fraction){
    # Pre-ATP | Normalized and non-normalized
    simRes_preATP_normalized <- mutate(simRes_preATP,
                                       m_pool = m_pool / back_derived['m_pool_ss_N2A'],
                                       m_signal = m_signal / back_derived['m_signal_ss_N2A']) %>%
      melt(id.vars='time')
    
    # plotting
    pltPreATP <- ggplot(simRes_preATP_normalized, aes(x=time, y=value*100))+
      geom_line()+facet_wrap(~variable, ncol=1)+theme_bw()+
      scale_x_continuous(breaks=seq(0, max(simRes_preATP$time), 60))+
      xlab("Time (minutes)")+ylab("%Saturation")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14, face='bold'),
            strip.text=element_text(size=12, face='bold'))
  }else{
    simRes_preATP <- melt(simRes_preATP, id.vars='time')
    pltPreATP <- ggplot(simRes_preATP, aes(x=time, y=value))+
      geom_line()+facet_wrap(~variable, ncol=1, scales='free_y')+theme_bw()+
      scale_x_continuous(breaks=seq(0, max(simRes_preATP$time), 60))+
      xlab("Time (minutes)")+ylab("2AG (mol)")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14, face='bold'),
            strip.text=element_text(size=12, face='bold'))
  }
  
  return(pltPreATP)
}
runSimulation_postATP <- function(multiplier, signal_fraction=1, simulation_tmax = 30, simulation_dt = 0.1, data_set = NULL, show_concentration=T){
  #parse values
  applied_pars <- pars$value * multiplier
  names(applied_pars) <- pars$names
  
  # back-derivation
  back_derived <- c(m_pool_ss_N2A = unname(applied_pars['k_production'] * applied_pars['m_DAG_N2A'] / applied_pars['k_metabolism']))
  back_derived['m_signal_ss_N2A'] <- unname(applied_pars['k_fill'] * back_derived['m_pool_ss_N2A'] / applied_pars['k_sink'])
  back_derived['f_burst'] <- unname(back_derived['m_signal_ss_N2A'] / back_derived['m_pool_ss_N2A'])
  
  # event matrix
  evMat_postATP <- data.frame(time = seq(0, simulation_tmax * 60 , simulation_dt*60))
  
  # Post-ATP simulation
  simRes_postATP <- postATP_model$solve(c(signal_fraction = signal_fraction, applied_pars, back_derived), evMat_postATP, omega=fitRes$omega, nSub=500) %>%
    data.frame()
  
  # post-processing - 2AG
  simRes_postATP_summary <- dplyr::select(simRes_postATP, sim.id, time, S) %>%
    dcast(time ~ sim.id, value.var='S')
  simRes_postATP_summary <- apply(simRes_postATP_summary[,2:ncol(simRes_postATP_summary)], 1,
                                  function(x) quantile(x, c(0.5, 0.05, 0.95))) %>% t() %>%
    data.frame() %>% rename("median"=1, "lower"=2, "upper"=3) %>%
    mutate(time = subset(simRes_postATP, sim.id==1)$time / 60)
  
  # plotting
  pltPostATP <- ggplot(simRes_postATP_summary, aes(x=time, y=median))+
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper), fill='lightblue', alpha=0.45)+
    geom_line()+theme_bw()+
    xlab("Time (minutes)")+ylab("Signal (a.u.)")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14, face='bold'))
  
  if(!is.null(data_set)){
    pltPostATP <- pltPostATP +
      geom_point(data=data_set, aes(x=time, y=dv), shape=1)
  }
  
  # concentration plotting
  if(show_concentration){
    simRes_postATP_summary_C <- dplyr::select(simRes_postATP, sim.id, time, C) %>%
      dcast(time ~ sim.id, value.var='C')
    simRes_postATP_summary_C <- apply(simRes_postATP_summary_C[,2:ncol(simRes_postATP_summary_C)], 1,
                                      function(x) quantile(x, c(0.5, 0.05, 0.95))) %>% t() %>%
      data.frame() %>% rename("median"=1, "lower"=2, "upper"=3) %>%
      mutate(time = subset(simRes_postATP, sim.id==1)$time / 60)
    
    pltPostATP_C <- ggplot(simRes_postATP_summary_C, aes(x=time, y=median / (10^-6)))+
      geom_ribbon(aes(x=time, ymin=lower/ (10^-6), ymax=upper/ (10^-6)), fill='lightblue', alpha=0.45)+
      geom_line()+theme_bw()+
      xlab("")+ylab("[2AG_HEK] (uM)")+
      theme(axis.text.y=element_text(size=12),
            axis.text.x=element_blank(),
            axis.title.y=element_text(size=14, face='bold'),
            axis.title.x=element_blank())
    plot_postATP_combined <- plot_grid(pltPostATP_C, pltPostATP, ncol=1, align='hv')
  }else{
    plot_postATP_combined <- pltPostATP
  }
  
  return(plot_postATP_combined)
}

#TROUBLESHOOTING----------------
# runSimulation_postATP(multiplier=1, signal_fraction=1, simulation_tmax = 30, simulation_dt = 0.1, data_set = NULL, show_concentration=T)