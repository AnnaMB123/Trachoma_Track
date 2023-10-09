library(iterators)
library(foreach)
library(doParallel)
library(truncdist)

source("Trachoma_Track_functions.R") #Standard trachoma model functions, note this is with children aged 0-9, can edit if need to fit for children aged 1-9
inits<-(read.csv("/Users/annaborlase/CIFF_trachoma/Inits.csv"))

parameters <- as.list(c(N=2500, #Population size
                        av_I_duration=2,av_ID_duration=28,min_ID=2,ag_ID=-0.00012, aq_ID=-0.039, #Parameters relating to duration of infection period, including ID period
                        av_D_duration=43,min_D=1, ag_D=-0.00098, aq_D=-0.078,  #Parameters relating to duration of disease period
                        v_1=1,v_2=2.6,phi=1.4,epsilon=0.5,#Parameters relating to lambda function- calculating force of infection
                        #Parameters relating to MDA
                        MDA_Cov=0.8, #This is standard, note-varied for Rombo analysis
                        MDA_Eff= 0.85, # Efficacy of MDA treatment
                        rho=0, #Vary in sensitivity analysis
                        nweeks_year=52,
                        babiesMaxAge=0.5, #Note this is years
                        youngChildMaxAge=9,#Note this is years
                        olderChildMaxAge=15, #Note this is years
                        b1=1,#this relates to bacterial load function
                        ep2=0.114)#this relates to bacterial load function
)


#Demography parameters
demog<-as.list(c(tau=1/(40*52),max_age=70*52,mean_age=20*52)) #tau= death rate in weeks^-1

nsim<-20 #Number of simulations
Target_TF_1_9<-20 ## This is just for setting initial number of infections etc, doesn't really matter if have long enough burnin i.e. more than 70years
Start_date<-as.Date("2018-01-01")
End_date<-as.Date("2030-12-31")
MDA_dates<-as.Date(c("2020-02-01","2021-05-01","2022-05-01")) #Can specify MDA dates here
burnin<-200*52
initialise_pop_and_seed(parameters = parameters, demog=demog, inits=inits, Target_TF_1_9 = Target_TF_1_9)
sim_params<-create_sim_params(nsim=nsim, burnin=burnin, Start_date=Start_date,End_date=End_date,MDA_dates=MDA_dates)
MDA_times<-MDA_dates_as_timestep(MDA_dates = MDA_dates, Start_date = Start_date, burnin = burnin)
Tx_mat<-Tx_matrix(params=parameters,sim_param=sim_params)

#MDA_times<- 99999999999 #Use this if want no MDA
set.seed(5)
bet<-runif(nsim,0.1,0.14) #If want to manually set beta

#Run simulations and store output

cores=detectCores()
cl<-makeCluster(cores[1]-1, type="FORK") #
#Parallel runs#
registerDoParallel(cl)

data_store_all_sim<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  
  data_store_all_sim[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:sim_params$timesim,demog=demog,bet=bet[i], MDA_times = MDA_times,
                                       Tx_mat = Tx_mat)
}
end_time<-Sys.time()
end_time-start_time
stopCluster(cl)

Summary_output_df<-extract_output_alltimes(data_store_all_sim = data_store_all_sim, sim_params = sim_params) #

library(ggplot2) 
Mean_inf_plot_children<- ggplot()+
  geom_line(aes(x=Summary_output_df$Time/52,y=Summary_output_df$Mean_inf_children),colour="red")+
  xlab("Time (years)")+
  ylab("Mean Infection Prev Children 0-9")
Mean_inf_plot_children

Mean_dz_plot_children<- ggplot()+
  geom_line(aes(x=Summary_output_df$Time/52,y=Summary_output_df$Mean_dz_children),colour="blue")+
  xlab("Time (years)")+
  ylab("Mean Disease Prev Children 0-9")
Mean_dz_plot_children


########################################################## Some alternative functions ##############################################

### Simulating topical treatment #####
### Not normally simulated as distinct from MDA (babies always given topical), but used for Tanzania analysis 

#Extra functions for doing topical MDA, determined by identifying people who have active trachoma (D==1) and efficacty of topical treatment
do_topical<-function(params,Age,IndD){
  Topical_Eff=0.73 #Bowman 2000
  Se_TF=0.915
  D<-which(IndD==1)
  treated<-D[which(runif(D)<Se_TF)]
  cured<-treated[which(runif(treated)< Topical_Eff)]
  treated_cured_topical<-c(cured)
  return(treated_cured_topical)
  
}

#This is time step in which topical treatment occurs#
Topical_timestep <- function(vals,params,demog){
  
  #Assign names to variables just so it's easier to track #
  
  IndI<-vals[[1]]
  IndD<-vals[[2]]
  No_Inf<-vals[[3]]
  T_latent<-vals[[4]]
  T_ID<-vals[[5]]
  T_D<-vals[[6]]
  Ind_latent<-vals[[7]]
  Ind_ID_period_base<-vals[[8]]
  Ind_D_period_base<-vals[[9]]
  bact_load<-vals[[10]]
  Age<-vals[[11]]
  
  #Identify who is treated AND cured
  treated_cured_topical<-do_topical(params=params,Age=Age,IndD=IndD)
  
  #Set treated/cured indivs infection status and bacterial load to 0
  IndI[treated_cured_topical]<-0 #clear infection they become I=0
  bact_load[treated_cured_topical]<-0 #clear disease they become D=0
  
  
  vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
  names(vals)<-c("IndI","IndD","No_Inf","T_latent","T_ID","T_D","Ind_latent","Ind_ID_period_base",
                 "Ind_D_period_base","bact_load","Age")
  
  
  return(vals)
} 

#This is alternative function to run simulation which includes separate topical treatment only
sim_Ind_MDA<-function(params,vals,timesim,demog,bet, MDA_times, Topical_times, Tx_mat){
  
  N = params$N
  N_Children_ages_0_9<-array()
  
  N_True_Infected_children_0_9<-array()
  True_Prev_Infection_children_0_9<-array()
  
  N_True_Diseased_children_0_9<-array()
  True_Prev_Disease_children_0_9<-array()
  
  True_Prev_Disease<-array()
  True_Prev_Infection<-array()
  
  #t<-array() #
  
  MDA_round<-array()
  
  for (i in timesim){
    
    if (i %in% MDA_times){
      
      MDA_round[i]<-which(MDA_times==i)
      
      vals<- MDA_timestep(vals,params,demog,MDA_round=MDA_round[i], Tx_mat)
      
      #t[i]<-i  # This no longer used but keep, was used when some MDA rounds targeted children only
      
      #N_Children_ages_1_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>=52))
      
      #N_True_Infected_children_1_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>=52)])
      #True_Prev_Infection_children_1_9[i]<-N_True_Infected_children_1_9[i]/N_Children_ages_1_9[i]
      
      #N_True_Diseased_children_1_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>=52)])
      #True_Prev_Disease_children_1_9[i]<-N_True_Diseased_children_1_9[i]/N_Children_ages_1_9[i]
      
      #True_Prev_Disease[i]<-sum(vals[[2]])/N
      #True_Prev_Infection[i]<-sum(vals[[1]])/N
      
    }
    
    if (i %in% Topical_times){
      
      vals<- Topical_timestep(vals,params,demog)
      
      #t[i]<-i  # This no longer used but keep, was used when some MDA rounds targeted children only
      
      #N_Children_ages_1_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>=52))
      
      #N_True_Infected_children_1_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>=52)])
      #True_Prev_Infection_children_1_9[i]<-N_True_Infected_children_1_9[i]/N_Children_ages_1_9[i]
      
      #N_True_Diseased_children_1_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>=52)])
      #True_Prev_Disease_children_1_9[i]<-N_True_Diseased_children_1_9[i]/N_Children_ages_1_9[i]
      
      #True_Prev_Disease[i]<-sum(vals[[2]])/N
      #True_Prev_Infection[i]<-sum(vals[[1]])/N
      
    }
    
    
    #{ 
    
    vals<- stepF_fixed(vals,params,demog,bet)
    
    ##t[i]<-i-1 #not needed currently, but keep, was used for when simulating children only
    
    N_Children_ages_0_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>=0))
    
    N_True_Infected_children_0_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>=0)])
    True_Prev_Infection_children_0_9[i]<-N_True_Infected_children_0_9[i]/N_Children_ages_0_9[i]
    
    N_True_Diseased_children_0_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>=0)])
    True_Prev_Disease_children_0_9[i]<-N_True_Diseased_children_0_9[i]/N_Children_ages_0_9[i]      
    
    True_Prev_Disease[i]<-sum(vals[[2]])/N
    True_Prev_Infection[i]<-sum(vals[[1]])/N
    
    #}
    
  }
  
  output<-(list(t,True_Prev_Infection_children_0_9,True_Prev_Disease_children_0_9,True_Prev_Disease,True_Prev_Infection))
  names(output)<-c("Time","True_Prev_Infection_children_0_9","True_Prev_Disease_children_0_9","True_Prev_Disease","True_Prev_Infection")
  return(output)
  
  
}

# Set topical treatment dates here then use MDA_dates_as_timestep function to convert to timestep as with MDA
Topical_dates<-as.Date(c("2023-01-01","2024-07-01","2025-01-01"))
Topical_times<-MDA_dates_as_timestep(MDA_dates = Topical_dates, Start_date = Start_date, burnin=burnin)

cores=detectCores()
cl<-makeCluster(cores[1]-6, type="FORK") #
data_store_all_sim<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()
#Parallel runs#
registerDoParallel(cl)
data_store_all_sim<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  data_store_all_sim[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:sim_params$timesim,demog=demog,bet=bet[i],MDA_times=MDA_times,Topical_times=Topical_times,Tx_mat=Tx_mat)
}

stopCluster(cl)

end_time<-Sys.time()
end_time-start_time

Summary_output_df<-extract_output_alltimes(data_store_all_sim = data_store_all_sim, sim_params = sim_params) #

library(ggplot2) 
Mean_inf_plot_children<- ggplot()+
  geom_line(aes(x=Summary_output_df$Time/52,y=Summary_output_df$Mean_inf_children),colour="red")+
  xlab("Time (years)")+
  ylab("Mean Infection Prev Children 0-9")
Mean_inf_plot_children

Mean_dz_plot_children<- ggplot()+
  geom_line(aes(x=Summary_output_df$Time/52,y=Summary_output_df$Mean_dz_children),colour="blue")+
  xlab("Time (years)")+
  ylab("Mean Disease Prev Children 0-9")
Mean_dz_plot_children


