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

nsim<-50 #Number of simulations
Target_TF_1_9<-20 
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





