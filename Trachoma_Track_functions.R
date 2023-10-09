###### UPDATED TRACHOMA FUNCTIONS as of October 2023 ####
#UPDATES= 
#New age immunity function with duration of D period a function of age and history of infection. 
#New age immunity function with duration ID period a function of age and history of infection
# Change to age group 0-9 to better fit with available Tanzania data (normally ages 1-9)
# Reset individuals now set bacterial load to zero
#Taken out scaling parameter gamma from get lambda function
#Output is true prevalence of infection/disease at each timestep in whole population (True_prev_Infection/True_Prev_Disease) or children aged 0-9 (True_Prev_disease_children_0_9/True_Prev_infection_children_0_9)
#This version originally called "Trachoma_functions_ID_09_redo"

sim_Ind_MDA<-function(params,vals,timesim,demog,bet, posterior_values, MDA_times, Tx_mat){
  
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
    
    
    #{ 
    
    vals<- stepF_fixed(vals,params,demog,bet)
    
    ##t[i]<-i-1 #not needed currently, but keep, was used for when simulating MDA for children only
    
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

stepF_fixed <- function(vals,params,demog,bet){
  #Assign names to variables just so it's easier to track # 
  IndI<-vals$IndI
  IndD<-vals$IndD
  No_Inf<-vals$No_Inf
  T_latent<-vals$T_latent
  T_ID<-vals$T_ID
  T_D<-vals$T_D
  Ind_latent<-vals$Ind_latent
  Ind_ID_period_base<-vals$Ind_ID_period_base
  Ind_D_period_base<-vals$Ind_D_period_base
  bact_load<-vals$bact_load
  Age<-vals$Age
  
  # Step 1: Identify individuals available for infection.
  Ss <- which(IndI==0 & IndD==0) 
  Ds <- which(IndI==0 & IndD==1 & T_latent==0) 
  
  # Step 2: Calculate infection pressure from previous time step and choose infected individuals  
  # This give a lambda for each individual dependent on age and disease status
  lambda_step<-1-exp(-getlambdaStep(params=params, Age=Age, demog=demog, bact_load=bact_load, 
                                    IndD=IndD, bet=bet))
  # new infections
  NewInf<-Ss[which(runif(length(Ss)) < lambda_step[Ss])] 
  NewReInf<-Ds[which(runif(length(Ds)) < lambda_step[Ds])]
  
  # Step 3: Identify transitions
  newDis <-which(T_latent==1) #Designated latent period for that individual is about to expire
  newClearInf <- which(T_ID==1) # Designated infectious period for that individual is about to expire
  newClearDis <- which(T_D==1) # Designated diseased period for that individual is about to expire
  newInfectious <- which(IndI==1 & T_latent==1) # Only individuals who have avoided MDA become infectious at end of latent period
  
  # Step 4: reduce counters
  T_latent[which(T_latent>0)]<-T_latent[which(T_latent>0)]-1 
  T_ID[which(T_ID>0)]<-T_ID[which(T_ID>0)]-1
  T_D[which(T_D>0)]<-T_D[which(T_D>0)]-1
  
  # Step 5: implement transitions
  # Transition: become diseased 
  IndD[newDis]<-1
  T_ID[newDis]<-ID_period_function(Ind_ID_period_base=Ind_ID_period_base[newDis],No_Inf = 
                                     No_Inf[newDis],params=params,Age=Age[newDis])
  # Transition: clear infection
  IndI[newClearInf]<-0
  T_D[newClearInf]<-D_period_function(Ind_D_period_base=Ind_D_period_base[newClearInf],No_Inf = 
                                        No_Inf[newClearInf],params=params,Age=Age[newClearInf])
  # Stop being infectious too
  bact_load[newClearInf]<-0
  # Transition: clear disease
  IndD[newClearDis]<-0
  # Transition: become infectious
  bact_load[newInfectious]<-bacterialLoad(No_Inf = No_Inf[newInfectious],params=params)
  
  # Step 6: implement infections
  # Transition: become infected
  IndI[NewInf]<-1
  IndI[NewReInf]<-1
  No_Inf[NewInf]<-No_Inf[NewInf]+1 
  No_Inf[NewReInf]<-No_Inf[NewReInf]+1
  T_latent[NewInf]<-Ind_latent[NewInf]
  T_latent[NewReInf]<-Ind_latent[NewReInf]
  T_ID[NewReInf]<-0 # prevent reinfections from recovering prematurely.
  T_D[NewReInf]<-0 # prevent reinfections from recovering prematurely.
  
  # Step 7: update ages
  #Update age, all age by 1w at each timestep, and reseting all "reset indivs" age to zero
  #Reset_indivs - Identify individuals who die in this timeset, either reach max age or random death rate
  Age<-Age+1
  reset_indivs<-Reset(Age=Age,demog=demog,params=params)
  Age[reset_indivs]<-0 
  #Resetting new params for all new individuals created 
  IndI[reset_indivs]<-0
  IndD[reset_indivs]<-0
  No_Inf[reset_indivs]<-0 
  T_latent[reset_indivs]<-0 
  T_ID[reset_indivs]<-0
  T_D[reset_indivs]<-0
  bact_load[reset_indivs]<-0 ######****
  
  vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
  names(vals)<-c("IndI","IndD","No_Inf","T_latent","T_ID","T_D","Ind_latent","Ind_ID_period_base",
                 "Ind_D_period_base","bact_load","Age")
  
  return(c(vals))
  
}


findYoungChildren = function(Age, youngChildMaxAge){
  return(which(Age< to_weeks(youngChildMaxAge)))
}

findAverageBactLoadYoungChildren <- function(bact_load, y_children){
  return(sum(bact_load[y_children])/length(y_children))
}


findBabies <- function(Age, babiesMaxAge){
  return(which(Age < to_weeks(babiesMaxAge)))
}


findOlderChildren = function(Age, youngChildMaxAge, olderChildMaxAge){
  return(which(Age>=to_weeks(youngChildMaxAge) & Age<to_weeks(olderChildMaxAge)))
}

findAverageBactLoadOlderChildren <- function(bact_load, o_children){
  return(sum(bact_load[o_children])/length(o_children))
}


findAdults = function(Age, olderChildMaxAge){
  return(which(Age>=to_weeks(olderChildMaxAge)))
}

findAverageBactLoadAdults <- function(bact_load, adults){
  return(sum(bact_load[adults])/length(adults))
}

#lambda scaled according to disease status
getlambdaStep = function(params,Age,demog,bact_load,IndD,bet){
  
  totalLoad <- c(0,0,0)
  y_children <- findYoungChildren(Age, params$youngChildMaxAge) #Young children
  totalLoad[1] <- sum(bact_load[y_children])/length(y_children)
  o_children <- findOlderChildren(Age, params$youngChildMaxAge, params$olderChildMaxAge)
  totalLoad[2] <- sum(bact_load[o_children])/length(o_children)
  adults <- findAdults(Age, params$olderChildMaxAge) #Adults
  totalLoad[3] <- sum(bact_load[adults])/length(adults)
  
  prevLambda <- bet * (params$v_1 * totalLoad + params$v_2 
                       * totalLoad ^ (params$phi + 1))
  
  demog_matrix <- t(matrix(rep(c(length(y_children),length(o_children),length(adults))/params$N,3), nrow=3, ncol=3))
  social_mixing <- (params$epsilon*diag(3) + (1 - params$epsilon))*demog_matrix #scales mixing with other groups 
  
  lambda <- as.vector(social_mixing%*%prevLambda) #note %*% is matrix multiplication
  
  return((lambda[findInterval(Age,c(0, to_weeks(params$youngChildMaxAge), to_weeks(params$olderChildMaxAge),
                                    to_weeks(demog$max_age)))])) 
  
  
}

#Function to identify individuals who either die due to background mortality, or who reach max age. 
Reset<-function(Age,demog,params){
  
  reset_indivs<-(which(runif(params$N)<1-exp(-demog$tau)| Age>demog$max_age))  
  return(reset_indivs)
  
}


#Function to decide when MDA occurs, return vector of timepoints when MDA carried out. 
#Can be annually or less freq/more freq

Set_t_MDA=function(sim_params,Freq_MDA){
  
  MDA_t<-c()
  for (i in 1: sim_params$N_MDA){
    MDA_t[i]<-sim_params$burnin+((i-1)*52*Freq_MDA) + 52 #Added 52 for if want first MDA to be 1 year after burnin, delete if doesn't matter
  }
  return(MDA_t)
}

#Create matrix to determine who gets treated at each MDA round, allowing for systematic non-compliance as specified by Dyson
Tx_matrix <- function(params,sim_params){
  
  N = params$N
  N_MDA = sim_params$N_MDA
  
  set.seed(5)
  ind_treat <- matrix(0,N,N_MDA)
  for(i in 1: N){
    z <- runif(N)
    ind_treat[,1]<-z < params$MDA_Cov #Randomly assign first treatment
    
    for (k in 2:N_MDA){
      ind_treat[i,k]<-rbinom(1,1,((params$MDA_Cov*(1-params$rho)+
                                     (params$rho*sum(ind_treat[i,1:k])))/(1+(k-2)*params$rho))) #Subsequent treatment probs function of previous treatments
      
    }
  }
  return(ind_treat)
  
}


#Decide who is cured during MDA- based on treatment matrix and probability of clearance given treated
doMDA<-function(params,Age,MDA_round, Tx_mat){
  
  babies<-which(Age<=26)
  MDA_round<-MDA_round
  treated_babies<-babies[which(Tx_mat[babies,MDA_round]==1)]
  cured_babies<-treated_babies[which(runif(treated_babies)<(params$MDA_Eff*0.5))]
  older<-which(Age>26)
  treated_older<-older[which(Tx_mat[older,MDA_round]==1)]
  cured_older<-treated_older[which(runif(treated_older)< params$MDA_Eff)]
  treated_cured<-c(cured_babies,cured_older)
  return(treated_cured)
  
}



#This is time step in which MDA occurs#

MDA_timestep <- function(vals,params,demog,MDA_round, Tx_mat){
  
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
  treated_cured<-doMDA(params=params,Age=Age,MDA_round=MDA_round, Tx_mat)
  
  #Set treated/cured indivs infection status and bacterial load to 0
  IndI[treated_cured]<-0 #clear infection they become I=0
  bact_load[treated_cured]<-0 #clear disease they become D=0
  
  
  vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
  names(vals)<-c("IndI","IndD","No_Inf","T_latent","T_ID","T_D","Ind_latent","Ind_ID_period_base",
                 "Ind_D_period_base","bact_load","Age")
  
  
  return(vals)
} 


#UPDATED FUNCTION to give duration of active infection

ID_period_function<-function(Ind_ID_period_base,No_Inf,params,Age){
  #T_ID<- round((Ind_ID_period_base-params$min_ID)*exp(-params$inf_red*(No_Inf-1))+params$min_ID)
  #ag<--0.00012 
  #aq<--0.039
  T_ID<- round(params$min_ID +(Ind_ID_period_base-params$min_ID)*exp(params$aq_ID*(No_Inf)+params$ag_ID*(Age)))
  return(T_ID)
}


# UPDATED FUNCTION includes age component of decay in disease duration, see "New_age_immunity" for alternative
D_period_function<-function(Ind_D_period_base,params,No_Inf,Age){
  #ag<--0.00098
  #aq<--0.078  
  T_D<- round(params$min_D +(Ind_D_period_base-params$min_D)*exp(params$aq_D*(No_Inf)+params$ag_D*(Age)))
  return(T_D)
}



#Function to scale bacterial load according to infection history
bacterialLoad = function(No_Inf,params){
  
  b1 = 1     ## params for bacterial load function, can put in general params
  ep2 = 0.114 ## params for bacterial load function
  bact_load=(b1*exp((No_Inf-1)*-ep2))
  return(bact_load) ## functional form for reduction in bacterial load with repeated infections
}


####Don't need to use this if have long enough burn in, but sets history of infection etc to more realistic values when initialise so stabilises faster
#Function to give category based on target TF, specified as a percentage, if simulating v high transmission may need to refine here
Target_TF_Cat<-function(Target_TF_1_9){
  Cat<-ifelse(Target_TF_1_9>=5&Target_TF_1_9<10,1,
              ifelse(Target_TF_1_9>=10&Target_TF_1_9<20,2,
                     ifelse(Target_TF_1_9>=20&Target_TF_1_9<30,3,
                            ifelse(Target_TF_1_9>=30&Target_TF_1_9<40,4,
                                   ifelse(Target_TF_1_9>=40&Target_TF_1_9<100,5,
                                          
                                          NA)))))
  return(Cat)
}

# Set Initial values: Updated inits function, assigns average duration of infection based on age given target TF, assumed simple linear relationship, plan to update so draw from poisson with this average as lambda
# Coefficients for assumed relationship between age and number of infection is stored in "inits" which needs to be loaded when starting simulations
# Note, targe TF is a % not a proportion here 

Set_inits<-function(params,demog, inits, Target_TF_1_9){
  N = params$N
  IndI <- rep(0,N) #Individual's infected status
  IndD <- rep(0,N) #Individual's disease status
  T_latent<-rep(0,N) #Duration of latent period (I), i.e. infected but not yet clinically diseased.
  T_ID<-rep(0,N) #Duration of current ID period,set when becomes infected and counts down with each time step
  T_D<-rep(0,N) #Duration individual spends diseased after clearing infection
  Ind_latent<-rep(params$av_I_duration,N) #Individual's latent period fixed for now
  # Ind_latent<-rpois(N, params$av_I_duration) #Individual's latent period fixed for now
  Ind_ID_period_base<-rpois(N,params$av_ID_duration) #Individual's baseline ID period (first infection)- need to check that this isn't less than the minimum
  Ind_D_period_base<-rpois(N,params$av_D_duration) #Individual's baseline diseased period (first infection)- need to check that this isn't less that the minimum
  bact_load<-rep(0,N) #Baseline bacterial load set to zero
  
  Cat=Target_TF_Cat(Target_TF_1_9)
  
  Age<-init_ages(params=parameters,demog=demog)
  
  No_Inf<-Assign_init_no_inf(Cat=Cat,inits=inits,Age=Age,N=N)
  
  vals <<- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
  names(vals)<<-c("IndI","IndD","No_Inf","T_latent","T_ID","T_D","Ind_latent","Ind_ID_period_base",
                  "Ind_D_period_base","bact_load","Age") 
  
  
  
}


Assign_init_no_inf<-function(Cat,inits,Age,N){
  No_Inf<-rep(0,N)
  for (i in 1:N){
    #No_Inf[i]=round(inits$Mean_Coeff1[Cat]+inits$Mean_Coeff2[Cat]*Age[i])
    No_Inf[i]=rpois(n=1, lambda=round(inits$Mean_Coeff1[Cat]+inits$Mean_Coeff2[Cat]*Age[i]))
  }
  return(No_Inf)
}

Seed_infection<-function(params,vals,Target_TF_1_9,inits){
  Cat=Target_TF_Cat(Target_TF_1_9) #Assign category based on target baseline TF
  Seed_percent=inits$Mean_Infection[Cat] #Assign seed % based on mean eqilibrium infect % for target TF
  vals$IndI[1:round(params$N*Seed_percent)] <<- 1 #set % infected dependent on 
  Init_infected <- which(vals$IndI==1)
  vals$T_latent[Init_infected]<<-vals$Ind_latent[Init_infected]# set latent period for those infected at start of simulation
  vals$No_Inf[Init_infected]<<-vals$No_Inf[Init_infected]+1 # Set number of infections to +1 for those infected at start of simulation.
}



#Initialise age distribution, draw from truncated geometric, needs trundist package 

#NOTE- ages are in weeks
init_ages<-function(params,demog){
  
  Ages<-rtrunc(params$N,spec="geom", prob=1-exp(-demog$tau), a=0, b=demog$max_age)
  
  return(Ages)
  
}

#Function if want to draw beta from prior, can still specify beta separately if want to 
bet_from_prior<-function(Target_TF_1_9,inits,sim_params){
  set.seed(5)
  Cat=Target_TF_Cat(Target_TF_1_9) #Assign category based on target baseline TF
  #bet<-rnorm(sim_params$nsim,mean=inits$Mean_Beta[Cat],sd=inits$S_D_Beta[Cat])
  bet<-rtrunc(sim_params$nsim,spec="norm", mean=inits$Mean_Beta[Cat],sd=inits$S_D_Beta[Cat],a=0.05)
  return(bet)
}

initialise_pop_and_seed<-function(parameters, demog, inits, Target_TF_1_9){
  set.seed(5)
  Set_inits(params = parameters, demog=demog, inits=inits, Target_TF_1_9 = Target_TF_1_9)
  Seed_infection(params = parameters, vals=vals, Target_TF_1_9=Target_TF_1_9, inits=inits)
}

create_sim_params<-function(nsim, burnin, Start_date,End_date,MDA_dates){
  post_burn<-ceiling(as.numeric(difftime(End_date, Start_date,units="weeks")))
  timesim<-burnin+post_burn
  sim_params<-as.list(c(timesim=timesim,burnin=burnin,N_MDA=length(MDA_dates),nsim=nsim))
}

MDA_dates_as_timestep<-function(MDA_dates, Start_date, burnin){
  MDA_times<-as.integer(difftime(MDA_dates, Start_date, units = "weeks") )+burnin 
}

to_weeks <- function(no_of_years){
  return (no_of_years * 52)
}

#Extracting output for plotting
extract_output_alltimes<-function(data_store_all_sim,sim_params){
  nsim<-sim_params$nsim
  timesim<-sim_params$timesim
  True_Prev_Infection_children<-matrix(0,nrow=timesim,ncol=nsim)
  True_Prev_Disease_children<-matrix(0,nrow=timesim,ncol=nsim)
  True_Prev_Infection<-matrix(0,nrow=timesim,ncol=nsim)
  True_Prev_Disease<-matrix(0,nrow=timesim,ncol=nsim)
  Time<-seq(0,timesim-1,1)
  
  #Fill with output, 
  for(i in 1:nsim){
    True_Prev_Infection_children[,i]<-data_store_all_sim[[i]]$True_Prev_Infection_children_0_9[1:timesim]
    True_Prev_Disease_children[,i]<-data_store_all_sim[[i]]$True_Prev_Disease_children_0_9[1:timesim]  
    True_Prev_Infection[,i]<-data_store_all_sim[[i]]$True_Prev_Infection[1:timesim]
    True_Prev_Disease[,i]<-data_store_all_sim[[i]]$True_Prev_Disease[1:timesim]  
    
  }
  
  
  Mean_dz_children<-apply(True_Prev_Disease_children,1,mean)
  Mean_inf_children<-apply(True_Prev_Infection_children,1,mean)
  Mean_dz<-apply(True_Prev_Disease,1,mean)
  Mean_inf<-apply(True_Prev_Infection,1,mean)
  
  
  
  Summary_output_df<-data.frame(Time, Mean_dz_children, Mean_inf_children,Mean_dz,Mean_inf) 
  
  
  return("summary_output_df"=Summary_output_df) # 
  
  
}
