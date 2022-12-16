##Solve ode "food" model for one drink for fitting gavage data
## April 9, 2020
## We just need to solve for one drink, but this requires two runs of the solver:
##  one run for during the drink, and one run for after the drink.
run_one_monkey_onedrink <- function(times,this_params){
  drink_sim = c()  #  initialize the matrix for storing output
  ## Set parameters for the Mgut function: parameter[25] is the drink number, parameter[26] says whether in a drink or after a drink
  this_params[25] = 1  # always one drink
  this_params[26]=  1 # start with "in drink" and then run again with 0 after end of drink.
  rate = this_params[15]/this_params[16]
  drink_data = list(rate=rate,start = 0,length = this_params[16])
  #solve DE for during drink input.
    t0 = 0
    dt = .1
    t_doseend = this_params[16]  # drink length
    this_times=seq(t0,t_doseend, by=dt)  
    yini=c(mu_V = 0, mu_T = 0, mu_L = 0,Food = 0, Mgut_Fun = 0) 
    sim=lsoda(y = yini, times = this_times, func = BEC_RHS, parms = list(params_ode=c(this_params,t0),drink_data=drink_data), jactype="fullint")
    # save simulated values
    drink_sim = sim
  # solve DE for time after drink
    this_params[15]=0 # volume of drink = 0
    this_times = c(t_doseend,times)
    # tell solver that we are not in a drink
    this_params[26] = 0
    yini = sim[length(sim[,1]),2:6]
    sim=lsoda(y = yini, times = this_times, func = BEC_RHS, parms = list(params_ode=c(this_params,t0),drink_data=drink_data), jactype="fullint")
    L=length(drink_sim[,1])
    drink_sim = drink_sim[-L,]   # get rid of last line of previous simulation so as not to count that time point twice
    drink_sim = rbind(drink_sim,sim)
    sim_output <- list(drink_sim = drink_sim,sim = sim)
return(sim_output)
}

