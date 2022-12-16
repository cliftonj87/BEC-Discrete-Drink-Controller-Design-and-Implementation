BEC_RHS <- function(t,y, params){
  with(as.list(y),{
  # params is a list.  
  # The first element is a vector of parameter values,
  # the second element is a matrix with drink data used in the Mgut DE.
  params_ode=params$params_ode
  drink_data = params$drink_data
  mu_V = y[1]
  mu_T = y[2]
  mu_L = y[3]
  Food = y[4]
  Mgut_Val = y[5]
  
  Vmax = params_ode[1]
  V_L = params_ode[2]
  V_T = params_ode[3]
  V_V = params_ode[4]
  R_A = params_ode[5]
  F_L = params_ode[6]
  F_P = params_ode[7]
  K_m = params_ode[8]
  k_VT = params_ode[9]
  k_TV = params_ode[10]
  k_HL = params_ode[11]
  k_PL = params_ode[12]
  k_LH = params_ode[13]
  beta  = params_ode[14]
  dose = params_ode[15]  ## use for one drink- dose in mg
  ## use next line if dose in ml at 4%
  #dose = params_ode[15]*.04*1000  # .04 g/mL*1000 mg/g -> want mg, dose given in ml of 4% ethanol
  dose_length = params_ode[16]

  # M_gut_multiplier parameters: shows how food slows down rate of alcohol going to the blood
  M_gut_min = params_ode[22]
  M_gut_s = params_ode[23]  # controls the steepness of the M_gut_multiplier function
  dig_rate = params_ode[24] # digestion rate
  drink_number = params_ode[25]  # this is the number of the drink, used in M_gut_function
  drinking = params_ode[26]  # this is equal to 1 if we are in the middle of a drink, used in M_gut_function
  t0 = drink_data$start[drink_number]
  r = dose/dose_length  # mg/sec - rate
  
  C_V = mu_V/V_V  # concentration in vasculature
  C_L = mu_L/V_L # concentration in liver
  C_T = mu_T/V_T # concentration in periheral tissue
  
  # rate equations: functions of time
  M_gut_multiplier = ((1-M_gut_min)/pi)*(-atan(Food-M_gut_s)+pi/2)+M_gut_min
  M_gut = M_gut_multiplier*Mgut_Val
  M_VT = k_VT*R_A*(1-F_L)*((C_V-C_T)>0)*(C_V-C_T)
  M_TV = k_TV*R_A*(1-F_L)*((C_T-C_V)>0)*(C_T-C_V)
  M_HL = k_HL*R_A*F_L*(1-F_P)*((C_V-C_L)>0)*(C_V-C_L)
  gamma = C_L  # to be consistent with Plawecki
  alpha = (1-F_P)*C_V - k_HL*(1-F_P)*((C_V-C_L)>0)*(C_V-C_L)
  + C_V*F_P + (1-k_PL)*M_gut -  k_PL*F_P*((C_V-C_L)>0)*(C_V-C_L)
  if (alpha < gamma){
    C_HV = 1/(1+k_LH)*(alpha + k_LH*gamma)
  } else {C_HV = alpha}  # concentration in hepatic vein
  M_LH = k_LH*R_A*F_L*((C_L - C_HV)>0)*(C_L-C_HV)
  M_PL = k_PL*R_A*F_L*F_P*((C_V + M_gut/(R_A*F_L*F_P) - C_L)>0)*(C_V + M_gut/(R_A*F_L*F_P) - C_L)
  M_metab = V_L*Vmax*C_L/(K_m + C_L)
  
  # differential equations: rate in - rate out
  dmu_Vdt = -M_HL - M_PL + M_LH + M_TV - M_VT + M_gut
  dmu_Ldt = M_PL + M_HL - M_LH - M_metab
  dmu_Tdt = M_VT - M_TV
  dFood_dt = - dig_rate*Food
 # during drink
  if (drinking ==1){
    # contribution from previous drinks
    if (drink_number == 1) {
      prev_drinks = 0
    }
    else{  
  prev_drinks = sum(drink_data$rate[1:(drink_number-1)]*beta*exp(-beta*t)*(exp(beta*(drink_data$start[1:(drink_number-1)]+drink_data$length[1:(drink_number-1)]))-exp(beta*drink_data$start[1:(drink_number-1)]))) 
   }
  dMgut_Val = r*beta*(1-exp(-beta*(t-t0))) - beta*Mgut_Val + prev_drinks
  }  # end if drinking ==1
  else { 
    prev_drinks = sum(drink_data$rate[1:(drink_number)]*beta*exp(-beta*t)*(exp(beta*(drink_data$start[1:(drink_number)]+drink_data$length[1:(drink_number)]))-exp(beta*drink_data$start[1:(drink_number)])))
    dMgut_Val =  - beta*Mgut_Val + prev_drinks
  }
  list(c(dmu_Vdt, dmu_Tdt, dmu_Ldt,dFood_dt,dMgut_Val))
    })
  }

