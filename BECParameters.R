#Parameter Values

Vmax = 3.0265 # Nelder-Mead 4d 04/22/21
V_L = 1.2863 #L volume of the liver # will be  reset later as a function  of monkey's weight
V_T = 13.404 #L volume of the peripheral tissue # will  also  be reset later
V_V = .64 #L volume of the vascular compartment #average for a 9kg male OR reset as a fraction  of weight
R_A = 0.0705 #L/sec volume flow rate to aorta
F_L = 0.26 #unitless fraction of input that exits; fraction of cardiac output directed to the liver
F_P = 0.75 #unitless fraction of liver directed cardiac output that is directed through the gut and portal vein
K_m = 16.55418 # Nelder-Mead-4d then Newton 6/19/2020
AllRate =  0.7931544 # Nelder-Mead 4d then Newton  6/19/2020
k_VT = AllRate #unitless fraction of the available etoh that is actually transported from peripheral artery to tissue
k_TV = AllRate #unitless fraction of the available etoh that is actually transported from tissue to peripheral artery
k_HL = AllRate #unitless fraction of the avilable etoh that is actually transported from hepatic artery to liver
k_PL = AllRate #unitless fraction of the avilable etoh that is actually transported from portal vein to the liver
k_LH = AllRate #unitless fraction of the avilable etoh that is actually transported from liver to hepatic vein
beta =  0.0007363716 # Nelder-Mead-4d then Newton 6/19/2020
dose = 3000 # mg total amount: take g/kg and multiply by weight
dose_length = 30 #  seconds length of dose administration -  30 seconds for gavage of 0.5g/kg of 4% ethanol and 60-90 seconds for 1.0 g/kg of 4% ethanol

# M_gut_multiplier parameters: shows how food slows down rate of alcohol going to the blood
# M_gut_min and M_gut_s were fit to INDUCTION data, Nov/Dec 2020
M_gut_min = .75   # fraction of maximum absorption rate when stomach is completely full ##FIT TO INDUCTION
M_gut_s = 1.5 # controls the steepness of the M_gut_multiplier function: mass of food (in grams) where absorption fraction is halfway to minimum
dig_rate = -log(0.1)/(2*3600) # digestion rate - per second.  Food is 90% digested in 2 hours, info from Kathy
M_gut_par_1 = 0 # this parameter and the next keep track of previous drink's rate and end time
M_gut_par_2 = 0

##Parameters
par1 = Vmax
par2 = V_L
par3 = V_T
par4 = V_V
par5 = R_A
par6 = F_L
par7 = F_P
par8 = K_m
par9 = k_VT
par10 = k_TV
par11 = k_HL
par12 = k_PL
par13 = k_LH
par14 =beta
par15 = dose
par16 = dose_length
par22 = M_gut_min
par23 = M_gut_s
par24 = dig_rate
par25 = M_gut_par_1  
par26 = M_gut_par_2


## Volume calculations based on Excel spreadsheet - December 12, 2018 - used averages
liver_to_body = .0148 # fraction  of body weight that is liver
liver_volume_to_weight = 1/1.13 # L liver volume per kilogram liver

fraction_blood = .064 # L/kg weight from biologist
fraction_bone =  .12 # from Tabble 1, Schroeter 2011
par17 = liver_to_body 
par18 = liver_volume_to_weight
par19 = fraction_blood # fraction of weight  that is blood
par20 = fraction_bone # fraction of weight that is bone
par21 = 1.01 # kg/L tissue  to  volume conversion

params = c(par1,par2,par3,par4,par5,par6,par7,par8,par9,par10,par11,par12,par13,par14,
           par15,par16,par17,par18,par19,par20,par21,par22,par23,par24,par25,par26)
param.names = c('Vmax','V_L','V_T','V_V','R_A','F_L','F_P','K_m','k_VT','k_TV','k_HL','k_PL','k_LH','beta',
                'dose','dose_length','liver_to_body','liver_volume_to_weight','fraction_blood',
                'fraction_bone','tissue_to_vol','M_gut_min','M_gut_s','dig_rate')
