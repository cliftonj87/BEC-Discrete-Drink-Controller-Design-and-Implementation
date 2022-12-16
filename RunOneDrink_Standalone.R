# Install the DE solver.  You only need to do this once, so comment out this line
# if you are running this multiple times.
install.packages("deSolve")
library(deSolve) 

# Read in three files.  We are assuming they are in the same folder as this file.
# read in parameters
source('BECParameters.R')
# read in ode file
source('BECode.R')
# read in code to simulate one monkey getting one dose
source('run_one_monkey_onedrink.R')

# Set time values for simulation
tf = 3*60*60  # time is in seconds
t0 = 100
dt = 1
tvals=seq(t0,tf, by=dt) #output every dt seconds
# Set the drink length to be one minute
params[16] = 60
# Call the function that solves the ODE
all_out = run_one_monkey_onedrink(tvals,params)
# extract the simulation output
out = all_out$drink_sim
# Turn the second column of the output into concentration appropriate for BEC
Concentration_V = out[,2]/params[4]/10
# Turn time into minutes
minute_times = out[,1]/60
# Plot
plot(minute_times,Concentration_V,col='purple',type='l',lwd=2,
     main = 'One Drink', ylab='Vascular BEC',xlab='Time in Minutes',ylim=c(0,2))