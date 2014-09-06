## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#  First, the parameters for our model:

#NOTE: unit of time is months
# Also not all pars are set. still missing betas and total pop data

pars <- c(m = 1.217 * 10^-3, # immigration rate of humans
          b = 2.417 * 10^-3, # per capita birth rate of humans
          B = 4.227, # per capita birth rate of mosquitos
          L = 3.0438, # rate of human movement from exposed to infectious
          u = 3.04375, # rate of mosquito movement from exposed to infectious
          c = 8.333 * 10^-2, # rate of loss of human immunity
          d = 1 * 10^-7, # rate of disease induced death in humans
          e = 1.2 * 10^-2, # proportion of bites that causes infection (mosquito->human)
          delta1 = 3.623, # density-independent death and emigration for mosquitoes
          delta2 = 6.722 * 10^-7, # density-dependent death and emigration for mosquitoes
          E = 4.7 * 10^-1, # proportion of bites that causes infection (infectious human->mosquito)
          Es = 2.35 * 10^-1, # proportion of bites that causes infection (recovered human->mosquito)
          F = 7.609, # average number of bites per mosquito
          D1= 4.808 * 10^-4, # density-independent death and emigration rate for humans
          D2 = 1.000 * 10^-5, # density-dependent death and emigration rate for humans
          q = 8.333 * 10^2, # rate of building effective immunity
          r = 5.558 * 10^-2)	# per capita rate of recovery

#	Then, the initial values of each state variable:

init.values <- c(Sh = .999, Eh = 0, Ih = .001, Rh = 0, Sm = 0, Em = 0, Im = 0)

#	The times we want to see

times <- seq(0, 100, by = 1)

#	Now we can define the differential equation model:

SEIR <- function(time, y.values, parameters) {
  with(as.list(c(y.values, parameters)), {
    
    dSh.dt = m + b * Nh + c * Rh - betaMH * (Im/Nh) * Sh + r * Ih - (D1 + D2 * Nh) * Sh
    dEh.dt = betaMH * (Im/Nh) * Sh - L * Eh - (D1 + D2 * Nh) * Eh
    dIh.dt = L * Eh - ((q * r)/(q + r)) * Ih - r * Ih - d * Ih - (D1 + D2 * Nh) * Ih
    dRh.dt = ((q * r)/(q + r)) * Ih - c * Rh -(D1 + D2 * Nh) * Rh
    dSm.dt = B * Nm - betaHM * (Ih/Nh) * Sm - betasHM * (Rh/Nh) * Sm - (delta1 + delta2 * Nm) * Sm
    dEm.dt = betaHM * (Ih/Nh) * Sm + betasHM * (Rh/Nh) * Sm - u * Em - (delta1 + delta2 * Nm) * Em
    dIm.dt = u * Em - (delta1 + delta2 * Nm) * Im
    
    return(list(c(dS.dt, dE.dt, dI.dt, dR.dt)))
  })
}

#	Having defined everything, now we ask the program ode
#	to actually solve the system:

out <- as.data.frame(ode(func = SEIR, y = init.values, 
                         parms = pars, times = times))

#	And now we can make a nice plot of our results:

matplot(out$time, out[ ,2:8], type = "l", xlab = "time", 
        ylab = "percent of population", main = "SEIR Malaria Model with Birth/Death", lwd = 2)

legend("topright", c("Sh", "Eh", "Ih", "Rh", "Sm", "Em", "Im"),
       col = 1:7, lty = 1:7)

tail(out)