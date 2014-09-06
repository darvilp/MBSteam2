## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

#For SIR all parms will be 0 except for contact rate and recovery rate
pars <- c(	
		contact	= .5,  	# contact rate
		recov 	= .1,   # recovery rate
		ldeath  = .0,
		hdeath	= .0,
		vac		= .0,
		aq		= .0,
		birth	= .0,
		civ		=.00)	

#	Then, the initial values of each state variable:

init.values <- c(S =.99,A=0,V=0,I =.01, R = 0,C=.00,H=0,L=0)

#	The times we want to see

times <- seq(0, 200, by = 1)

#	Now we can define the differential equation model:
SAVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				dS.dt =	birth*S-aq*S-vac*S-contact*I*S-ldeath*S
				dA.dt = aq*S-C*A-hdeath*A+birth
				dV.dt = vac*S-ldeath*V
				dI.dt = contact*I*S-recov*I
				dR.dt = recov*I-hdeath*R
				dC.dt = civ
				dH.dt = hdeath*(R+I+A)
				dL.dt = ldeath*(S+V)
				
				return(list(c(dS.dt,dA.dt,dV.dt, dI.dt, dR.dt,dC.dt,dH.dt,dL.dt)))
			})
}

#	Having defined everything, now we ask the program ode
#	to actually solve the system:

out <- as.data.frame(ode(func = SAVIR, y = init.values, 
				parms = pars, times = times))

#	And now we can make a nice plot of our results:
matplot(out$time, out[ ,c(2,5,6)], type = "l", xlab = "time", 
		ylab = "percent of population", main = "SIR Model", lwd = 2)

legend("topright", c("Susceptible", "Infectious", "Removed"),
		col = 1:3, lty = 1:3)

#matplot(out$time, out[ ,7:8], type = "l", xlab = "time", 
# ylab = "num", main = "Deaths", lwd = 2)

#legend("topright", c("High","Low"),
#		col = 1:2, lty = 1:2)

tail(out)