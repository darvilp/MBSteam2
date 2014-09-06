## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

pars <- c(	
		contact	= .9,  	
		recov 	= .01,   
		ldeath  = .00,
		hdeath	= .00,
		vac		= .01,
		aq		= .03,
		birth	= .00,#This can be lower if Civ vaccination is removed
		civ		=.00)	
init.values <- c(S =.5,A=.49,V=.005,I =.005, R = 0,C=.0,H=0,L=0)
times <- seq(0, 1000, by = 1)
SAVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				dS.dt =	birth*(S)-aq*S-(vac)*S-contact*I*S-ldeath*S#-C/75*S
				dA.dt = aq*S-C*A-hdeath*A
				dV.dt = (vac)*S-ldeath*V#+C/75*S
				dI.dt = contact*I*S-hdeath*I-recov*I
				dR.dt = recov*I-hdeath*R
				dC.dt = 0
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
x11()
matplot(out$time, out[ ,2:6], type = "l", xlab = "time", 
		ylab = "percent of population", main = "SAVIR Model with Birth/Death", lwd = 2)

legend("topright", c("Susceptible","Aqcuired Immunity","Vacinnated", "Infectious", "Removed"),
		col = 1:5, lty = 1:5)
#x11()
#matplot(out$time, out[ ,8:9], type = "l", xlab = "time", 
#		ylab = "num", main = "Deaths", lwd = 2)
#
#legend("top", c("High","Low"),
#		col = 1:2, lty = 1:2)

tail(out)