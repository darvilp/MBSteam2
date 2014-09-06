## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

pars <- c(	
		contact	= .9,  	
		recov 	= .1,   
		ldeath  = .1,
		hdeath	= .2,
		vac		= .07,
		aq		= .3,
		birth	= .5,
		civ		=.01)	
init.values <- c(S =.5,A=.49,V=.005,I =.005, R = 0,C=.001,H=0,L=0)
times <- seq(0, 300, by = 1)
SAVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				dS.dt =	birth*S-aq*S-vac*S-contact*I*S-ldeath*S
				dA.dt = aq*S-C*A-hdeath*A
				dV.dt = vac*S-ldeath*V
				dI.dt = contact*I*S-hdeath*I-recov*I
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
x11()
matplot(out$time, out[ ,2:6], type = "l", xlab = "time", 
		ylab = "percent of population", main = "SAVIR Model with Birth/Death", lwd = 2)

legend("topright", c("Susceptible","Aqcuired Immunity","Vacinnated", "Infectious", "Removed"),
		col = 1:5, lty = 1:5)
#x11()
#matplot(out$time, out[ ,7:8], type = "l", xlab = "time", 
# ylab = "num", main = "Deaths", lwd = 2)

#legend("topright", c("High","Low"),
#		col = 1:2, lty = 1:2)

tail(out)