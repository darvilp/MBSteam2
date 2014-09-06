####################################################
#This is a basic skeleton of a script from what
#I remembered of the compartment diagram. 
####################################################
library(deSolve)
rm(list=ls())

#Parameters
pars <- c( 
	vecdeath=.01,
	bite=.01,
	vecexposetime=.01,
	
	birth=.01,
	sdeath=.01,
	edeath=.01,
	exposetime= rnorm(1,mean=.01), #rnorm(1,mean=whatever)in the future. Also could put it in the function
	ideath=.01,
	recovery=rnorm(1,mean=.01), #rnorm(1,mean=whatever)
	rdeath=.01,
	reinfect=.01
			)
#Initial Values
init.values <- c(
		Sv=1,
		Ev=0,
		Iv=0,
		S = 100, E =0,
		I = 0, R = 0
		)

#Times for Graphs
times <- seq(0, 15*365, by = 1)

SEIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				dSv.dt = vecdeath*(Ev+Iv)-Sv*bite*I*R
				dEv.dt = -vecdeath*Ev    +Sv*bite*I*R -vecexposetime*Ev
				dIv.dt = -vecdeath*Iv		        +vecexposetime*Ev
				
				dS.dt =  birth*(S+E+I+R)-sdeath*S -S*bite*Iv 
				dE.dt = 		      -edeath*E +S*bite*Iv -exposetime*E             +reinfect*bite*Iv*R	
				dI.dt = 			-ideath*I    	   +exposetime*E -recovery*I
				dR.dt = 			-rdeath*R				     +recovery*I +reinfect*bite*Iv*R

				return(list(c(dSv.dt,dEv.dt,dIv.dt,
				dS.dt,dE.dt,dI.dt,dR.dt)))
				

			})
}


out = as.data.frame(ode(func = SEIR, y = init.values, 
				parms = pars, times = times))

#---------------------------------------------------
#	Plots
x11()
par(mfrow = c(2,1))

matplot(out$time/365, out[c(2:4)], type = "l", xlab = "Time (years)", 
		ylab = "Percent of Original Population", main = "Mosquitos", lwd = 3)

legend("topright", c("Susceptible","Exposed", 
				"Infectious"),
		col = 1:5, lty = 5)

#People Plot
matplot(out$time/365, out[c(2:4)], type = "l", xlab = "Time (years)", 
		ylab = "Number", main = "Peoples", lwd = 3)

legend("topright", c("Susceptible","Exposed", 
				"Infectious", "Recovered"),
		col = 1:4, lty = 5)
head(out) # last six values of all states
tail(out)

