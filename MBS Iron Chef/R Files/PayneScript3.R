####################################################
#This is based on the parasitic compartment model 
#I have one arrow different because I think it was copied wrong on mistake
#So Ruptured Blood Cell goes to T instead of Susceptible Blood Cell
####################################################
library(deSolve)
rm(list=ls())

#Parameters
pars <- c( 
	#mosquito stuff
	vecdeath=.13,
	bites=.022,
	bitel=.048,
	biteh=.48,
	vecexposetime=.091,
	#birth and death
	birth=0,
	sdeath=0,
	edeath=0,
	ildeath=.000009,
	ihdeath=.000009,
	#class transfers
	etime=5.4,
	iltime=.01,	
	ihtoil=.0035
	
	#exposetime= rnorm(1,mean=.01), #rnorm(1,mean=whatever)in the future. Also could put it in the function
	#recovery=rnorm(1,mean=.01), #rnorm(1,mean=whatever)
	
			)
#Initial Values
init.values <- c(
		Sv=1,
		Ev=0,
		Iv=0,
		
		S = 100, E =0,
		IL = 15, IH=0,
		
		D=0
		)

#Times for Graphs
times <- seq(0, 1*365, by = 1)

SEIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				dSv.dt =  vecdeath*(Ev+Iv)-Sv*bitel*IL-Sv*biteh*IH
				dEv.dt = -vecdeath*Ev     +Sv*bitel*IL+Sv*biteh*IH -vecexposetime*Ev
				dIv.dt = -vecdeath*Iv		        			             +vecexposetime*Ev
				
				#        Deaths    Moving to other states		birth/other
				dS.dt = -sdeath*S						          -bites*S*Iv		+birth*(S+E+IL+IH)
				dE.dt = -edeath*E 	+bites*S*Iv 	 	 -etime*E
				dIL.dt= -ildeath*IL	+etime*E+ihtoil*IH  -iltime*IL
				dIH.dt= -ihdeath*IH	+iltime*IL			 -ihtoil*IH			
				dD.dt =  ildeath*IL	+ihdeath*IH
				
				return(list(c(dSv.dt,dEv.dt,dIv.dt,
				dS.dt,dE.dt,dIL.dt,dIH.dt,dD.dt)))
				

			})
}


out = as.data.frame(ode(func = SEIR, y = init.values, 
				parms = pars, times = times))

#---------------------------------------------------
#	Plots
x11()
par(mfrow = c(2,1))

matplot(out$time, out[c(2:4)], type = "l", xlab = "Time (days)", 
		ylab = "Percent of Original Population", main = "Mosquitos", lwd = 3)

legend("topright", c("Susceptible","Exposed", 
				"Infectious"),
		col = 1:5, lty = 5)

#People Plot
matplot(out$time, out[c(5:8)], type = "l", xlab = "Time (days)", 
		ylab = "Number", main = "People", lwd = 3)

legend("topright", c("S","E","IL","IH"),
		col = 1:4, lty = 5)
head(out) # last six values of all states
tail(out)

