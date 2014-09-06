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
	vecdeath=.01,
	bite=.01,
	vecexposetime=.01,
	#birth and death
	birth=.01,
	sdeath=.01,
	edeath=.01,
	ildeath=.01,
	rldeath=.01,
	sbdeath=.01,
	tdeath=.01,
	scdeath=.01,
	rsdeath=.01,
	gdeath=.01,
	#class transfers
	etime=.01,
	rltime=.01,
	sbtime=.01,
	sb2t=.01,
	rs2t=.01,
	t2sc=.01,
	t2g=.01,
	rstime=.01
	
	
	#exposetime= rnorm(1,mean=.01), #rnorm(1,mean=whatever)in the future. Also could put it in the function
	#recovery=rnorm(1,mean=.01), #rnorm(1,mean=whatever)
	
			)
#Initial Values
init.values <- c(
		Sv=1,
		Ev=0,
		Iv=0,
		
		S = 100, E =0,
		IL = 15, RL=0,
		
		SB=0,T=0,SC=0,RS=0,G=0
		)

#Times for Graphs
times <- seq(0, 15*365, by = 1)

SEIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				dSv.dt = vecdeath*(Ev+Iv)-Sv*bite*G
				dEv.dt = -vecdeath*Ev    +Sv*bite*G -vecexposetime*Ev
				dIv.dt = -vecdeath*Iv		        +vecexposetime*Ev
				
				#        Deaths    Moving to other states		birth/other
				dS.dt = -sdeath*S						-bite*S*Iv	+birth*(S+E+IL+RL+SB+T+SC+RS+G)
				dE.dt = -edeath*E 	+bite*S*Iv 	 		-etime*E
				dIL.dt= -ildeath*IL	+etime*E			-rltime*IL
				dRL.dt= -rldeath*RL	+rltime*IL			-sbtime*RL
				dSB.dt= -sbdeath*SB	+sbtime*RL 			-sb2t*SB
				dT.dt = -tdeath*T	+sb2t*SB +rs2t*RS	-t2sc*T-t2g*T
				dSC.dt= -scdeath*SC	+t2sc*T				-rstime*SC
				dRS.dt= -rsdeath*RS +rstime*SC			-rs2t*RS
				dG.dt = -gdeath*G	+t2sc*T

				return(list(c(dSv.dt,dEv.dt,dIv.dt,
				dS.dt,dE.dt,dIL.dt,dRL.dt,dSB.dt,dT.dt,dSC.dt,dRS.dt,dG.dt)))
				

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
matplot(out$time/365, out[c(5:13)], type = "l", xlab = "Time (years)", 
		ylab = "Number", main = "People", lwd = 3)

legend("topright", c("Garbage"),
		col = 1:4, lty = 5)
head(out) # last six values of all states
tail(out)

