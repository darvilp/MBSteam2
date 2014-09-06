#This was the emailed copy

## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

#Set of pars from early March
#pars <- c(
#		contacty= 0.190,  	
#		recovy 	= 0.028,   
#		ldeathy  = 6.8/365/1000, #disease independent death
#		hdeathy	= .004633/1000, #disease induced death for young
#		vacy		= .05,
#		aqy		= .03, #Aqcuired immunity rate
#		birthy	= .0, #This can be lower if Civ vaccination is Recovered
#		civy		= .0001, #Modernization factor
#		
#		aging = 1/365,
#		dehyd = .1, #percentage of vaccines that do not work
#		
#		contact	= .142,  	
#		recov 	= .028,   
#		ldeath  = 6.8/365/1000, #disease independent death
#		hdeath	= .00567/1000, #disease induced death for old
#		vac		= .07,
#		aq		= .03,
#		birth	= 24.3/365/1000, #This can be lower if Civ vaccination is Recovered; dS.dt
#		civ		= .0001)	

#pars to produce an SIR model that shows stability eventually with S and I constant and all "births" going into R
#Consider this the current case in pakistan w/o vac
#SIR (novac)
pars <- c(
		contacty= 0.190/1,  	
		recovy 	= 0.028,   
		ldeathy  = 6.8/365/1000, #disease independent death
		hdeathy	= .004633/1000, #disease induced death for young
		vacy		= 0,
		aqy		= .0, #Aqcuired immunity rate
		birthy	= .0, #This can be lower if Civ vaccination is Recovered
		civy		= .000, #Modernization factor
		
		aging = 1/365/14,
		dehyd = 0, #percentage of vaccines that do not work
		
		contact	= .142/1,  	
		recov 	= .028,   
		ldeath  = 6.8/365/1000, #disease independent death
		hdeath	= .00567/100, #disease induced death for old
		vac		= 0,
		aq		= .0,
		birth	= 24.3/365/1000, #This can be lower if Civ vaccination is Recovered; dS.dt
		civ		= .0)
#SVIR w/vac
pars <- c(
		contacty= 0.190/1,  	
		recovy 	= 0.028,   
		ldeathy  = 6.8/365/1000, #disease independent death
		hdeathy	= .004633/100, #disease induced death for young
		vacy		= .001,
		aqy		= .0, #Aqcuired immunity rate
		birthy	= .0, #This can be lower if Civ vaccination is Recovered
		civy		= .000, #Modernization factor
		
		aging = 1/365/14,
		dehyd = 0, #percentage of vaccines that do not work
		
		contact	= .142/1,  	
		recov 	= .028,   
		ldeath  = 6.8/365/1000, #disease independent death
		hdeath	= .00567/100, #disease induced death for old
		vac		= 0,
		aq		= .0,
		birth	= 24.3/365/1000, #This can be lower if Civ vaccination is Recovered; dS.dt
		civ		= .0)

#SVIR w/vac+extras(Dehyd Modern)
pars <- c(
		contacty= 0.190/1,  	
		recovy 	= 0.028,   
		ldeathy = 6.8/365/1000, #disease independent death
		hdeathy	= .004633/100, #disease induced death for young. Approx twice as likely to die
		vacy	= .001, #.0002
		aqy		= .0, #Aqcuired immunity rate
		birthy	= .0, #This can be lower if Civ vaccination is Recovered
		civy	= .0000002, #Modernization factor
		
		aging 	= 1/365/14,
		dehyd 	=  .01, #percentage of vaccines that do not work
		
		contact	= .142/1,  	
		recov 	= .028,   
		ldeath  = 6.8/365/1000, #disease independent death
		hdeath	= .00567/100, #disease induced death for old approx twice as likely to die
		vac		= 0,
		aq		= .0,
		birth	= 24.3/365/1000, #This can be lower if Civ vaccination is Recovered; dS.dt
		civ		= .0,

		vaccycle=floor(365/4),#currently about 4 rounds of vac. per year in pakistan-> 60 total from 93 to 2007
		vaccyclelength=5 )

#200 cases of polio in 2011-> taking into account Pakistan's pop and with .3% of cases reported:
#I=~.0004=4E-4 currently. 
#So I=~ 1E-6 should be eradication condition

youngperc=.37
oldperc=1-youngperc
init.values <- c(
		Sy = .495 *youngperc, Vy = .00*youngperc,
		Iy = .01*youngperc, Ry = .495*youngperc,
		
		Cy = .0, Hy = 0, Ly = 0,
		
		S = .495*oldperc, V = .00*oldperc,
		I = .01*oldperc, R = .495*oldperc,
		
		C = .0, H = 0, L = 0,M=0)
#Tested with rough equilibria values and it looks weird. Mostly because we start with no vaccinated. 
times <- seq(0, 15*365, by = 1)

SVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				
				if (time%%vaccycle <=vaccyclelength){
					dSy.dt = birth*(S+V+I+R)+Sy*(birthy-vacy-contacty*(Iy+I)-ldeathy-aging)+dehyd*Vy+(Cy)*Ry 
					dVy.dt = vacy*Sy-ldeathy*Vy-aging*Vy-dehyd*Vy
					dM.dt  = vacy*(Sy+Iy+Ry)	
				}
				else{
					dSy.dt = birth*(S+V+I+R)+Sy*(birthy-vacy-contacty*(Iy+I)-ldeathy-aging)+(Cy)*Ry
					dVy.dt = -ldeathy*Vy-aging*Vy
					dM.dt  = 0	
				}
				dCy.dt = civy
				dIy.dt = contacty*(I+Iy)*Sy-(hdeathy+hdeathy*dehyd)*Iy-recovy*Iy-aging*Iy
				dRy.dt = recovy*Iy-hdeathy*Ry-aging*Ry-Cy*Ry
				dHy.dt = hdeathy*(Ry+(1+dehyd)*Iy)
				dLy.dt = ldeathy*(Sy+Vy)
				
				dS.dt =	-aq*S-(vac)*S-contact*(I+Iy)*S-ldeath*S-C*S+aging*Sy
				dV.dt = (vac)*S-ldeath*V+aging*Vy
				dI.dt = contact*(I+Iy)*S-hdeath*I*(1+dehyd)-recov*I+aging*Iy
				dR.dt = recov*I-hdeath*R+aging*Ry+(aq-C)*S
				dC.dt = civ
				dH.dt = hdeath*(R+I*(1+dehyd))
				dL.dt = ldeath*(S+V)
				
				return(list(c(dSy.dt,dVy.dt,dIy.dt,dRy.dt,dCy.dt,
										dHy.dt,dLy.dt,dS.dt,dV.dt, dI.dt, dR.dt,
										dC.dt,dH.dt,dL.dt,dM.dt)))
				

			})
}


#	Having defined everything, now we ask the program ode
#	to actually solve the system:

out = as.data.frame(ode(func = SVIR, y = init.values, 
				parms = pars, times = times))

#---------------------------------------------------
#	Plot of SVIR Model of Young Population
x11()
par(mfrow = c(2,2))

matplot(out$time/365, out[c(2:5)], type = "l", xlab = "Time (years)", 
		ylab = "Percent of Original Population", main = "SVIR Model of Young with 
				Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Vacinnated", 
				"Infectious", "Recovered"),
		col = 1:5, lty = 1:5)

#---------------------------------------------------
#  Plot of SVIR Model of Old Population


matplot(out$time/365, out[c(9:12)], type = "l", xlab = "Time (years)", 
		ylab = "Percent of Original Population", main = "SVIR Model of Old 
				with Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Vacinnated", 
				"Infectious", "Recovered"),
		col = 1:5, lty = 1:5)
#---------------------------------------------------
#  Plot of Death of Young Population


matplot(out$time/365, out[c(7,8)], type = "l", xlab = "Time (years)", 
		ylab = "num", main = "Deaths of Young", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent"),
		col = 1:2, lty = 1:2)
#--------------------------------------------------
#  Plot of Death of Old Population


matplot(out$time/365, out[c(14,15,16)], type = "l", xlab = "Time (years)", 
		ylab = "num", main = "Deaths of Old", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent, Vaccinations given"),
		col = 1:3, lty = 1:3)

#x11()
#matplot(out$time/365, out[c(16)], type = "l", xlab = "Time (years)", 
#		ylab = "num", main = "Vaccinations given", lwd = 2)
#
#legend("topleft", c("Vaccinations given"),
#		col = 1, lty = 1)

head(out) # last six values of all states
tail(out)

