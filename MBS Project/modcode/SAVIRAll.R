#This was the emailed copy

## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

pars <- c(
		contacty= 0.1,  	
		recovy 	= 0.03,   
		ldeathy  = .03, #disease independent death
		hdeathy	= .06, #disease induced death for young
		vacy		= .05,
		aqy		= .03, #Aqcuired immunity rate
		birthy	= .0, #This can be lower if Civ vaccination is removed
		civy		= .0001, #Modernization factor
		
		aging = .01,
		dehyd = .1, #percentage of vaccines that do not work
		
		contact	= .05,  	
		recov 	= .02,   
		ldeath  = .06, #disease independent death
		hdeath	= .08, #disease induced death for old
		vac		= .07,
		aq		= .03,
		birth	= .1, #This can be lower if Civ vaccination is removed; dS.dt
		civ		= .0001)	

init.values <- c(Sy = .3, Ay = .3, Vy = .005, Iy = .05, Ry = 0,
		Cy = .0, Hy = 0, Ly = 0, S = .6, A = .3, V = .005,
		I = .005, R = 0, C = .0, H = 0, L = 0)
times <- seq(0, 1000, by = 1)
SAVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				dSy.dt = birth*(S+A+V+I+R)+birthy*Sy-aqy*Sy-(vac)*Sy-contacty*Iy*Sy-ldeathy*
						Sy-Cy/75*Sy-aging*Sy+dehyd*Vy
				dAy.dt = aqy*Sy-Cy*Ay-hdeathy*Ay-aging*Ay
				dVy.dt = (vacy)*Sy-ldeathy*Vy+Cy/75*Sy-aging*Vy-dehyd*Vy
				dIy.dt = contacty*(I+Iy)*Sy-hdeathy*Iy-recovy*Iy-aging*Iy
				dRy.dt = recovy*Iy-hdeathy*Ry-aging*Ry
				dCy.dt = civy
				dHy.dt = hdeathy*(Ry+Iy+Ay)
				dLy.dt = ldeathy*(Sy+Vy)
				
				dS.dt =	birth*S-aq*S-(vac)*S-contact*I*S-ldeath*S-C/75*S+
						aging*Sy+dehyd*V
				dA.dt = aq*S-C*A-hdeath*A+aging*Ay
				dV.dt = (vac)*S-ldeath*V+C/75*S+aging*Vy-dehyd*V
				dI.dt = contact*(I+Iy)*S-hdeath*I-recov*I+aging*Iy
				dR.dt = recov*I-hdeath*R+aging*Ry
				dC.dt = civ
				dH.dt = hdeath*(R+I+A)
				dL.dt = ldeath*(S+V)
				
				return(list(c(dSy.dt,dAy.dt,dVy.dt,dIy.dt,dRy.dt,dCy.dt,
										dHy.dt,dLy.dt,dS.dt,dA.dt,dV.dt, dI.dt, dR.dt,
										dC.dt,dH.dt,dL.dt)))
			})
}

#	Having defined everything, now we ask the program ode
#	to actually solve the system:

out = as.data.frame(ode(func = SAVIR, y = init.values, 
				parms = pars, times = times))

#---------------------------------------------------
#	Plot of SAVIR Model of Young Population

x11()
matplot(out$time, out[c(2:6)], type = "l", xlab = "Time", 
		ylab = "Percent of Original Population", main = "SAVIR Model of Young with 
				Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Aqcuired Immunity","Vacinnated", 
				"Infectious", "Removed"),
		col = 1:5, lty = 1:5)

#---------------------------------------------------
#  Plot of SAVIR Model of Old Population

x11()
matplot(out$time, out[c(10:14)], type = "l", xlab = "Time", 
		ylab = "Percent of Original Population", main = "SAVIR Model of Old 
				with Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Aqcuired Immunity","Vacinnated", 
				"Infectious", "Removed"),
		col = 1:5, lty = 1:5)
#---------------------------------------------------
#  Plot of Death of Young Population

x11()
matplot(out$time, out[c(8,9)], type = "l", xlab = "time", 
		ylab = "num", main = "Deaths of Young", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent"),
		col = 1:2, lty = 1:2)
#---------------------------------------------------
#  Plot of Death of Old Population

x11()
matplot(out$time, out[c(16,17)], type = "l", xlab = "time", 
		ylab = "num", main = "Deaths of Old", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent"),
		col = 1:2, lty = 1:2)


tail(out) # last six values of all states
