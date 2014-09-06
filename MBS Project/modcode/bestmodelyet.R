#This was the emailed copy

## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())

#	First, the parameters for our model:

pars <- c(
		contacty= 0.190,  	
		recovy 	= 0.028,   
		ldeathy  = (8/1000)/365, #disease independent death
		hdeathy	= (.004633)/365, #disease induced death for young
		vacy		= 0,
		aqy		= .03, #Aqcuired immunity rate
		civy		= .0001, #Modernization factor
		
		aging = (1/365)*(1/14),
		dehyd = .1, #percentage of vaccines that do not work
		
		contact	= .142,  	
		recov 	= .028,   
		ldeath  = (8/1000)/365, #disease independent death
		hdeath	= (.00567)/365, #disease induced death for old
		vac		= 0,
		aq		= .03,
		birth	= (31/1000)/365, #This can be lower if Civ vaccination is removed; dS.dt
		civ		= .0001)	

init.values <- c(Sy = .245, Vy = .005, Iy = .05, Ry = .05,
		Cy = .0, Hy = 0, Ly = 0, S = .60, V = 0,
		I = 0, R = .05, C = .0, H = 0, L = 0)
times <- seq(0, 1000, by = 1)
SVIR <- function(time, y.values, parameters) {
	with(as.list(c(y.values, parameters)), {
				
				dSy.dt = birth*(S+V+I+R)-(vac)*Sy-contacty*Iy*Sy-ldeathy*Sy-Cy*Sy-aging*Sy+dehyd*Vy
				dVy.dt = (vacy)*Sy-ldeathy*Vy+Cy*Sy-aging*Vy-dehyd*Vy
				dIy.dt = contacty*(I+Iy)*Sy-(hdeathy+ldeath)*Iy-recovy*Iy-aging*Iy
				dRy.dt = recovy*Iy-(hdeathy+ldeath)*Ry-aging*Ry+(-Cy)*Sy
				dCy.dt = civy
				dHy.dt = hdeathy*(Ry+Iy)
				dLy.dt = ldeathy*(Sy+Vy)
				
				dS.dt =	-aq*S-(vac)*S-contact*I*S-ldeath*S-C*S+aging*Sy+dehyd*V
				dV.dt = (vac)*S-ldeath*V+C*S+aging*Vy-dehyd*V
				dI.dt = contact*(I+Iy)*S-(hdeath+ldeath)*I-recov*I+aging*Iy
				dR.dt = recov*I-(hdeath+ldeath)*R+aging*Ry+(aq-C)*S
				dC.dt = civ
				dH.dt = hdeath*(R+I)
				dL.dt = ldeath*(S+V)
				
				return(list(c(dSy.dt,dVy.dt,dIy.dt,dRy.dt,dCy.dt,
										dHy.dt,dLy.dt,dS.dt,dV.dt, dI.dt, dR.dt,
										dC.dt,dH.dt,dL.dt)))
			})
}

#	Having defined everything, now we ask the program ode
#	to actually solve the system:

out = as.data.frame(ode(func = SVIR, y = init.values, 
				parms = pars, times = times))

#---------------------------------------------------
#	Plot of SAVIR Model of Young Population

x11()
matplot(out$time, out[c(2:5)], type = "l", xlab = "Time", 
		ylab = "Percent of Original Population", main = "SAVIR Model of Young with 
				Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Vacinnated", 
				"Infectious", "Removed"),
		col = 1:5, lty = 1:5)

#---------------------------------------------------
#  Plot of SAVIR Model of Old Population

x11()
matplot(out$time, out[c(9:12)], type = "l", xlab = "Time", 
		ylab = "Percent of Original Population", main = "SAVIR Model of Old 
				with Birth/Death and Dehydration", lwd = 2)

legend("topright", c("Susceptible","Vacinnated", 
				"Infectious", "Removed"),
		col = 1:5, lty = 1:5)
#---------------------------------------------------
#  Plot of Death of Young Population

x11()
matplot(out$time, out[c(7,8)], type = "l", xlab = "time", 
		ylab = "num", main = "Deaths of Young", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent"),
		col = 1:2, lty = 1:2)
#-----------	----------------------------------------
#  Plot of Death of Old Population

x11()
matplot(out$time, out[c(14,15)], type = "l", xlab = "time", 
		ylab = "num", main = "Deaths of Old", lwd = 2)

legend("topleft", c("Disease Dependent","Disease Independent"),
		col = 1:2, lty = 1:2)


tail(out) # last six values of all states
