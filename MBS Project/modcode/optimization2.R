#This was the emailed copy

## Load in the DiffEq solver

library(deSolve)

# Clear the memory

rm(list=ls())
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
times <- seq(0, 10*365, by = 1)

bestM=9000000

	print(i)
SVIR <- function(time, y.values, parameters) {
	
	with(as.list(c(y.values, parameters)), {
				
				
				if((Iy+I)/(S+V+I+R+Sy+Vy+Iy+Ry)>.01 & time>55 & done==0){
					break
					done=0
					
					if(M<bestM){
					bestM=M
					bestvac=(i)/1000
					print(bestvac)
					print(bestM)
					print(time)
					print('endline')
						}
					}
				else{print('h')}
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
				dH.dt = hdeath*(R+I(1+dehyd))
				dL.dt = ldeath*(S+V)
				
				return(list(c(dSy.dt,dVy.dt,dIy.dt,dRy.dt,dCy.dt,
										dHy.dt,dLy.dt,dS.dt,dV.dt, dI.dt, dR.dt,
										dC.dt,dH.dt,dL.dt,dM.dt)))			
				})
}


for(i in 1:3){
	pars <- c(
			done=1,
			contacty= 0.190/1,  	
			recovy 	= 0.028,   
			ldeathy = 6.8/365/1000, #disease independent death
			hdeathy	= .004633/100, #disease induced death for young. Approx twice as likely to die
			vacy	= i/1000,
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
	
print(ode(func = SVIR, y = init.values, 
					parms = pars, times = times))
	
}


