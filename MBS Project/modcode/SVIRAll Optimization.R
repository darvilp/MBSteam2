#This was the emailed copy

## Load in the DiffEq solver
require(compiler) 
enableJIT(3)
library(deSolve)

# Clear the memory

rm(list=ls())
bestM=40
bestvacy=40
besttime=5000000
bestH=800000000
vacyloopmax=100000#will go through values 1/vacyloomax to 1 for vacy
vacyseq=seq(.00005,.0005,by=.00005)
vaccyclelengthseq=seq(1,28,1)
vaccycleseq=seq(1,12,1)

for(j in vaccycleseq){
for(k in vaccyclelengthseq){
for(i in vacyseq){	
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
					
					dS.dt =	-aq*S-(vac)*S-contact*(I+Iy)*S-ldeath*S-0*S+aging*Sy
					dV.dt = (vac)*S-ldeath*V+aging*Vy
					dI.dt = contact*(I+Iy)*S-hdeath*I*(1+dehyd)-recov*I+aging*Iy
					dR.dt = recov*I-hdeath*R+aging*Ry+(aq-0)*S
					dC.dt = civ
					dH.dt = hdeath*(R+I*(1+dehyd))
					dL.dt = ldeath*(S+V)
					
					return(list(c(dSy.dt,dVy.dt,dIy.dt,dRy.dt,dCy.dt,
											dHy.dt,dLy.dt,dS.dt,dV.dt, dI.dt, dR.dt,
											dC.dt,dH.dt,dL.dt,dM.dt)))		
				})
	}
	
	#SVIR w/vac+extras(Dehyd Modern)
	pars <- c(
			contacty= 0.190/1,  	
			recovy 	= 0.028,   
			ldeathy = 6.8/365/1000, #disease independent death
			hdeathy	= .004633/100, #disease induced death for young. Approx twice as likely to die
			vacy	= i,
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
			
			vaccycle=floor(365/j),#currently about 4 rounds of vac. per year in pakistan-> 60 total from 93 to 2007
			vaccyclelength=k )
	
	#200 cases of polio in 2011-> taking into account Pakistan's pop and with .3% of cases reported:
	#I=~.0004=4E-4 currently. 
	#So I=~ 1E-6 should be eradication condition
	
	youngperc=.37
	oldperc=1-youngperc
	init.values <- c(Sy = .495 *youngperc, Vy = .00*youngperc, Iy = .01*youngperc, Ry = .495*youngperc,
			Cy = .0, Hy = 0, Ly = 0,
			
			S = .495*oldperc, V = .00*oldperc,
			I = .01*oldperc, R = .495*oldperc, C = 9000000000, H = 0, L = 0,
			
			M=0)
	years=15
	times <- seq(0, years*365, by = 1)
	root <- function(t, y.values, parms) {
			
			return((((y.values[3]+y.values[10])/(y.values[1]+y.values[2]+y.values[3]+y.values[4]+y.values[8]+y.values[9]+y.values[10]+y.values[11]))-.0000000005))#When Iy is at -num then root is found #.000000005,0000000005 with 0 to 1 by .001 shows inflection
		}
	
		
	#Event happens when root is found. Set everything to zero so the equation immediately solves. 
	#Keep M because that's what we want
	
	eventfun <- function(time, y.values, parameters){ 
		with (as.list(c(y.values, parameters)),{
					if(time>0 & time<C ){
					C=time
					}
					Sy = 0
					Vy = 0
					Iy=0
					Ry=0
					Cy=0
					Hy=Hy+Ry
					Ly=0
					S=0
					V=0
					I=0
					R=0
					H=H+R
					L=0
					M=M					
	return(c(Sy,Vy,Iy,Ry,Cy,Hy,Ly,S,V,I,R,C,H,L,M))
				})
	}
	out = as.data.frame(ode(func = SVIR, y = init.values, 
					parms = pars, times = times,events = list(func = eventfun, root = TRUE), rootfun = root),
					)
	


##---------------------------------------------------
#	#Plot of SVIR Model of Young Population
#	x11()
#	par(mfrow = c(2,2))
#	
#	matplot(out$time/365, out[c(2:5)], type = "l", xlab = "Time", 
#			ylab = "Percent of Original Population", main = "SVIR Model of Young with 
#					Birth/Death and Dehydration", lwd = 2)
#	
#	legend("topright", c("Susceptible","Vacinnated", 
#					"Infectious", "Recovered"),
#			col = 1:5, lty = 1:5)
#	
#	#---------------------------------------------------
#	#  Plot of SVIR Model of Old Population
#	
#	
#	matplot(out$time/365, out[c(9:12)], type = "l", xlab = "Time", 
#			ylab = "Percent of Original Population", main = "SVIR Model of Old 
#					with Birth/Death and Dehydration", lwd = 2)
#	
#	legend("topright", c("Susceptible","Vacinnated", 
#					"Infectious", "Recovered"),
#			col = 1:5, lty = 1:5)
#	#---------------------------------------------------
#	#  Plot of Death of Young Population
#	
#	
#	matplot(out$time/365, out[c(7,8)], type = "l", xlab = "time", 
#			ylab = "num", main = "Deaths of Young", lwd = 2)
#	
#	legend("topleft", c("Disease Dependent","Disease Independent"),
#			col = 1:2, lty = 1:2)
#	#-----------	----------------------------------------
#	#  Plot of Death of Old Population
#	
#	
#	matplot(out$time/365, out[c(14,15,16)], type = "l", xlab = "time", 
#			ylab = "num", main = "Deaths of Old", lwd = 2)
#	
#	legend("topleft", c("Disease Dependent","Disease Independent, Vaccinations given"),
#			col = 1:3, lty = 1:3)
	roottime=out$C[years*365]
	endH=out$H[years*365]+out$Hy[years*365]
	endM=out$M[years*365]
	write.table(cbind(roottime,i,j,k,endH,endM), "model.csv",append=TRUE,col.names=FALSE,sep =',')
	
	
	if(out$M[years*365]<bestM & roottime>0 & roottime<(years*365-5) ){	
		bestM=out$M[years*365]
		bestim=i
		bestjm=j
		bestkm=k
		besttm=roottime
		besthm=endH
		print(cat('BESTM',bestM,'was obtained with ijk',i,j,k, 'at', roottime))
	}
	
	if(roottime<besttime & roottime>0 & roottime<(years*365-5) ){	
		
		besttime=roottime
		bestit=i
		bestjt=j
		bestkt=k
		bestmt=out$M[years*365]
		bestht=endH
		print(cat('BESTTIME!',out$M[years*365],'was obtained with ijk',i,j,k, 'at', roottime,'with H',endH))
	}
	
	if((endH)<bestH & roottime>0 & roottime<(years*365-5)){
		
		besth=out$H[years*365]+out$Hy[years*365]
		bestiH=i
		bestjH=j
		bestkH=k
		bestmH=out$M[years*365]
		besttH=roottime
		print(cat('BESTH',out$M[years*365],'was obtained with ijk',i,j,k,'at',roottime,'with H',endH))
	}
	else{
		if(roottime>0){
			print(cat(out$M[years*365],'didnt make the cut with ijk',i,j,k,'at',roottime,'with H',endH))}
	}
}
}
}
print(cat('BESTM',bestM,'was obtained with ijk',bestim,bestjm,bestkm, 'at', besttm,'with',besthm))
print(cat('BESTTIME',bestmt,'was obtained with ijk',bestit,bestjt,bestkt, 'at', besttime,'with',bestht))
print(cat('BESTH',bestmH,'was obtained with ijk',bestiH,bestjH,bestkH, 'at', besttH,'with H',endH))

#x11()
#matplot(out$time/365, out[c(16)], type = "l", xlab = "time", 
#		ylab = "num", main = "Vaccinations given", lwd = 2)
#
#legend("topleft", c("Vaccinations given"),
#		col = 1, lty = 1)

#print(out[c(9:12)]) # last six values of all states
#tail(out)

