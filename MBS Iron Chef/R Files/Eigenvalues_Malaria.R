#	Using R to find solutions to systems of ODEs and to 
#	check those solutions for stability.
#
#	Chris Leary, 4 February 2013

#	WARNING:  The solvers can be glitchy!


#Install the rootSolve package

#install.packages("rootSolve")

#Load the rootSolve package

library(rootSolve)


# Clear the variables

rm(list=ls())

# 	Values of the parameters

vecdeath=.01
bites=.01
bitel=.01
biteh=.01
vecexposetime=.01
#birth and death
birth=0
sdeath=.0
edeath=.0
ildeath=.0
ihdeath=.0
#class transfers
etime=.01
iltime=.01	
ihtoil=.01

#	We define the right hand side of the system that we want to solve

rhs <- function(p)	{
  vecdeath=.13
  bites=.022
  bitel=.04
  biteh=.48
  vecexposetime=.091
  #birth and death
  birth=0
  sdeath=0
  edeath=0
  ildeath=.0009
  ihdeath=.0009
  #class transfers
  etime=5.4
  iltime=.01	
  ihtoil=.0035
	Sv <- p[1]
	Ev <- p[2]
	Iv <- p[3]
  S <- p[4]
  E <- p[5]
  IL <- p[6]
  IH <- p[7]
  D <- p[8]
	r <- rep(NA, length(p))
	r[1] <- vecdeath*(Ev+Iv)-Sv*bitel*IL-Sv*biteh*IH
	r[2] <- -vecdeath*Ev    +Sv*bitel*IL+Sv*biteh*IH -vecexposetime*Ev
	r[3] <- -vecdeath*Iv  	        			  +vecexposetime*Ev
  r[4] <- -sdeath*S  					          -bites*S*Iv		+birth*(S+E+IL+IH)
  r[5] <- -edeath*E   +bites*S*Iv 	 	 -etime*E
  r[6] <- -ildeath*IL  +etime*E+ihtoil*IH  -iltime*IL
  r[7] <- -ihdeath*IH  +iltime*IL			 -ihtoil*IH	
  r[8] <- ildeath*IL  +ihdeath*IH
	r
	}
	
p0 <- c(Sv = 1, Ev = 2.17*12/10, Iv = 1.52+.152*2, S = 0, E = 0, IL = 1, IH = 1, D = 0)	# An initial guess as to the solution

ans <- multiroot(f = rhs, start = p0)

ans

############

#	Define the Jacobian

Svstar = ans$root[1]
Evstar = ans$root[2]
Ivstar = ans$root[3]
Sstar  = ans$root[4]
Estar  = ans$root[5]
ILstar = ans$root[6]
IHstar = ans$root[7]
Dstar  = ans$root[8]

J = matrix( c(-bitel*ILstar-biteh*IHstar, vecdeath, vecdeath, 0, 0, -bitel*Svstar, -biteh*Svstar, 0,
              bitel*ILstar + bites*IHstar, -vecdeath - vecexposetime, 0, 0, 0, bitel*Svstar, biteh*Svstar, 0,
              0, vecexposetime, vecdeath, 0, 0, 0, 0, 0,
              0, 0, -bites*Sstar, -sdeath-bites+birth, birth, birth, birth, 0,
              0, 0, bites*Sstar, bites * Ivstar, -edeath-etime, 0, 0, 0,
              0, 0, 0, 0, etime, -iltime-ildeath, ihtoil, 0,
              0, 0, 0, 0, 0, iltime, -ihdeath-ihtoil, 0,
              0, 0, 0, 0 , 0, ildeath, ihdeath, 0), nrow = 8, byrow = TRUE)
			
eigen(J)

#	Remember:  The solution is stable if and only if every eigenvalue has a negative real part.