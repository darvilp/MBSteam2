Models:
Basic SIR
-Just the basic SIR with a contact rate and a recovery rate

SAVIR
-Adds a Vaccinated class and an Acquired immunity class that Susceptible can move into
-The system is still closed so SAV will drain into R

SAVIR BCD
-SAVIR with Births, 'Civilization' and Deaths
-The Civilization factor is a parameter that changes over time
-The C factor reduces the rate of acquired immunity to 0 eventually
-The C factor could add more people to the Vaccinated class
-Two death rates, a high one for anyone who has contracted the disease (AIR)
and a low death rate for anyone who has not (SV)

AgeSAVIR 
-Separates the SAVIR model into two age classes with separate sets of parms for each
-Has Births and deaths to move people from age groups

AgeSAVIR BCD
-Separates SAVIR BCD into two age groups with separate sets of parms for each


All models with deaths keep track of how people died in the variables H and L for 
high death and low deathrates.

Things that need to be done:
Births and Deaths models:
-Currently births are dS.dt= birth*S when they should be birth*(S+A+V+I+R)
need to figure out how to do this and get correct deaths parameters
Aging Models:
-Currently no interaction between age groups, so young only infect young, old only infect old. Should be a quick fix if no
problems come up. 
*Fixed
-Parameters need to be figured out for AgeSAVIR and then can be exported out into the other models
