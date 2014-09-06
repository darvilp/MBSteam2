library(deSolve)
x = array(12)
gamma = matrix(nrow=1, ncol=3)
alpha = matrix(nrow=1, ncol=3)
betaYY = 
betaYM = 
betaYO = 
betaMM = 
betaMY = 
betaMO = 
betaOO = 
betaOY = 
betaOM = 
maxtime = 60
SYo = .332
SMo = .505
SOo = .063
IYo = 1
IMo = 0
IOo = 0
RYo = 0
RMo = 0
ROo = 0
VYo = 0
VMo = 0
VOo = 0

Tmax = 100 # number of times steps (+1)
mod = function(t,x,parms) {
  SY = x[1]
  SM = x[2]
  SO = x[3]
  IY = x[4]
  IM = x[5]
  IO = x[6]
  RY = x[7]
  RM = x[8]
  RO = x[9] 
  with(as.list(parms) , {
    dSY = - (betaYY*SY*IY) - (betaYM*SY*IM) - (betaYO*SY*IO) - (alpha[1]*SY)
    dSM = - (betaMM*SM*IM) - (betaMY*SM*IY) - (betaMO*SM*IO) - (alpha[2]*SM)
    dSO = - (betaOO*SO*IO) - (betaOY*SO*IY) - (betaOM*SO*IM) - (alpha[3]*SO)
    dIY = (betaYY*SY*IY) + (betaYM*SY*IM) + (betaYO*SY*IO) - (gamma[1]*IY)
    dIM = (betaMM*SM*IM) + (betaMY*SM*IY) + (betaMO*SM*IO) - (gamma[2]*IM)
    dIO = (betaOO*SO*IO) + (betaOY*SO*IY) + (betaOM*SO*IM) - (gamma[3]*IO)
    dRY = gamma[1]*IY
    dRM = gamma[2]*IM
    dRO = gamma[3]*IO
    dVY = alpha[1]*SY
    dVM = alpha[2]*SM
    dVO = alpha[3]*SO
    res=c(dSY,dSM,dSO,dIY,dIM,dIO,dRY,dRM,dRO,dVY,dVM,dVO)
    list(res)
  })}
times = seq(0, 30, by = .1)
parms = c(gamma,alpha,betaYY,betaYM,betaYO,betaMM,betaMY,betaMO,betaOO,betaOY,betaOM)
xstart = c(SY = SYo, SM = SMo, SO = SOo, IY = IYo, IM = IMo, IO = IOo, RY = RYo, RM= RMo, RO = ROo, VY = VYo, VM = VMo, VO = VOo)
output = as.data.frame(ode(xstart, times, mod, parms))
gSY = output$SH
gSM = output$SL
gSO = output$IH
gIY = output$IY
gIM = output$IM
gIO = output$IO
gRY = output$RY
gRM = output$RM
gRO = output$RO
gVY = output$VY
gVM = output$VM
gVO = output$VO


plot(gSY,xlim=c(0,maxtime),ylim=c(0,10),type="l",lwd=2,col="green",
     ylab="Percent of Population",xlab="Time",
     cex.lab=1.5)
leg.txt=c("I High","I Low")
legend("topleft",leg.txt,lty=1,lwd=2,
       col=c("green","red"))
lines(gSM,xlim=c(1,maxtime),lwd=2,col="red")
lines(gSO,xlim=c(1,maxtime),lwd=2,col="red")
lines(gIY,xlim=c(1,maxtime),lwd=2,col="red")
lines(gIM,xlim=c(1,maxtime),lwd=2,col="red")
lines(gIO,xlim=c(1,maxtime),lwd=2,col="red")