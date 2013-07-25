nsim <- 20
popn <- numeric(nsim)
for (k in 1:nsim) {
environ.sim.dat<-nupoint.env.simulator(pars=c(60,10,50),
                                       z.mat=NULL,
                                       xlim=c(0,200),ylim=c(0,100),
                                       grid.resolution=1,grad.type='NORM',det.type='HNORM',
                                       observer.coords=c(100,0),nbr.targets=1000,
                                       environment.simulator.control=
                                         list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
                                       mask.mat=NULL,mask.ang=0,plot=FALSE,
                                       perp.lines=NULL,n=NULL)

test <- truncate(sightings=environ.sim.dat$sightings)
trunc.dist <- test$trunc.radius

# replace sightings inside environ.sim.dat with truncated sightings
environ.sim.dat[[1]] <- test$sightings
# parameter estimation
sim.norm.fit<-nupoint.env.fit(pars=c(50,10,40),
                              z=environ.sim.dat$sightings$z, 
                              rd=environ.sim.dat$sightings$d,
                              dzdy=environ.sim.dat$sightings$dzdy,
                              z.mat=environ.sim.dat$z.mat,
                              dzdy.mat=environ.sim.dat$zGradmat,
                              rd.mat=environ.sim.dat$rd.mat,
                              minz=min(environ.sim.dat$z.mat),
                              wx=environ.sim.dat$settings$xlim[2],
                              wy=environ.sim.dat$settings$ylim[2],
                              wz=max(environ.sim.dat$z.mat),
                              grad.type=environ.sim.dat$settings$grad.type,
                              det.type=environ.sim.dat$settings$det.type,
                              n=NULL,lower.b=rep(1,3),upper.b=rep(100,3))
#  estimate P for HT
#  truncate the grid at the truncation distance 
new.rdmat <- environ.sim.dat$rd.mat
new.zmat <- environ.sim.dat$z.mat
new.zgrad <- environ.sim.dat$zGradmat
for (i in 1:100) {
  for (j in 1:200) {
    if (new.rdmat[i,j]>trunc.dist) {
      new.rdmat[i,j] <- NA
      new.zmat[i,j] <- NA
      new.zgrad[i,j] <- NA
  }
  }
}
mat.g <- detectF(new.rdmat[!is.na(new.rdmat)], "HNORM", sim.norm.fit$par[3])
mat.pi <- pi.z.f("NORM", pars=sim.norm.fit$par[1:2], z=new.zmat[!is.na(new.zmat)], 
                 z.lim=c(min(new.zmat, na.rm=TRUE), max(new.zmat, na.rm=TRUE)))
#  Abundance within truncation zone
Nhat.a <- dim(environ.sim.dat$sightings)[1]/sum(mat.g*mat.pi*abs(new.zgrad[!is.na(new.zgrad)])/(1*environ.sim.dat$settings$xlim[2]))
#   Scale Nhat.a to entire study area by dividing by integral pi(x,y) in region a
divisor <- sum(mat.pi*abs(new.zgrad[!is.na(new.zgrad)])/(1*environ.sim.dat$settings$xlim[2]))
Nhat.region <- Nhat.a / divisor
popn[k] <- Nhat.region
}
