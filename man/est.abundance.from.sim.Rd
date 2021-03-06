\name{est.abundance.from.sim}
\alias{est.abundance.from.sim}
\title{Estimate animal abundance only in the circumstance where data come from simulation}
\usage{
  est.abundance.from.sim(environ.sim.dat, trunc.prop = 0.9,
    det.model = NULL, gradient.model = NULL)
}
\arguments{
  \item{fit.obj}{fitted object}

  \item{truncation}{distance proportion (default 0.9) such
  that sightings beyond 0.9*max.r are deleted}

  \item{det.model}{detection model under which to perform
  estimation, if NULL (default), then estimation is
  performed under true model}

  \item{gradient.model}{gradient model under which to
  perform estimation, if NULL (default), estimation is
  performed under true model}
}
\value{
  list containing abundance estimate within covered region
  and abundance estimate for entire study area (assuming
  grid cells are unit square in area)
}
\description{
  This function estimates animal abundance within the study
  area (grid) by calculating density \eqn{\pi (z(x,y))} as
  a function of covariate for each grid cell.
}
\details{
  Calls to appropriate distribution (normal, lognormal,
  beta, uniform, mixture of normals) in association with
  the parameters estimated by the likelihood routine
  (\code{nupoint.env.fit}) are summed to produce estimate.

  Should your grid cell sizes not be unit square, then
  multiply the value returned by this function by the grid
  cell size to produce abundance estimate in the units
  appropriate for your study.
}
\examples{
mysim <- nupoint.env.simulator(pars=c(60,10,50),z.mat=NULL,xlim=c(0,200),ylim=c(0,100),
           grid.resolution=1,grad.type='NORM',det.type='HNORM',
           observer.coords=c(100,0),nbr.targets=350,
           environment.simulator.control=list(c(X=50,Y=10,sd=60),c(X=90,Y=0,sd=30)),
           mask.mat=NULL,mask.ang=0,plot=TRUE,
           perp.lines=NULL,n=NULL)
 correct.model <- est.abundance.from.sim(mysim, trunc.prop=0.9, det.model=NULL, gradient.model=NULL)
 wrong.uni <- est.abundance.from.sim(mysim, trunc.prop=0.9, det.model=NULL, gradient.model="UNIFORM")
 wrong.beta <- est.abundance.from.sim(mysim, trunc.prop=0.9, det.model=NULL, gradient.model="BETA")
 wrong.lognorm <- est.abundance.from.sim(mysim, trunc.prop=0.9, det.model=NULL, gradient.model="LOGNORM")
#  compare AIC scores among the 4 models; which is lowest?
}
\author{
  Eric Rexstad
}
\references{
  M.J. Cox, D.L. Borchers, D.A. Demer, G.R. Cutter, and
  A.S. Brierley. 2011. Estimating the density of Antarctic
  krill (Euphausia superba) from multi-beam echo-sounder
  observations using distance sampling methods. Journal of
  the Royal Statistical Society: Series C (Applied
  Statistics), 60(2):301-316.

  M.J. Cox, D.L. Borchers and N. Kelly. 2013. nupoint: An R
  package for density estimation from point transects in
  the presence of non-uniform animal density Methods in
  Ecology and Evolution 4(6):589-594

  Marques, T.A. , Buckland, S.T. , Borchers, D.L. , Tosh,
  D. and McDonald, R.A. 2010.  Point transect sampling
  along linear features.  Biometrics 66(4):1247-1255.
}

