library(secr)
library(rgdal)

rm(list = ls())

#### import data ####
# traps
CVtraps1 <- read.traps(data = read.csv('traps1.csv', row.names = 1),
                       detector = 'count')
CVtraps2 <- read.traps(data = read.csv('traps2.csv', row.names = 1),
                       detector = 'count')
CVtraps3 <- read.traps(data = read.csv('traps3.csv', row.names = 1),
                       detector = 'count')
CVtraps4 <- read.traps(data = read.csv('traps4.csv', row.names = 1),
                       detector = 'count')

# capture history (rank 1-3)
CVcaptures2107 <- read.csv('CV_caps1_r3.csv')
CVdata1 <- make.capthist(captures = CVcaptures2107, traps = CVtraps1, fmt = 'trapID') 

CVcaptures2201 <- read.csv('CV_caps2_r3.csv')
CVdata2 <- make.capthist(captures = CVcaptures2201, traps = CVtraps2, fmt = 'trapID') 

CVcaptures2207 <- read.csv('CV_caps3_r3.csv')
CVdata3 <- make.capthist(captures = CVcaptures2207, traps = CVtraps3, fmt = 'trapID') 

CVcaptures2301 <- read.csv('CV_caps4_r3.csv')
CVdata4 <- make.capthist(captures = CVcaptures2301, traps = CVtraps4, fmt = 'trapID') 

# summary
summary(CVdata1)
summary(CVdata2)
summary(CVdata3)
summary(CVdata4)

# closure test
# if p > 0.05, closure can be assumed
closure.test(CVdata1)
closure.test(CVdata2)
closure.test(CVdata3)
closure.test(CVdata4)

# movement parameters
CVdist01 <- RPSV(CVdata1, CC = TRUE)
CVdist02 <- RPSV(CVdata2, CC = TRUE)
CVdist03 <- RPSV(CVdata3, CC = TRUE)
CVdist04 <- RPSV(CVdata4, CC = TRUE)
CV_RPSV <- c(CVdist01, CVdist02, CVdist03, CVdist04)

CVdist1 <- MMDM(CVdata1)
CVdist2 <- MMDM(CVdata2)
CVdist3 <- MMDM(CVdata3)
CVdist4 <- MMDM(CVdata4)
CV_MMDM <- c(CVdist1, CVdist2, CVdist3, CVdist4)

#####
# standard mask
CVmask0 <- make.mask(CVtraps1, buffer = 4*max(CV_RPSV), type ='trapbuffer')

# add elevation as a covariate
elev <- rast('ELEV.asc')
plot(elev)

elevation <- addCovariates(CVmask0, elev)
head(covariates(elevation))


#### models with standard mask ####
# period 1
CVcounts_hn01 <- secr.fit(CVdata1, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0) # no covariate

CVcounts_hn01 <- secr.fit(CVdata1, model=list(D~ELEV, g0~1, sigma~1), 
                          detectfn = 'HN', mask = elevation) # elevation as covariate

coefficients(CVcounts_hn01)
CVpred01 <- predict(CVcounts_hn01)
CVnumb01 <- region.N(CVcounts_hn01)

# period 2
CVcounts_hn02 <- secr.fit(CVdata2, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)

CVcounts_hn02 <- secr.fit(CVdata2, model=list(D~ELEV, g0~1, sigma~1), 
                          detectfn = 'HN', mask = elevation)

coefficients(CVcounts_hn02)
CVpred02 <- predict(CVcounts_hn02)
CVnumb02 <- region.N(CVcounts_hn02)

# period 3
CVcounts_hn03 <- secr.fit(CVdata3, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)

CVcounts_hn03 <- secr.fit(CVdata3, model=list(D~ELEV, g0~1, sigma~1), 
                          detectfn = 'HN', mask = elevation)

coefficients(CVcounts_hn03)
CVpred03 <- predict(CVcounts_hn03)
CVnumb03 <- region.N(CVcounts_hn03)

# period 4
CVcounts_hn04 <- secr.fit(CVdata4, model=list(D~1, g0~1, sigma~1), 
                         detectfn = 'HN', mask = CVmask0)

CVcounts_hn04 <- secr.fit(CVdata4, model=list(D~ELEV, g0~1, sigma~1), 
                          detectfn = 'HN', mask = elevation)

coefficients(CVcounts_hn04)
CVpred04 <- predict(CVcounts_hn04)
CVnumb04 <- region.N(CVcounts_hn04)

#####

# import habitat
hbt1 <- readOGR('.', 'CV_mask') # buffer = 0.5*MMDM
CVmask1 <- make.mask(traps(CVdata1), type='polybuffer', 
                     poly = hbt1, poly.habitat = T)
plot(CVmask1)

#### models with SDM_derived mask ####
# period 1
CVcounts_hn1 <- secr.fit(CVdata1, detectfn = 'HN', 
                         mask = CVmask1, start = autoini(CVdata1,CVmask1))

coefficients(CVcounts_hn1)
CVpred1 <- predict(CVcounts_hn1)
CVnumb1 <- region.N(CVcounts_hn1)

CVsur1 <- fx.total(CVcounts_hn1)
plot(CVsur1, covariate = 'D.sum', scale = 1000)


# period 2
CVcounts_hn2 <- secr.fit(CVdata2, detectfn = 'HN', 
                         mask = CVmask1, start = autoini(CVdata2,CVmask1))

coefficients(CVcounts_hn2)
CVpred2 <- predict(CVcounts_hn2)
CVnumb2 <- region.N(CVcounts_hn2)

CVsur2 <- fx.total(CVcounts_hn2)
plot(CVsur2, covariate = 'D.sum', scale = 1000)


# period 3
CVcounts_hn3 <- secr.fit(CVdata3, detectfn = 'HN', 
                         mask = CVmask1, start = autoini(CVdata3,CVmask1))

coefficients(CVcounts_hn3)
CVpred3 <- predict(CVcounts_hn3)
CVnumb3 <- region.N(CVcounts_hn3)

CVsur3 <- fx.total(CVcounts_hn3)
plot(CVsur3, covariate = 'D.sum', scale = 1000)


# period 4
CVcounts_hn4 <- secr.fit(CVdata4, detectfn = 'HN', 
                         mask = CVmask1, start = autoini(CVdata4,CVmask1))

coefficients(CVcounts_hn4)
CVpred4 <- predict(CVcounts_hn4)
CVnumb4 <- region.N(CVcounts_hn4)

CVsur4 <- fx.total(CVcounts_hn4)
plot(CVsur4, covariate = 'D.sum', scale = 1000)
