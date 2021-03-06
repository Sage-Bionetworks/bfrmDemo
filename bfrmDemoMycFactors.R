## bfrmDemoMycFactors.R

##########
# RUN BFRM IN THE SPARSE ANOVA MODE
##########
require(synapseClient)
require(bfrm)
require(Biobase)

nMycEnt <- loadEntity('syn464157')
nMycEset <- nMycEnt$objects$nMycEset
mycPheno <- pData(nMycEset)

anovaResult <- bfrm(exprs(nMycEset), 
                    design = ifelse(mycPheno$treatment == 'MYC', 1, 0))
mPPib <- anovaResult@results$mPostPib
topProbeLogical <- mPPib[ , 2] >= 0.99
topProbeInd <- grep("TRUE", topProbeLogical)

##########
# RUN BFRM IN THE FACTOR DISCOVERY MODE
##########

mycEvolveFactor <- evolve(exprs(nMycEset), 
                           init = as.numeric(topProbeInd),
                           priorpsia = 2,
                           priorpsib = 0.005,
                           varThreshold = 0.85,
                           facThreshold = 0.95,
                           maxVarIter = 30,
                           minFacVars = 10,
                           maxFacVars = length(topProbeInd),
                           maxFacs = 50,
                           maxVars = length(topProbeInd)
                           )


##########
# VISUALIZE THE RESULT
##########

image(t(mycEvolveFactor@results$mF[ -1, ]))

## Each row represents the factor scores across all 18 samples (columns).
## An internal reality-check is to now take your factors and project them back
## on the original data and see how they compare to your result mF matrix

mFHat <- projection(mycEvolveFactor, exprs(nMycEset))

## Look at the projection heatmap

image(mFHat)


