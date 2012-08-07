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

bCatEvolveFactor <- evolve(exprs(nMycEset), 
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