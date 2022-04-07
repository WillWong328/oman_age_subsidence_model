library(modifiedBChron)
library(rethinking)
source("subsidence_analysis.R")







##--------------------------------------------------------------------UPDATING AGE-DEPTHS----------------------------------------------------------------------##




# Load data
df <- read.csv(file = "data_set.csv")
boundary_ages <- read.csv(file = "important_depths.csv")$depths


# Display and plot the original data set
print(df)

#pdf("OriginalDataset.pdf")                                    # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
plot(df$age,df$position*-1,col="blue",lwd=2,pch=19,
    xlab="Ages [Ma]",ylab="Depth [m]",main="Correlated Ages for the MIQRAT-1 Well of Oman",
    ylim=c(max(df$position*-1)+50,min(df$position*-1)-100),xlim=c(max(df$age)+10,min(df$age)-10))
arrows(x0=df$age ,x1=df$age ,y0=(df$position-df$thickness)*-1 ,y1=(df$position+df$thickness)*-1 , code=3, angle=90, length=0.05, col="blue", lwd=1) # Vertical error bars -> 2sigma uncertainty plotted
arrows(x0=df$age-df$ageSds, x1=df$age+df$ageSds ,y0=df$position*-1, y1=df$position*-1, code=3, angle=90, length=0.05, col="blue", lwd=1)            # Horizontal error bars -> 2sigma uncertainty plotted
legend("topleft", legend=c("Correlated Age-Depths",expression(paste("2",sigma," Uncertainty"))),                                                    # Legend titles
       col=c("blue", "blue"), pch=c(19,3))
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****





# Define our Bayesian inference functions for updating ages and depths based on the law of superposition.
likelihood_updateAD <- function(params) {
    ages <- params[,1]
    depths <- params[,2]

    singlelikelihoodsAges <- dnorm(df$age, mean = ages, sd = df$ageSds/2, log = T)
    singlelikelihoodsDepths <- dnorm(df$position*-1, mean = depths, sd = df$thickness/2, log = T)

    sumll <- sum(singlelikelihoodsAges) + sum(singlelikelihoodsDepths)
    return(sumll)
}

prior_updateAD <- function(params) {
    ages <- params[,1]
    depths <- params[,2]

    if(identical(order(ages),order(depths))) {
        return(log(1))
    } else {
        return(log(0))
    }
}

posterior_updateAD <- function(params) {
    return(likelihood_updateAD(params) + prior_updateAD(params))
}

proposal_updateAD <- function(params) {
    len <- length(params[,1])
    ages <- params[,1]
    depths <- params[,2]

    propAges <- rnorm(len, mean = ages, sd = df$ageSds/2)
    propDepths <- rnorm(len, mean = depths, sd = df$thickness/2)

    propADs <- cbind(propAges,propDepths)
    return(propADs)
}

metropolis_updateAD <- function(startvalues, iterations) {
    len <- length(startvalues[,1])
    chain <- array(dim = c(len,2,iterations+1))
    chain[,,1] <- startvalues

    for(i in seq(1,iterations)) {
        proposal <- proposal_updateAD(chain[,,i])

        probab <- exp(posterior_updateAD(proposal) - posterior_updateAD(chain[,,i]))
        if(is.nan(probab)) {
            chain[,,i+1] <- chain[,,i]
        } else if (runif(1) < probab) {
            chain[,,i+1] <- proposal
        } else {
            chain[,,i+1] <- chain[,,i]
        }
    }
    return(chain)
}



# MCMC
startvalues <- cbind(df$age,df$position*-1)
chain <- metropolis_updateAD(startvalues, 500000)



# Set up updated dataset
burnIn <- 10000
len <- length(chain[,1,1])

post_ages <- array(dim = c(len,1))
post_depths <- array(dim = c(len,1))
post_ages_HPDI <- array(dim = c(len,2))
post_depths_HPDI <- array(dim = c(len,2))
post_ages_uncertainties <- array(dim = c(len,2))
post_depths_uncertainties <- array(dim = c(len,2))
post_ages_sd <- array(dim = c(len,1))
post_depths_sd <- array(dim = c(len,1))

for(i in seq(1,len)) {
    post_ages[i] <- median(chain[i,1,-(1:burnIn)])
    post_depths[i] <- median(chain[i,2,-(1:burnIn)])
    post_ages_HPDI[i,] <- HPDI(chain[i,1,-(1:burnIn)],prob=0.95)
    post_depths_HPDI[i,] <- HPDI(chain[i,2,-(1:burnIn)],prob=0.95)
}

for(i in seq(1,len)) {
    post_ages_uncertainties[i,1] <- round(post_ages[i] - post_ages_HPDI[i,1], 2)
    post_ages_uncertainties[i,2] <- round(post_ages_HPDI[i,2] - post_ages[i], 2)
    post_ages_sd[i] <- median(post_ages_uncertainties[i,])
    post_depths_uncertainties[i,1] <- round(post_depths[i] - post_depths_HPDI[i,1], 2)
    post_depths_uncertainties[i,2] <- round(post_depths_HPDI[i,2] - post_depths[i], 2)
    post_depths_sd[i] <- median(post_depths_uncertainties[i,])
}

dfUpdated <- data.frame(
    ids <- df$ids,
    ages <- post_ages,
    ageSds <- post_ages_sd,
    position <- post_depths*-1,
    thickness <- post_depths_sd,
    distTypes <- df$distType
)



# Plot parameter histograms and markov chains

#pdf("HistUpdatedAges.pdf")                                     # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(4,3))
for(i in seq(1,len)) {
    hist(chain[i,1,-(1:burnIn)], nclass=30, xlab="Age [Ma]", main=df$ids[i], ylab=NULL, yaxt="n")
    mtext(paste(round(post_ages[i],2)," +",post_ages_uncertainties[i,2],"/ -",post_ages_uncertainties[i,1]),
        side=3, at=post_ages[i], cex=0.65, col="red")
    abline(v = df$age[i], col="blue", lwd=1.5)
    abline(v = post_ages[i], col="red", lwd=1.5)
    abline(v = post_ages_HPDI[i,1], col="red", lwd=1.5, lty=2)
    abline(v = post_ages_HPDI[i,2], col="red", lwd=1.5, lty=2)
}
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****




#pdf("HistUpdatedDepths.pdf")                                    # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(4,3))
for(i in seq(1,len)) {
    hist(chain[i,2,-(1:burnIn)], nclass=30, xlab="Depth [m]", main=df$ids[i], ylab=NULL, yaxt="n")
    mtext(paste(round(post_depths[i],2)," +",post_depths_uncertainties[i,2],"/ -",post_depths_uncertainties[i,1]),
        side=3, at=post_depths[i], cex=0.65, col="red")
    abline(v = df$position[i]*-1, col="blue", lwd=1.5)
    abline(v = post_depths[i], col="red", lwd=1.5)
    abline(v = post_depths_HPDI[i,1], col="red", lwd=1.5, lty=2)
    abline(v = post_depths_HPDI[i,2], col="red", lwd=1.5, lty=2)
}
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****





#pdf("ChainUpdatedAges.pdf")                                   # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(4,3))
for(i in seq(1,len)) {
    plot(chain[i,1,-(1:burnIn)], type="l", ylab="Age [Ma]", main=df$ids[i], xlab="Markov Chain Number")
    mtext(paste(round(post_ages[i],2)," +",post_ages_uncertainties[i,2],"/ -",post_ages_uncertainties[i,1]),
        side=3, cex=0.65, col="red")
    abline(h = df$age[i], col="blue", lwd=2)
    abline(h = post_ages[i], col="red", lwd=2)
    abline(h = post_ages_HPDI[i,1], col="red", lwd=2, lty=2)
    abline(h = post_ages_HPDI[i,2], col="red", lwd=2, lty=2)
}
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****





#pdf("ChainUpdatedDepths.pdf")                                  # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(4,3))
for(i in seq(1,len)) {
    plot(chain[i,2,-(1:burnIn)], type="l", ylab="Depth [m]", main=df$ids[i], xlab="Markov Chain Number")
    abline(h = df$position[i]*-1, col="blue", lwd=3.5)
    mtext(paste(round(post_depths[i],2)," +",post_depths_uncertainties[i,2],"/ -",post_depths_uncertainties[i,1]),
        side=3, cex=0.65, col="red")
    abline(h = post_depths[i], col="red", lwd=2)
    abline(h = post_depths_HPDI[i,1], col="red", lwd=2, lty=2)
    abline(h = post_depths_HPDI[i,2], col="red", lwd=2, lty=2)
}
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****






# Plot the updated dataset

#pdf("UpdatedAgeDepths.pdf")                                  # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
plot(df$age,df$position*-1,col="blue",lwd=2,pch=19,
    xlab="Ages [Ma]",ylab="Depth [m]",main="Updated Ages for the MIQRAT-1 Well of Oman",
    ylim=c(max(df$position*-1)+50,min(df$position*-1)-100),xlim=c(max(df$age)+10,min(df$age)-10))
points(post_ages,post_depths,col="red",lwd=2,pch=19)
arrows(x0=df$age ,x1=df$age ,y0=(df$position-df$thickness)*-1 ,y1=(df$position+df$thickness)*-1 , code=3, angle=90, length=0.05, col="blue", lwd=1)                          # Vertical error bars -> 2sigma uncertainty plotted
arrows(x0=df$age-df$ageSds, x1=df$age+df$ageSds ,y0=df$position*-1, y1=df$position*-1, code=3, angle=90, length=0.05, col="blue", lwd=1)                                     # Horizontal error bars -> 2sigma uncertainty plotted
arrows(x0=post_ages,x1=post_ages,y0=post_depths-post_depths_uncertainties[,1],y1=post_depths+post_depths_uncertainties[,2], code=3, angle=90, length=0.05, col="red", lwd=1) # Vertical error bars
arrows(x0=post_ages-post_ages_uncertainties[,1],x1=post_ages+post_ages_uncertainties[,2],y0=post_depths,y1=post_depths, code=3, angle=90, length=0.05, col="red", lwd=1)     # Horizontal error bars
legend("topleft", legend=c("Correlated Age-Depths [Mean]",expression(paste("2",sigma," Uncertainty")), NA, "Updated Age-Depths [Median]", "95% HDPI"),                       # Legend titles
       col=c("blue", "blue", NA, "red", "red"), pch=c(19,3, NA, 19,3))
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****











##----------------------------------------------------------------------MODIFIED BCHRON------------------------------------------------------------------------##




# Run Modified BChron
age_model <- ageModel(ages = dfUpdated$ages,
                  ageSds = dfUpdated$ageSds/2,
                  positions = dfUpdated$position,
                  ids = dfUpdated$ids,
                  positionThicknesses = dfUpdated$thickness/2,
                  distTypes = dfUpdated$distType,
                  predictPositions = seq(-4248, -3100, by = 1),
                  MC = 20000,      # Number of iterations
                  burn = 2000)     # Number of iterations to discard




# Plot the new chronology

#pdf("BayesianAgeModel.pdf")                                    # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
modelPlot(model = age_model,
           scale = 1,  # Changes the height of the probability distributions
           ylim = c(-4248,-3100),
           xlim = c(650,550),
           main = "Bayesian Age Model for the MIQRAT-1 Well of Oman",
           xlab = "Age [Ma]",
           ylab = "Stratigraphic Position [m]")
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****




# Estimate the ages at lithological boundaries
age_predictions <- agePredict(model = age_model,
           newPositions = boundary_ages,
           newPositionThicknesses = matrix(5,1,length(boundary_ages)))




# Re-plot the age model with boundary ages listed

#pdf("BayesianAgeModel.pdf")                                    # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
modelPlot(model = age_model,
            agePredictOutput = age_predictions,
            scale = 0.1,
            ylim = c(-4248,-3200),
            xlim = c(680,550),
            main = "Bayesian Age Model for the MIQRAT-1 Well of Oman",
            xlab = "Age [Ma]",
            ylab = "Stratigraphic Position [m]")
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****


# Display extracted ages and uncertainties
extracted_ages <- data.frame(
    Depth <- age_predictions$HDI[1]*-1,
    Median <- round(age_predictions$HDI[3],2),
    "95% Min" <- round(age_predictions$HDI[3]-age_predictions$HDI[2],4),
    "95% Max" <- round(age_predictions$HDI[4]-age_predictions$HDI[3],4)
)
numAges <- length(extracted_ages[,1])
extracted_ages





# Estimate the onset, nadir, and termination age of the SE
age_predictions_SE <- agePredict(model = age_model,
           newPositions = c(-3823,-3730,-3353),
           newPositionThicknesses = c(5,20,5))



# Re-plot the age model with SE estimates listed
SE_ages <- age_predictions_SE
SE_ages$HDI[,2:4] <- round(age_predictions_SE$HDI[,2:4],2)

#pdf("BayesianAgeModel.pdf")                                    # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
modelPlot(model = age_model,
            agePredictOutput = SE_ages,
            scale = 0.5,
            ylim = c(-4248,-3200),
            xlim = c(650,550),
            main = "Bayesian Age Model for the MIQRAT-1 Well of Oman",
            xlab = "Age [Ma]",
            ylab = "Stratigraphic Position [m]")
#dev.off()                                                      # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****




##----------------------------------------------------------------------SUBSIDENCE MODEL------------------------------------------------------------------------##




# Set up the upper and lower boundary ages/depths and uncertainties
baseAges <- age_predictions$HDI[2:16,3]
upperAges <- age_predictions$HDI[1:15,3]
baseDepths <- age_predictions$HDI[2:16,1]*-1
upperDepths <- age_predictions$HDI[1:15,1]*-1
decompAgesMin <- age_predictions$HDI[16:1,2]
decompAgesMax <- age_predictions$HDI[16:1,4]


# Create data frame for decompaction and subsidence models
strat_units <- c("Buah Dolostones", "Buah Limestones", "Shuram Limestones", "Shuram Interbedded-Limestones", "Shuram Interbedded-Shales",
                 "Shuram Shales", "Khufai Dolostones 1", "Khufai Shales 1", "Khufai Dolostones 2", "Khufai Limestones", "Khufai Dolostones 3",
                 "Khufai Shales 2", "Khufai Dolostones 4", "Masirah Bay-Shales", "Hadash Dolostones")

decomp_sub_data <- data.frame(
    "Stratigraphic Unit" = strat_units,
    "Surface Porosity" = c(0.2, 0.4, 0.4, 0.4, 0.63, 0.63, 0.2, 0.63, 0.2, 0.4, 0.2, 0.63, 0.2, 0.63, 0.2),
    "Porosity-Depth Coefficient" = c(0.6, 0.6, 0.6, 0.6, 0.51, 0.51, 0.6, 0.51, 0.6, 0.6, 0.6, 0.51, 0.6, 0.51, 0.6),
    "Grain Density" = c(2870, 2710, 2710, 2715, 2715, 2720, 2870, 2720, 2870, 2710, 2870, 2720, 2870, 2720, 2870),
    "Age (Upper) [Ma]" = upperAges,
    "Age (Base) [Ma]" = baseAges,
    "Present Depth (Upper) [m]" = upperDepths,
    "Present Depth (Base) [m]" = baseDepths,
    "Water Depth Correction [m]" = c(5, 10, 15, 15, 20, 30, 5, 10, 15, 20, 30, 40, 40, 40, 20),
    #"Water Depth Correction [m]" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    "Eustatic Correction [m]" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

# Display the input data
decomp_sub_data



# Decompact and perform backstripping on sedimentary column
results <- backstripping(decomp_sub_data)


# Display decompaction and backstripping results
results[1] # Decompacted thickness
results[2] # Age horizon thicknesses
results[3] # Average porosities
results[4] # Bulk densities
results[5] # Bulk density of column
results[6] # Display backstripped tectonic subsidence (values)



# Plot water-loaded tectonic subsidence
decompacted_thickness <- matrix(unlist(results[1]), nrow=16, ncol=15)                 # nrow and ncol need to have the same dimensions as the dataframe "results[1]"
decomp_ages <- c(decomp_sub_data[nrow(decomp_sub_data),6],rev(decomp_sub_data[,5]))
decomp_depths <- c(0,decompacted_thickness[nrow(decompacted_thickness),])
yw <- unlist(results[6])
yw_depths <- c(0,yw)


#pdf("SubsidencePlot.pdf")          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
plot(decomp_ages,yw_depths, type="o", col="blue",                                                       # Plot type and colour
     ylim=rev(range(decomp_depths)),
     xlim=rev(range(decomp_ages)),
     xlab="Age [Ma]", ylab="Depth [m]",                                                                 # Axis labels
     main="Subsidence Profiles of the MIQRAT-1 Well")                                                   # Title of plot
lines(decomp_ages,decomp_depths, type="o", col="red")
legend("bottomleft", legend=c("Water-Loaded Tectonic Subsidence","Total Subsidence"),                   # Legend titles
       col=c("blue","red"), lty=c(1,1))                                                                 # Colour and icon type for legend
#dev.off()                          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****



##-------------------------------------------------------------BAYESIAN THERMAL SUBSIDENCE MODEL----------------------------------------------------------------##



# Redo the backstripping, but only for the portion of the data suspected of undergoing thermal subsidence
numDataPoints <- 4                                            # number of points used for fitting the thermal subsidence curve
decomp_sub_data1 <- tail(decomp_sub_data,n=numDataPoints-1)   # n value starting from the base of the column
results1 <- backstripping(decomp_sub_data1)

## Determine the stretching factor using the inversion method (Allen and Allen, 2013)
decomp_ages1 <- c(decomp_sub_data1[nrow(decomp_sub_data1),6],rev(decomp_sub_data1[,5]))
yw1 <- unlist(results1[6])
beta1 <- beta_factor(decomp_ages1,yw1)
B1 <- unlist(beta1[2])

# Display the stretching factor
B1




## We will now attempt to better quantify the tectonic subsidence history by doing a Bayesian curve fit to the water-loaded subsidence data.
## The parameter of interest here is the beta factor.

# First, we'll need to define a function that calculates the likelihood of each yw data point based on the predictor
likelihood_yw <- function(params, duration, tsub_data, index) {
    B <- params[1]
    sd <- params[2]
    riftDriftOffset <- params[3]
    index <- index + 1 + riftDriftOffset
    tsub <- thermal_subsidence(B,duration)
    pred <- tsub[2,index]

    singlelikelihoods = dnorm(tsub_data[2,], mean = pred, sd = sd, log = T)
    sumll = sum(singlelikelihoods)
    return(sumll)
}



# Here, we shall define our prior function. We know from Le Guerroue et al (2006) that B must be less than 1.41
prior_yw <- function(params){
    B <- params[1]
    sd <- params[2]
    riftDriftOffset <- params[3]

    if(B >= 1.12 & B <= 1.41 & riftDriftOffset >= 0) {
        Bprior <- dnorm(B, mean = B1, sd = 0.02, log = T)                                       # sd can be adjusted accordingly to properly explore parameter space
        sdPrior <- dunif(sd, min = 0, max = 150, log = T)                                        # min and max values can be adjusted accordingly
        riftDriftOffsetPrior <-dnorm(riftDriftOffset, mean=startAgeOffset, sd = 0.5, log =T)    # startAgeOffset is always the max uncertainty value for the posterior age extrations. This places the start age at estimate for the base of the Hadash Fm.
        return(Bprior + sdPrior + riftDriftOffsetPrior)
    } else {
        return(log(0))
    }
}

posterior_yw <- function(params, duration, tsub_data, index){
   return (likelihood_yw(params, duration, tsub_data, index) + prior_yw(params))
}

proposal_yw <- function(params) {
    return(rnorm(3, mean = params, sd = c(0.01,5,0.5)))                             # sd can be adjusted accordingly
}

metropolis_yw <- function(startvalues, duration, tsub_data, index, iterations){
    chain = array(dim = c(iterations+1,3))
    chain[1,] = startvalues
    for (i in seq(1,iterations)){
        proposal = proposal_yw(chain[i,])

        probab = exp(posterior_yw(proposal, duration, tsub_data, index) - posterior_yw(chain[i,], duration, tsub_data, index))
        if (is.nan(probab)) {
            chain[i+1,] = chain[i,]
        } else if (runif(1) < probab){
            chain[i+1,] = proposal
        } else {
            chain[i+1,] = chain[i,]
        }
    }
    return(chain)
}



# Prepare data for MCMC
agePoints <- decomp_ages[1] - decomp_ages
tsub_data <- rbind(round(agePoints[1:numDataPoints]),yw_depths[1:numDataPoints])  # The index numbers refer to the points in which thermal subsidence is occuring
index <- tsub_data[1,]                                                            # These are the values in the first row of tsub_data
duration <- range(round(decomp_ages))[2] - range(round(decomp_ages))[1] + 20      # This value must be greater than the total duration estimate of thermal subsidence
startAgeOffset <- extracted_ages[numAges,4]
startAge <- age_predictions$HDI[numAges,4]


# MCMC
startvalues = c(B1,50,startAgeOffset)
chain_yw = metropolis_yw(startvalues, duration, tsub_data, index, 100000)


# Plot parameter histograms and markov chains
burnIn = 10000

pdf("ParameterHistAndChains_n6.pdf")          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(2,3))
hist(chain_yw[-(1:burnIn),1],nclass=30, , main="Stretching Factor", xlab="Stretching Factor", ylab=NULL, yaxt="n")
mtext(paste(round(median(chain_yw[-(1:burnIn),1]),2)," +",round(HPDI(chain_yw[-(1:burnIn),1],0.95)[2]-median(chain_yw[-(1:burnIn),1]),2),"/ -", round(median(chain_yw[-(1:burnIn),1])-HPDI(chain_yw[-(1:burnIn),1],0.95)[1],2)),
    side=3, at=median(chain_yw[-(1:burnIn),1]), cex=0.65, col="red")
abline(v = median(chain_yw[-(1:burnIn),1]), col="red", lwd=2)
abline(v = HPDI(chain_yw[-(1:burnIn),1],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(chain_yw[-(1:burnIn),1],0.95)[1], col="red", lwd=2, lty=2)
hist(startAge-chain_yw[-(1:burnIn),3],nclass=30, main="Rift-Drift Transition", xlab="Year [Ma]", ylab=NULL, yaxt="n")
mtext(paste(round(median(startAge-chain_yw[-(1:burnIn),3]),2)," +",round(median(chain_yw[-(1:burnIn),3])-HPDI(chain_yw[-(1:burnIn),3],0.95)[1],2),"/ -", round(HPDI(chain_yw[-(1:burnIn),3],0.95)[2]-median(chain_yw[-(1:burnIn),3]),2)),
    side=3, at=median(startAge-chain_yw[-(1:burnIn),3]), cex=0.65, col="red")
abline(v = median(startAge-chain_yw[-(1:burnIn),3]), col="red", lwd=2)
abline(v = HPDI(startAge-chain_yw[-(1:burnIn),3],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(startAge-chain_yw[-(1:burnIn),3],0.95)[1], col="red", lwd=2, lty=2)
hist(chain_yw[-(1:burnIn),2],nclass=30, main="sd", xlab="sd", ylab=NULL, yaxt="n")
mtext(paste(round(median(chain_yw[-(1:burnIn),2]),2)," +",round(HPDI(chain_yw[-(1:burnIn),2],0.95)[2]-median(chain_yw[-(1:burnIn),2]),2),"/ -", round(median(chain_yw[-(1:burnIn),2])-HPDI(chain_yw[-(1:burnIn),2],0.95)[1],2)),
    side=3, at=median(chain_yw[-(1:burnIn),2]), cex=0.65, col="red")
abline(v = median(chain_yw[-(1:burnIn),2]), col="red", lwd=2)
abline(v = HPDI(chain_yw[-(1:burnIn),2],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(chain_yw[-(1:burnIn),2],0.95)[1], col="red", lwd=2, lty=2)
plot(chain_yw[-(1:burnIn),1], type = "l", xlab="Markov Chain Number", ylab="Stretching Factor", main = "Stretching Factor")
mtext(paste(round(median(chain_yw[-(1:burnIn),1]),2)," +",round(HPDI(chain_yw[-(1:burnIn),1],0.95)[2]-median(chain_yw[-(1:burnIn),1]),2),"/ -", round(median(chain_yw[-(1:burnIn),1])-HPDI(chain_yw[-(1:burnIn),1],0.95)[1],2)),
    side=3, cex=0.65, col="red")
abline(h = median(chain_yw[-(1:burnIn),1]), col="red", lwd=2)
abline(h = HPDI(chain_yw[-(1:burnIn),1],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(chain_yw[-(1:burnIn),1],0.95)[1], col="red", lwd=2, lty=2)
plot(startAge-chain_yw[-(1:burnIn),3], type = "l", xlab="Markov Chain Number" , ylab="Year [Ma]", main = "Rift-Drift Transition", )
mtext(paste(round(median(startAge-chain_yw[-(1:burnIn),3]),2)," +",round(median(chain_yw[-(1:burnIn),3])-HPDI(chain_yw[-(1:burnIn),3],0.95)[1],2),"/ -", round(HPDI(chain_yw[-(1:burnIn),3],0.95)[2]-median(chain_yw[-(1:burnIn),3]),2)),
    side=3, cex=0.65, col="red")
abline(h = median(startAge-chain_yw[-(1:burnIn),3]), col="red", lwd=2)
abline(h = HPDI(startAge-chain_yw[-(1:burnIn),3],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(startAge-chain_yw[-(1:burnIn),3],0.95)[1], col="red", lwd=2, lty=2)
plot(chain_yw[-(1:burnIn),2], type = "l", xlab="Markov Chain Number" , ylab="sd", main = "sd", )
mtext(paste(round(median(chain_yw[-(1:burnIn),2]),2)," +",round(HPDI(chain_yw[-(1:burnIn),2],0.95)[2]-median(chain_yw[-(1:burnIn),2]),2),"/ -", round(median(chain_yw[-(1:burnIn),2])-HPDI(chain_yw[-(1:burnIn),2],0.95)[1],2)),
    side=3, cex=0.65, col="red")
abline(h = median(chain_yw[-(1:burnIn),2]), col="red", lwd=2)
abline(h = HPDI(chain_yw[-(1:burnIn),2],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(chain_yw[-(1:burnIn),2],0.95)[1], col="red", lwd=2, lty=2)
dev.off()                                  # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****



# Plot the Bayesian thermal subsidence model
chain_yw <- chain_yw[-(1:burnIn),]
t_subsidence_HPDI <- matrix(0,2,duration+1)
t_ages_HPDI <- matrix(0,2,duration+1)

t_subsidence_med <- thermal_subsidence(median(chain_yw[,1]),duration)
t_subsidence_HPDI[1,] <- thermal_subsidence(HPDI(chain_yw[,1],0.95)[2],duration)[2,]
t_subsidence_HPDI[2,] <- thermal_subsidence(HPDI(chain_yw[,1],0.95)[1],duration)[2,]
t_ages_med <- startAge - median(chain_yw[,3]) - t_subsidence_med[1,]
t_ages_HPDI[1,] <- startAge - HPDI(chain_yw[,3],0.95)[1] - t_subsidence_med[1,]
t_ages_HPDI[2,] <- startAge - HPDI(chain_yw[,3],0.95)[2] - t_subsidence_med[1,]

shadeRange <- t_ages_HPDI[1,1]-t_ages_HPDI[2,1]


#pdf("ThermalSubsidencePlot_n6.pdf")          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
plot(NULL,
     ylim=rev(range(yw_depths)),
     xlim=c(640,550),
     xlab="Age [Ma]", ylab="Depth [m]",                                                                                                         # Axis labels
     main="Thermal Subsidence Curve")                                                                                                           # Title of plot
for(i in seq(0,shadeRange,by=0.1)) {
    shade(t_subsidence_HPDI, t_ages_HPDI[1,]-i, col="gray")
}
shade(t_subsidence_HPDI, t_ages_med, col="gray")
shade(t_subsidence_HPDI, t_ages_HPDI[1,], col="gray")
arrows(x0=decompAgesMax, x1=decompAgesMin ,y0=yw_depths, y1=yw_depths, code=3, angle=90, length=0.05, col="blue", lwd=0.5)
lines(t_ages_med,t_subsidence_med[2,],type="l",col="black", lwd=2)
points(decomp_ages,yw_depths, type="o", col="blue", lwd=1.25)
legend("topright", legend=c("Water-Loaded Tectonic Subsidence","Thermal Subsidence Curve [Median]", "Thermal Subsidence Curve [95% HPDI]"),     # Legend titles
       col=c("blue","black","gray"), pch=c(1,NA,15), lty=c(NA,1,NA))
#dev.off()                          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****






## Create composite of n=3 and n=4 configurations

#chain_yw3 <- chain_yw[-(1:burnIn),]                    # These need to be un-commented and saved after each configuration run
#chain_yw4 <- chain_yw[-(1:burnIn),]                    # These need to be un-commented and saved after each configuration run

chain_ywComposite <- rbind(chain_yw3,chain_yw4)



# Histograms and Chain Plots of composite

#pdf("ParameterHistAndChains_Composite.pdf")          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
par(mfrow = c(2,3))
hist(chain_ywComposite[,1],nclass=30, , main="Stretching Factor", xlab="Stretching Factor", ylab=NULL, yaxt="n")
mtext(paste(round(median(chain_ywComposite[,1]),2)," +",round(HPDI(chain_ywComposite[,1],0.95)[2]-median(chain_ywComposite[,1]),2),"/ -", round(median(chain_ywComposite[,1])-HPDI(chain_ywComposite[,1],0.95)[1],2)),
    side=3, at=median(chain_ywComposite[,1]), cex=0.65, col="red")
abline(v = median(chain_ywComposite[,1]), col="red", lwd=2)
abline(v = HPDI(chain_ywComposite[,1],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(chain_ywComposite[,1],0.95)[1], col="red", lwd=2, lty=2)
hist(startAge-chain_ywComposite[,3],nclass=30, main="Rift-Drift Transition", xlab="Year [Ma]", ylab=NULL, yaxt="n")
mtext(paste(round(median(startAge-chain_ywComposite[,3]),2)," +",round(median(chain_ywComposite[,3])-HPDI(chain_ywComposite[,3],0.95)[1],2),"/ -", round(HPDI(chain_ywComposite[,3],0.95)[2]-median(chain_ywComposite[,3]),2)),
    side=3, at=median(startAge-chain_ywComposite[,3]), cex=0.65, col="red")
abline(v = median(startAge-chain_ywComposite[,3]), col="red", lwd=2)
abline(v = HPDI(startAge-chain_ywComposite[,3],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(startAge-chain_ywComposite[,3],0.95)[1], col="red", lwd=2, lty=2)
hist(chain_ywComposite[,2],nclass=30, main="sd", xlab="sd", ylab=NULL, yaxt="n")
mtext(paste(round(median(chain_ywComposite[,2]),2)," +",round(HPDI(chain_ywComposite[,2],0.95)[2]-median(chain_ywComposite[,2]),2),"/ -", round(median(chain_ywComposite[,2])-HPDI(chain_ywComposite[,2],0.95)[1],2)),
    side=3, at=median(chain_ywComposite[,2]), cex=0.65, col="red")
abline(v = median(chain_ywComposite[,2]), col="red", lwd=2)
abline(v = HPDI(chain_ywComposite[,2],0.95)[2], col="red", lwd=2, lty=2)
abline(v = HPDI(chain_ywComposite[,2],0.95)[1], col="red", lwd=2, lty=2)
plot(chain_ywComposite[,1], type = "l", xlab="Markov Chain Number", ylab="Stretching Factor", main = "Stretching Factor")
mtext(paste(round(median(chain_ywComposite[,1]),2)," +",round(HPDI(chain_ywComposite[,1],0.95)[2]-median(chain_ywComposite[,1]),2),"/ -", round(median(chain_ywComposite[,1])-HPDI(chain_ywComposite[,1],0.95)[1],2)),
    side=3, cex=0.65, col="red")
abline(h = median(chain_ywComposite[,1]), col="red", lwd=2)
abline(h = HPDI(chain_ywComposite[,1],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(chain_ywComposite[,1],0.95)[1], col="red", lwd=2, lty=2)
plot(startAge-chain_ywComposite[,3], type = "l", xlab="Markov Chain Number" , ylab="Year [Ma]", main = "Rift-Drift Transition", )
mtext(paste(round(median(startAge-chain_ywComposite[,3]),2)," +",round(median(chain_ywComposite[,3])-HPDI(chain_ywComposite[,3],0.95)[1],2),"/ -", round(HPDI(chain_ywComposite[,3],0.95)[2]-median(chain_ywComposite[,3]),2)),
    side=3, cex=0.65, col="red")
abline(h = median(startAge-chain_ywComposite[,3]), col="red", lwd=2)
abline(h = HPDI(startAge-chain_ywComposite[,3],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(startAge-chain_ywComposite[,3],0.95)[1], col="red", lwd=2, lty=2)
plot(chain_ywComposite[,2], type = "l", xlab="Markov Chain Number" , ylab="sd", main = "sd", )
mtext(paste(round(median(chain_ywComposite[,2]),2)," +",round(HPDI(chain_ywComposite[,2],0.95)[2]-median(chain_ywComposite[,2]),2),"/ -", round(median(chain_ywComposite[,2])-HPDI(chain_ywComposite[,2],0.95)[1],2)),
    side=3, cex=0.65, col="red")
abline(h = median(chain_ywComposite[,2]), col="red", lwd=2)
abline(h = HPDI(chain_ywComposite[,2],0.95)[2], col="red", lwd=2, lty=2)
abline(h = HPDI(chain_ywComposite[,2],0.95)[1], col="red", lwd=2, lty=2)
#dev.off()



# Plot composite
t_subsidence_HPDI <- matrix(0,2,duration+1)
t_ages_HPDI <- matrix(0,2,duration+1)

t_subsidence_med <- thermal_subsidence(median(chain_ywComposite[,1]),duration)
t_subsidence_HPDI[1,] <- thermal_subsidence(HPDI(chain_ywComposite[,1],0.95)[2],duration)[2,]
t_subsidence_HPDI[2,] <- thermal_subsidence(HPDI(chain_ywComposite[,1],0.95)[1],duration)[2,]
t_ages_med <- startAge - median(chain_ywComposite[,3]) - t_subsidence_med[1,]
t_ages_HPDI[1,] <- startAge - HPDI(chain_ywComposite[,3],0.95)[1] - t_subsidence_med[1,]
t_ages_HPDI[2,] <- startAge - HPDI(chain_ywComposite[,3],0.95)[2] - t_subsidence_med[1,]

shadeRange <- t_ages_HPDI[1,1]-t_ages_HPDI[2,1]


#pdf("ThermalSubsidencePlot_Composite.pdf")          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
plot(NULL,
     ylim=rev(range(yw_depths)),
     xlim=c(640,550),
     xlab="Age [Ma]", ylab="Depth [m]",                                                                                                         # Axis labels
     main="Thermal Subsidence Curve")                                                                                                           # Title of plot
for(i in seq(0,shadeRange,by=0.1)) {
    shade(t_subsidence_HPDI, t_ages_HPDI[1,]-i, col="gray")
}
shade(t_subsidence_HPDI, t_ages_med, col="gray")
shade(t_subsidence_HPDI, t_ages_HPDI[2,], col="gray")
arrows(x0=decompAgesMax, x1=decompAgesMin ,y0=yw_depths, y1=yw_depths, code=3, angle=90, length=0.05, col="blue", lwd=0.5)
lines(t_ages_med,t_subsidence_med[2,],type="l",col="black", lwd=2)
points(decomp_ages,yw_depths, type="o", col="blue", lwd=1.25)
legend("topright", legend=c("Water-Loaded Tectonic Subsidence","Thermal Subsidence Curve [Median]", "Thermal Subsidence Curve [95% HPDI]"),     # Legend titles
       col=c("blue","black","gray"), pch=c(1,NA,15), lty=c(NA,1,NA))
#dev.off()                          # ***** UN-COMMENT THIS LINE TO SAVE PLOT *****
