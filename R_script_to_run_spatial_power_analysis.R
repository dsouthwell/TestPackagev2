
###################################################################################################################
#SPATIAL POWER ANALYSIS MODELLING IN R (SPAMR)
###################################################################################################################

# Darren Southwell 
# Email: darren.southwell@unimelb.edu.au

# This package calculates the statistical power to detect simulated occupancy trends for multiple species in spatially explicit landscapes
# The package requires an occupancy raster layer for each species loaded into R as a raster stack. 
# Simulations also require detectability raster layers for each species. Up to four independent detection methods can be simulated for each species separately. 
# Users must specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to sites for each detection methods.
# Pre-determined monitoring sites can be loaded into the package for analysis. 
# Alternatively, monitoring sites can be selected randomly at the start of each simulation, be stratified across environmental layers, or be placed on cells with the highest relative species richness.
# Occupancy and detectability raster layers can remain static over time or be dynamic in response to deterministic or stochastic disturbance events.
# To model a deterministic disturbance, users must load in a time series of raster layers, specifying where the disturbance will occur and its proportional effect on occupancy
# To model a stochastic disturbance, users must load in a history of disturbances, a function specifying the relationship between a disturbance event and time since a disturbance, and functions relating occupancy and detectability with time since a disturbance for each species. 
# Users simulate either an increasing or decreasing trend in occupancy, the magnitude of which is given by the effect size.
# Simulations can be run in series or in parallel. Users choose the alpha significance level and either a one-tailed or two-tailed significance test
# After running the simulations, the package returns the proportion of times (i.e. statistical power) a significant trend in occupancy was detected in the simulated detection histories 
# Statistical power can be estimated acorss the entire landscape, and within smaller-level management units

###################################################################################################################
#STEP 1 - Load package and set working directory
###################################################################################################################

#rm(list=ls()) #Remove any pre-existing files in the global environment
#require(devtools)
#install_github("dsouthwell/TestPackagev2")
#require(TestPackage2) #Loading this package will automatically load the other required packages 
#setwd("~/Kakadu/Paper2/Stacks/") #Set the working directory were all files are stored

###################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species
###################################################################################################################

# Simulations require raster layers of occupancy and detectability to be loaded into R for each species. These raster layers need to be built beforehand and 
# loaded into R as raster stacks. All cell values should be between 0 and 1. Cells not considered for monitoring should be set to NA.
# Occupancy values can be constant or vary across space in response to environmental covariates
# The package contains sample occupancy raster layers for 10 species in Kakadu and Nitmiluk National Park, in northern Australia. To view these, type 'occ', or load in your own below.

#occ <- stack("Rep_occ_restricted.tif") #Load in raster stack of species occupancy raster layers 
#plot(occ[[1]]) #Plot occupancy stack

# Simulations also require detectability raster layers for each species. All cell values should be between 0 and 1. Detectability layers should have the same dimensions as the occupancy layers.   
# Detectability values can be constant or vary across space in response to environmental covariates.
# Any combination of 4 separate detection methods are allowed for each species. For example, species A can be detected with methods 1 and 2, while species B can be detected with methods 3 and 4, and so fourth.
# Importantly, 4 raster stacks (one for each detection method) must be loaded in regardless of whether they are all used. 
# If a particular method does not apply to a species, all cells in the correspnding raster layer should be equal to zero. 
# The package contains a sample detectability raster layers for 10 species in Kakadu and Nitmiluk National Parks. Type det1.method1 to view the first detection method or load in your own below.

#det.method1 <- stack("Rep_search.tif") #Load in species detectability layers for method 1 
#det.method2 <- stack("Rep_spot.tif") #Load in species detectability layers for method 2 
#det.method3 <- stack("Rep_spot.tif") #Load in species detectability layers for method 3 
#det.method4 <- stack("Rep_spot.tif") #Load in species detectability layers for method 4 

###################################################################################################################
#STEP 3 - Load a table specifying which detection method applies to each species
###################################################################################################################

# Simulations require that the user loads in a table (species.list) specifying which detection method is relevant to each species.
# The rows in this table lists each species, the columns refer to each of the four detection methods. A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
# If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. Otheriwse, load in your own below from a csv file.

#species.list <- read.csv("E:/Kakadu_temp_files/Species_list_package.csv",stringsAsFactors=FALSE) #Load in a table specifying which methods are relevant to each species

###################################################################################################################
#STEP 4 - Decide whether or not to model deterministic or stochastic disturbances 
###################################################################################################################

# The package can simulate a static landscape over time or a landscape that is effected by disturbance events
# To model disturbance events, the user specifies model.disturbance <- TRUE, otherwise they specify model.disturbance <- FALSE 

model.disturbance <- TRUE #Choose whether disturbances are simulated at monitoring sites

# Two types of disturbance events can be simulated - 1) Determistic disturbance events or 2) Stochastic disturbance events.
# For deterministic disturbances, the location of the disturbance and its effect on occupancy is known beforehand.
# For stoachstic disturbances, the incidence of a disturbance is probabilistic, but its effect on occupancy and detectability is known.

disturbance.type <- "stochastic" #Specify a 'stochastic' or 'deterministic' disturbance event

# For deterministic disturbances, a raster stack must be loaded in, with a raster layer for each year of the simulation period specifyig the location and effect of the disturbance.
# Cells not subject to the disturbance must have values of 1, while cells disturbed should contain a value specifying the proportional change in occupancy due to the disturbance
# A disturbance layer should be loaded in regardless of whether model.disturbance is set to TRUE. If the user does not want to model a deterministic disturbance, all cells should be set to 1. 
# Load in a disturbance raster stack or type 'disturbance' to look at an example.

#disturbance <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in deterministic disturbance event

# For stochastic disturbances, the user must load in a raster stack of the disturbance history (fire.hist). 
# The disturbance history stack should contain a raster layer for each proceeding year. Cells should indicate whether a disturbance occurred (1) or not (0).
# The package then uses this stack to calculate the 'time since a disturbance' and 'number of disturbance' at each cell for the proceeding period
# If the user is not modelling a stochastic disturbance, then a dummy disturbance history layer should be loaded anyway. 
# Load in the disturbance history raster stack or type 'fire.hist' to look at an example.

#fire.hist <- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") #Load in disturbance history raster stack

# For stochastic disturbances, the user must also load a raster stack of the environmental covariates that influence occupancy and detectability (covariates).  
# The covariates stack is used to update occupancy and detectability raster layers in response to stochastic disturbance events 
# There can be any number of layers in the covariates stack - it depends on what is thought to influence species occupancy and detectability
# If the user is not modelling a stochastic disturbance, then a dummy disturbance history layer should be loaded anyway.
# Load in the covariates stack below or type 'covariates' to look at an example. 
# In the example dataset, the raster layers correspond to: 1=Elevation, 2=Distance to creek, 3=Time since last fire, 4=Fire frequency, 5=Maximum temperature, 
# 6=Fire extent, 7=Fire patchiness, 8=Soil type, 9=Terrain ruggedness, 10=Annual rainfall, 11=Projected foliage cover

#covariates <- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") #Load in covariate raster stack

# For stochastic disturbances, a vector (pr.dist) must be specified relating the probability of a disturbance event with the time since the last disturbance.
# This function should be a vector with length equal to the number of layers in the fire.hist stack. All values should lie between 0 and 1.
# The first element specifies the probability of a disturbance event in one year since a disturbance, the second element is the probability of a disturbance event after 2 years without a disturbance, and so forth.
# An example function is specified below.

Pr.dist <- c(0.48,0.53,0.42,0.40,0.39,0.26,0.32,0.25,0.28,0.30,0.26,0.21,0.22,0.26,0.05) #Probability of a disturbance for each time period since the last disturbance

# For stochastic disturbances, users must specify how each species will respond to disturbances in terms of occupancy and detectability.
# Users must load in a table (species.occ) that defines how species occupancy is related to environmental covariates.
# In this table, the species are listed by rows, and the environmental covariates are listed by columns. An extra column is added for the intercept term.
# The elements of the species.occ table therefore specifies the species-specific relationships to environmental covariates 
# Note, the number of columns is equal to the number layers in the covariates raster layer (not inlcuding the intercept or species column). The order of covariates must also match the order of layers the covariate raster stack.
# If a particular covariate is not relevant to a species, then it should contain a 0 value. 
# Load in the species.occ table or type 'species.occ' to look at an example. 

#species.occ <- read.csv("E:/Kakadu_temp_files/Species_covariates_package_OccupancyModel.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

# The package must be told which layers are dynamic (i.e. are updated by simulating the stochastic disturbance events)
# Dist.occ.time identifies the 'time since disturbance' layer in the covariates raster stack. 
# Dist.occ.freq identifies the 'number of disturbances' layer in the covariates raster stack.
# Specify them below or enter the following for the example dataset.

dist.occ.time <- 3 #Points out that the time since distrubance layer is the 3rd layer in the covariates raster stack
dist.occ.freq <- 4 #Points out that the number of disturbance layer is the 4rd layer in the covariates raster stack

# If simulating stochastic disturbances, a table (species.det) must also be specified that relates detectability with environmental covariates for each species.  
# Note, this table has 4 intercept columns for the 4 possible detection methods.
# Once again, the number and order of environmental covariates should match the number of layers in the 'covariates' raster stack 
# Load in the species.det table or type 'species.det' to look at an example. 

#species.det <- read.csv("E:/Kakadu_temp_files/Species_covariates_package_DetectionModel.csv", header=TRUE,stringsAsFactors=FALSE)

###################################################################################################################
#STEP 4 - Define HOW MANY sites to monitor
###################################################################################################################

# Users can either: 1) load in existing monitoring sites; 2) select sites at random; 3) stratify sites across environmental layers; 4) position sites on cells with the highest relative species richness
# If assessing existing monitoring sites, set load.sites <- TRUE, FALSE otherwise.

load.sites <- TRUE #Set to TRUE to monitor existing sites, FALSE otherwise. 

# If users do wish to assess existing monitoring sites, they must provide the XY coordinates of these sites in a csv file. 
# Importantly the X coordinate should have the heading POINT_X and the Y coordinate should have the heading POINT_Y
# If load.sites is TRUE, load in the XY-coordinates. Type 'sites' to see an example.

#sites<-read.csv("~/Kakadu/Sites/Sites_Three_Parks.csv",head=TRUE, sep=",") #If you specifying TRUE above, load in the XY coordinates of the pre-selected monitoring sites. Note, the column header should be POINT_X and POINT_Y. Type 'site' to see sample data

# Users might wish to monitor only a subset of the loaded sites, or alternatively, all of the loaded sites plus other sites selected randomly throughout the landscape
# Users should select all.loaded.sites <- TRUE to ONLY monitor the loaded sites, or all.loaded.sites <- FALSE if they wish to monitor more or less than what is specified. 

all.loaded.sites <- FALSE #Set to TRUE to monitor all loaded sites, otherwise FALSE to selected a subset at random during each simulation 

#Specify the number of sites to simulate monitoring. 
#Importantly, if load.sites <- TRUE and all.loaded.sites <- TRUE, then n.sites must equal the number of rows in the sites table. 

n.sites <- 100 #If you're monitoring a subset of the pre-selected sites, enter the number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. Any additional sites are selected randomly across the entire landscape

###################################################################################################################
#STEP 5 - Decide WHERE to monitor 
###################################################################################################################

# Instead of loading in existing sites, they can be stratified across environmental layers
# To do this, users must load in a raster layer (stratify) specifying the environmental layers in which to sample. 
# Each environmental layer should be identified with a different integer value.
# The package automatically distributes an equal number of sites randomly within each of the layers
# Load in a stratify raster layer or type 'stratify' to look at an example. 

#stratify <- raster("~/Kakadu/Paper2/Mapsv2/Range.tif") #Load in environmental strata

# Alternatively, users can randomly select a proportion of sites in 'remote' and 'non-remote' areas. 
# For example, a remote area might include all cells that exceed a certain distance from roads etc.
# In this raster layer, 1's must be used to indicate remote cells, 0s must indicate non-remote cells.
# Load in a remote layer or type 'remote' to look at an example...

#remote <- raster("~/Kakadu/Paper2/Mapsv2/remote.tif") #Load in remote layer

# Users can then define the proportion of sites that are randomly selected in remote or remote areas
# Users not wanted to sepearte sites in remote and non-remote areas can make all cells in the remote layer equal to 1, and set the ratio equal to 1. 

R <- 1 #the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas. 

# If users have selected load.sites <- FALSE, they must now specify how to position the new sites in the landscape during each simulation.
# There are 3 options
    # 1) "random" - n.sites are selected at random at the start of each simulation
    # 2) "stratified" - an equal number of sites are randomly selected within each layer in the 'stratified' raster layer
    # 3) "maxocc" - n.sites are positioned on the cells with the highest relative species richness; that is, the highest summed occupancy values

new.site.selection <- "random" #Define how new sites are selected

# The package doesn't plot the location of sites during each step of the simulation. 
# However, users can plot the site locations for the first simulation to check whether sites are positioned as expected.

plot.sites <- TRUE #Set to TRUE to plot monitoring sites at the start of the simulation

###################################################################################################################
#Step 6: Define WHEN to monitor
###################################################################################################################

# The length of the monitoring program must be specified.
# Note, it should match the number of layers in the 'disturbance' raster stack.

Tmax <- 10 #Define length of monitoring program 

# The years in which monitoring occurs must also be specified.
# It is assumed that all sites are monitored each survey year
# Note, surveys must occur on the last year; that is, be equal to Tmax.

s.years <- c(1,5,Tmax) #the years when monitoring occurs

# Specify the number of repeat visits at a site durign a survey year for each of the four detection methods
# The first element in the vector below (n.method) corresponds to to the number of repeat visits using detection method 1.
# The second element in the vector below (n.method) corresponds to to the number of repeat visits using detection method 2, and so forth.
# Note, the number of repeat visits can be equal to 0 or should exceed 1 (i.e. not be equal to 1). 

n.method <- c(3,3,3,5) #Number of repeat visits to sites corresponding to methods 1,2,3,4

###################################################################################################################
#Step 7: set up parameters for power analysis
###################################################################################################################

# Users must now decide whether calculate power just across the landscape or also within nested smaller-level management units

park.level <- TRUE #Set to TRUE to estimate power at a park level, FALSE to just estimate power at a landscape-level 

# If power is to be calculated within smaller management units (park.level == TRUE), users must idenitify where these parks are within the landscape. 
# A parks layer should be loaded into R as a raster layer. Cells within each park should be identified by a separate integer, starting at 1.  
# Load in a parks layer or type 'parks' for an example.

#parks <- raster("~/Kakadu/Paper2/Mapsv2/parks.tif") #Load in parks layer

# There is an option to model either a decreasing or increasing trend in occupancy through time. 
# Users should set trend <- "decreasing" to model a decreasing trend, or trend <- "increasing" to model an increasing trend 

trend <- "decreasing" #Choose from "increasing" or "decreasing"

# A vector effect.size must be defined, which specifies the effect sizes to be simulated. 
# All elements must be between 0 and 1. There is no restriction on the length of the vector (i.e. any number of effect sizes can be tested)  
# If trend == decreasing, the effect sizes should be interpreted as the proportional reduction in initial occupancy of a cell
# If trend == increasing, the effect sizes should be intepreted as the proportional increase in the difference between initial occupancy of a cell and 1.

effect.size <- c(0.1,0.4,0.6,0.9) #Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes

# Users must set the significance level of the statistical test. Choose from 3 options -0.01, 0.05 or 0.1.

alpha <- 0.05 #Set significance level 

# Users must decide on a one-tailed or two-tailed significance test. A one-tailed test should be selected if the direction of change in occupancy is known.

two.tailed <- TRUE #Set to TRUE for a two tailed test, FALSE for a one-tailed test

# Define the number of simulations from which to calculate statistical power. Power is the proportion of these simulations in which there is a significant trend in the trend parameter

nsims <- 5 #Define the number of simulations used to calculate power


###################################################################################################################
#Step 8: Create empty arrays 
###################################################################################################################

# The code below creates empty arrays to store the results of simulations. This is done outside of the loop to speed up computations.

n.species <- nlayers(occ) #Calculates the number of species
n.park <- unique(parks) #Calculates the number of parks
det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays

#source("C:/Darren/12_Postdoc/2_Kakadu/Package/9_Functions_test_package_4Methods_random_XYsites.R") #Make sure the most recent source file is loaded

###################################################################################################################
#Step 9: Select monitoring sites
###################################################################################################################

# This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
# Occupancy and detectability values are then extracted from the sites for further processing. 
# If all loaded.sites == TRUE, or if new.site.selection =='maxocc', the select.sites function is not called again 

xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites) #Define monitoring sites

park.ID <- parks[cellFromXY(parks, xy.sites)] #Extract parks values at monitoring sites
xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites 
for (ss in 1:n.species){
  if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} #Extract detectability values from sites using method 1
  if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} #Extract detectability values from sites using method 2
  if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} #Extract detectability values from sites using method 3
  if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
}

###################################################################################################################
#Step 10: Check input layers and tables and run simulations
###################################################################################################################

# The last step before running the power analysis is to check that all raster layers have the same dimensions.
# A warning message will be displayed if there is a mismatch between rasters.
# Users should run the function below. If no warnings are returned, the simulations are ready to run.

check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK 

# This type of analysis can take a long time to run on a desktop depending on the size of the raster layers, the umber of species, whetehr disturbances are modelled etc. 
# Simulations can be sped up by running them across multiple cores.
# If running in parallel, select the number of cores to run. Set to 1 if running in series 

n.cores <- detectCores()-1 #Define the number of cores to run in parallel. If you son't want to run in parallel, set to 1
cl <- makeCluster(n.cores) #initiate clusters
registerDoParallel(cl) #initiate clusters

# Now, run the simulations by copying this function into R.
# NOTE - TO RUN IN PARALLEL USING N.CORES, MAKE SURE ITS SET TO "%DOPAR%"
# IF RUNNING IN SERIES, SET TO %DO%"
# Also, you cannot track progress if running in parallel. If running in series, a text output will keep track of the simulations

start.time<-proc.time() #Record start time
pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("TestPackage2"))%dopar%{ #Run the power analysis. Set %dopar% to run in parallel, %do% otherwise. 
  sapply(effect.size, run.power, 
         nsims=nsims, 
         alpha=alpha, 
         Tmax=Tmax, 
         s.years=s.years, 
         trend=trend,
         sites=sites, 
         model.disturbance=model.disturbance, 
         disturbance.type=disturbance.type,
         species.list=species.list, 
         dist.occ.time=dist.occ.time,
         dist.occ.freq=dist.occ.freq,
         park.level=park.level, 
         xy.sites=xy.sites, 
         R=R,
         two.tailed=two.tailed, 
         plot.sites=plot.sites, 
         n.sites=n.sites, 
         n.species=n.species, 
         n.park=n.park, 
         all.loaded.sites=all.loaded.sites,
         load.sites=load.sites, 
         new.site.selection=new.site.selection,
         occ.time=occ.time, 
         det.method1.time=det.method1.time, 
         det.method2.time=det.method2.time, 
         det.method3.time=det.method3.time, 
         det.method4.time=det.method4.time,
         n.method=n.method,
         occ=occ,
         det.method1=det.method1,
         det.method2=det.method2,
         det.method3=det.method3,
         det.method4=det.method4,
         stratify=stratify,
         remote=remote,
         parks=parks,
         Pr.dist=Pr.dist,
         disturbance=disturbance)
}
stopCluster(cl)
time.elapsed<-proc.time()-start.time #Record end time
time.elapsed #Print elapsed time

###################################################################################################################
#Step 11: Plot results
###################################################################################################################

#After running the simulations, calculate power and plot with respect to the effect size

Results <- plot.power(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power 

###################################################################################################################
#Step 12: Save results
###################################################################################################################

setwd("~/Kakadu/Paper2/Resultsv3")
save(pwr, file = paste("Rep_No_Fire_increase",".RData", sep = ""))


