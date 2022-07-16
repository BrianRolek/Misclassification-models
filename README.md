# Misclassification-models
Bayesian dynamic occupancy models that account for misclassifications (false positive detections) and detection probability (false negative detections) using R and the nimble package.

Supplemental materials for: 
Larned, A. F., B. W. Rolek, K. Silaphone, S. Pruett, R. Bowman, B. Lohr. 2022. Habitat use and dynamics of an Endangered passerine, the Florida Grasshopper Sparrow (Ammodramus savannarum floridanus), after prescribed fires.
 
Contact author: brianrolek at gmail.com.
 
Metadata associated with the file "data\data-final.Rdata". Once loaded into R, this file contains the object "dat.conv", a list containing data for input into the JAGS model. Data fields include: 

+ Y = Observed states: 1 = no detection, 2 = uncertain detection, 3 = certain detection. A matrix of dimensions nsite x nvisit x nyear.
+ nsite = The number of sites surveyed at Avon Park and retained for analyses.
+ nvisit = The maximum number of survey visits within breeding seasons.
+ nyear = The total number of years of surveys. 
+ date = The date during each survey visit. A matrix of dimensions nsite x nvisit x nyear.
+ hr = Time after civil dawn during the start of each survey visit. A matrix of dimensions nsite x nvisit x nyear.

