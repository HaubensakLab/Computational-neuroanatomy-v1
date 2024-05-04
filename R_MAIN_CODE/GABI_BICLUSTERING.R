# Copyright (C) 2021 VRVis.
# All rights reserved.
# Contact: VRVis Forschungs-GmbH (office@vrvis.at)
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. All advertising materials mentioning features or use of this software
#    must display the following acknowledgement:
#    This product includes software developed by the VRVis Forschungs-GmbH.
# 4. Neither the name of the VRVis Forschungs-GmbH nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#

# libFolder<-"/home/imp/florian.ganglberger/Rlib"
setwd("C:\\Users\\ganglberger\\Documents\\Evolution\\evo_paper_v3_olga_data\\R_MAIN_CODE")
libFolder <- NULL #library folder of R


if(sum(!(c("GABi","foreach","doParallel","png","gplots","hash","matrixStats") %in% installed.packages(lib.loc=libFolder)[,"Package"]))>0){
  install.packages("GABi", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("foreach", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("doParallel", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("hash", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("gplots", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("png", repos = "http://cran.at.r-project.org/", lib=libFolder)
  install.packages("matrixStats", repos = "http://cran.at.r-project.org/", lib=libFolder)
}


library(GABi, lib.loc=libFolder)
library(gplots, lib.loc=libFolder) 
library(foreach, lib.loc=libFolder)
library(doParallel, lib.loc=libFolder)
library(matrixStats, lib.loc=libFolder)


# The function uses a genetic algorithm to find biclusters with constant columns. Every variable class is
# weighted equally, independend of its size. The consistency of the columns of a bicluster ("how constant they should be")
# is given by the variableClassConsistencyFactors for each variable class, given as the maximum standard deviation (consistency factor=0)
# or minimum mean (consistency factor=1), a bicluster column should have.  
# 
# dataset:  (rows=samples,column=variables), variableClasses (vector the size of dataset columns)
#             that represents the variable class by numbers (1 to amount of classes)
#
# variableClassConsistencyFactors: (vector the size of the amount of variable classes) that is the maximum standard 
#                                  deviation (consistency factor=0) or minimum mean (consistency factor=1) a bi-cluster column 
#                                  should have (depending on the variable distribution)
#
#                                  e.g. to have a real value column with the values of 0.1,0.13,0.1,0.17,0.2 as suitable 
#                                  bicluster (standard deviation is 0.0441) one needs to have a variableClassConsistencyFactor
#                                  bigger than 0.0441 such as 0.05 
#                                  e.g to have a binary column with the values of 1,1,1,0,1 as suitable bicluster (4/5 values are 1
#                                  --> 80% similar values) one needs to have a variableClassConsistencyFactor of 0.4
#                                  Here is a consistency factors cheat table FOR BINARY DATA
#                                            0.43 --> 75% similar binary values
#                                            0.40 --> 80% similar binary values
#                                            0.30 --> 90% similar binary values
#                                            0.21 --> 95% similar binary values
#                                            0.15 --> 97.5% similar binary values
#                                            0.1 -->  99% similar binary values
#
# variableClassConsistencyFactorsTypes: (a vector the size of the amount of variable classes) defines which type of consistency factor 
#                                       should be used for variableClassConsistencyFactors 
#                                       0 --> maximum standard deviation (bicluster columns will have a lower standard deviation than the consistency factor)
#                                       1 --> minimum mean (bicluster columns will have a higher mean than the consistency factor)
#                                       2 --> minimum value (bicluster columns will have higher values than the consistency factor)
#
# variableClassAllowZeroSum:  (vector the size of the amount of variable classes) which states if a bicluster 
#                             column can have a majority of zeroes (at least half). This is for the reason if someone includes a column with 
#                             a lot of zeros which would cluster together   
#
# variableClassNotOnTabuList: (vector the size of the amount of variable classes) which states if a bicluster 
#                             column should land on the tabu-list after it is part of a bicluster or not. If the
#                             variable/feature lands on the tabu-list, it won't count for fitness (leads to more diverse results)
#
# minimumAmountOfVariableClasses: A number. Minimum amount of variable classes that should be found in a bicluster
#
# amountOfBicluserSearches: a number which defines the amount of biclusters that the biclustering will be found (=k). GABi will be run 
#                           k-times, finding the best solution that has not yet been found for every run (note that a bicluster can have several
#                           variations, so it varies by a few samples - those will be also found and included! So the amount of biclusters found
#                           will be much more then k!)
# clusteringName: filename where the clustering saves intermediate results. If it is "" or NULL, no intermediate results will be saved
#
# useCores:   Amount of CPU cores that will be used for biclustering
# 
# nsubpops:   Amount of sub-populations for biclustering. Depends on the data. The smaller and less distinct the biclusters are, the more 
#             subpopulations are needed. Should be between the amount of CPU cores (since these are parallelized) and the amount of sample points
# 
# popsize:   Size of population for biclustering. Depends on the data. The more samples/rows your data has and less distinct the biclusters are, the larger the population is
#            needed.
#
# return value: A list of biclusters. Every element in the list is a bicluster, which consists of 3 variables: 
#                                                  samples: a vector with the (row) indizes of the samples
#                                                           that are part of the bicluster
#                                                  features: a vector with the (column) indices of the variables 
#                                                            that are part of the bicluster
#                                                  score: the score (how consistent is the bicluster). Higher = better
#
weighted_GABi <-
  function(dataSet,
           variableClasses,
           variableClassConsistencyFactors,
           variableWeightSimilar,
           variableClassAllowZeroSum,
           variableClassNotOnTabuList,
           variableClassConsistencyFactorsTypes,
           minimumAmountOfVariableClasses = 1,
           amountOfBicluserSearches = 5,
           clusteringName = "",
           useCores = 1,
           popsize,
           nsubpops) {
    #Converting parameters to work with GABis fitness function arguments
    fitnessArgs <- list()
    fitnessArgs$weightingsVector <-
      rep(1, length(variableClasses))
    fitnessArgs$consistency <- rep(0, length(variableClasses))
    fitnessArgs$allowZeroSum <- rep(0, length(variableClasses))
    fitnessArgs$variableClassNotOnTabuList <-
      rep(0, length(variableClasses))
    fitnessArgs$variableClasses <- variableClasses
    fitnessArgs$minimumAmountOfVariableClasses <-
      minimumAmountOfVariableClasses
    
    similarWeightings<- rep(1, length(variableClasses))
    
    for (i in unique(variableClasses)) {
      fitnessArgs$weightingsVector[variableClasses == i] <-
        length(variableClasses) / sum(variableClasses == i) / length(unique(variableClasses))
      fitnessArgs$consistency[variableClasses == i] <-
        variableClassConsistencyFactors[i]
      fitnessArgs$similarWeights[variableClasses == i] <-
        variableWeightSimilar[i]
      fitnessArgs$allowZeroSum[variableClasses == i] <-
        variableClassAllowZeroSum[i]
      fitnessArgs$consistencyType[variableClasses == i] <-
        variableClassConsistencyFactorsTypes[i]
      fitnessArgs$notOnTabuList[variableClasses == i] <-
        variableClassNotOnTabuList[i]
      
      if(variableWeightSimilar[i]==TRUE){
        if(fitnessArgs$allowZeroSum[i]==TRUE){
          similarWeightings[variableClasses == i] <- nrow(dataSet[,variableClasses == i]) / (colSums(t(t(dataSet[,variableClasses == i])>=0),na.rm=TRUE))
        }
        if(fitnessArgs$allowZeroSum[i]==FALSE){
          similarWeightings[variableClasses == i] <- nrow(dataSet[,variableClasses == i]) / (colSums(t(t(dataSet[,variableClasses == i])>0),na.rm=TRUE))
        }
      }
    }
    similarWeightings[similarWeightings > quantile(similarWeightings, 0.9)] <-
      quantile(similarWeightings, 0.9)
    sumOfSimilarWeightings <- similarWeightings/sum(similarWeightings)
    
    fitnessArgs$weightingsVector <-
      fitnessArgs$weightingsVector * sumOfSimilarWeightings * length(variableClasses)
    
    names(sumOfSimilarWeightings)<-colnames(dataSet)
    print(sumOfSimilarWeightings)
    names(fitnessArgs$weightingsVector)<-colnames(dataSet)
    print(fitnessArgs$weightingsVector)
    
    weighted.fitness <- function(chr) {
      # we don't use tabu$samples and tabu$features like in the standard fitness functions, since if we find a bicluster with 2 features, but
      # nearly all samples, all those samples can not be used anymore
      #
      # e.g. find this bicluster (1)
      #    0 0 0 0 0 0
      #    0 1 1 0 0 0
      #    0 1 1 0 2 2
      #    0 1 1 0 2 2
      #    0 1 1 0 0 0
      #    0 0 0 0 0 0
      # then this space is blocked for other solutions (2)
      #    0 X X 0 0 0
      #    X X X X X X
      #    X X X X X X
      #    X X X X X X
      #    X X X X X X
      #    0 X X 0 0 0
      # this is blocked if we check the identity for samples and features separately
      #    0 0 0 0 0 0
      #    0 X X 0 0 0
      #    0 X X 0 2 2
      #    0 X X 0 2 2
      #    0 X X 0 0 0
      #    0 0 0 0 0 0
      # therefore 2 can still be found!
      #
      #
      #
      #
      #
      
      
      score <- 0
      if (sum(chr == as.raw(1)) > 1) {
        featuresBool <-
          rep(FALSE, nrow(x)) #vector of which features meet the consistency criterium
        
        #for all columns with a consistency type = 0 (maximum standard deviation of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 0) > 0)
          featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                   0] <-
            featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] == 0] |
            (matrixStats::rowSds(x[diffRows, chr == as.raw(1)][fitnessArgs$consistencyType[diffRows] ==
                                                                 0,], na.rm = TRUE) <= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                           0])
        #for all columns with a consistency type = 1 (minimum mean of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 1) > 0)
          featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                   1] <-
            featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] == 1] |
            (rowMeans(x[diffRows, chr == as.raw(1)][fitnessArgs$consistencyType[diffRows] ==
                                                      1,], na.rm = TRUE) >= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                1])
        #for all columns with a consistency type = 2 (minimum value of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 2) > 0)
          featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                   2] <-
            featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] == 2] |
            (matrixStats::rowMins(x[diffRows, chr == as.raw(1)][fitnessArgs$consistencyType[diffRows] ==
                                                                  2,], na.rm = TRUE) >= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                            2])
        #for all columns with a consistency type = 3 (maximum value of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 3) > 0)
          featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                   3] <-
            featuresBool[diffRows][fitnessArgs$consistencyType[diffRows] == 3] |
            (matrixStats::rowMaxs(x[diffRows, chr == as.raw(1)][fitnessArgs$consistencyType[diffRows] ==
                                                                  3,], na.rm = TRUE) <= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                            3])
        
        featuresBool[is.na(featuresBool)] <-
          FALSE #if a feature has only NAs, than it is clearly not better than the consistency value
        
        if (sum(fitnessArgs$allowZeroSum == FALSE &
                featuresBool) > 1) {
          featuresBool[fitnessArgs$allowZeroSum == FALSE &
                         featuresBool] <-
            rowSums(x[fitnessArgs$allowZeroSum == FALSE &
                        featuresBool, chr == as.raw(1)] == 0, na.rm = TRUE) == 0.0 #no zeroes per feature allowed
        } else{
          if (sum(fitnessArgs$allowZeroSum == FALSE & featuresBool) > 0) {
            featuresBool[fitnessArgs$allowZeroSum == FALSE &
                           featuresBool] <-
              sum(x[fitnessArgs$allowZeroSum == FALSE &
                      featuresBool, chr == as.raw(1)] == 0, na.rm = TRUE) == 0.0 #no zeroes per feature allowed
          }
        }
        
        if (sum(featuresBool) > 1) {
          featuresBool[featuresBool] <-
            rowSums(is.na(x[featuresBool, chr == as.raw(1)]), na.rm = TRUE) == 0.0 #no NAs per feature allowed --> they can't be allowed in general!
        } else{
          if (sum(featuresBool) > 0) {
            featuresBool[featuresBool] <-
              sum(is.na(x[featuresBool, chr == as.raw(1)]), na.rm = TRUE) == 0.0 #no NAs per feature allowed --> they can't be allowed in general!
          }
        }
        
        if (sum(featuresBool) > 0) {
          score <-
            (
              length(setdiff(
                which(chr == as.raw(1)), fitnessArgs$tabu$samples
              )) *
                sum(fitnessArgs$weightingsVector[featuresBool &
                                                   !(
                                                     is.element(1:nrow(x), fitnessArgs$tabu$features) &
                                                       !fitnessArgs$notOnTabuList
                                                   )]) *
                (
                  length(unique(fitnessArgs$variableClasses[featuresBool &
                                                              !(
                                                                is.element(1:nrow(x), fitnessArgs$tabu$features) &
                                                                  !fitnessArgs$notOnTabuList
                                                              )])) >= fitnessArgs$minimumAmountOfVariableClasses
                )
            )
        }
      }
      score
    }
    
    #this is called at the end of a loop
    weighted.selection <- function(cols) {
      features <- 0
      if (length(cols) > 1) {
        featuresBool <- rep(FALSE, nrow(x[diffRows, cols]))
        
        #for all columns with a consistency type = 0 (maximum standard deviation of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 0) > 0)
          featuresBool[fitnessArgs$consistencyType[diffRows] == 0] <-
            featuresBool[fitnessArgs$consistencyType[diffRows] == 0] |
            (
              matrixStats::rowSds(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] ==
                                                      0,], na.rm = TRUE) <= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                0] &
                (
                  fitnessArgs$allowZeroSum[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                       0] |
                    rowSums(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] == 0,] ==
                              0, na.rm = TRUE) == 0.0
                )
            )
        #for all columns with a consistency type = 1 (minimum mean of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 1) > 0)
          featuresBool[fitnessArgs$consistencyType[diffRows] == 1] <-
            featuresBool[fitnessArgs$consistencyType[diffRows] == 1] |
            (
              rowMeans(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] == 1,], na.rm =
                         TRUE) >= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                      1] &
                (
                  fitnessArgs$allowZeroSum[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                       1] |
                    rowSums(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] == 1,] ==
                              0, na.rm = TRUE) == 0.0
                )
            )
        #for all columns with a consistency type = 2 (minimum value of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 2) > 0)
          featuresBool[fitnessArgs$consistencyType[diffRows] == 2] <-
            featuresBool[fitnessArgs$consistencyType[diffRows] == 2] |
            (
              matrixStats::rowMins(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] ==
                                                       2,], na.rm = TRUE) >= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                 2] &
                (
                  fitnessArgs$allowZeroSum[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                       2] |
                    rowSums(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] == 2,] ==
                              0, na.rm = TRUE) == 0.0
                )
            )
        #for all columns with a consistency type = 3 (maximum value of feature/column)
        if (sum(fitnessArgs$consistencyType[diffRows] == 3) > 0)
          featuresBool[fitnessArgs$consistencyType[diffRows] == 3] <-
            featuresBool[fitnessArgs$consistencyType[diffRows] == 3] |
            (
              matrixStats::rowMaxs(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] ==
                                                       3,], na.rm = TRUE) <= fitnessArgs$consistency[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                                                                                 3] &
                (
                  fitnessArgs$allowZeroSum[diffRows][fitnessArgs$consistencyType[diffRows] ==
                                                       3] |
                    rowSums(x[diffRows, cols][fitnessArgs$consistencyType[diffRows] == 3,] ==
                              0, na.rm = TRUE) == 0.0
                )
            )
        
        featuresBool[is.na(featuresBool)] <-
          FALSE #if a feature has only NAs, than it is clearly not better than the consistency value
        
        if (sum(featuresBool) > 1) {
          featuresBool[featuresBool] <-
            rowSums(is.na(x[diffRows,][featuresBool, cols]), na.rm = TRUE) == 0.0 #no NAs per feature allowed --> they can't be allowed in general!
        } else{
          if (sum(featuresBool) > 0) {
            featuresBool[featuresBool] <-
              sum(is.na(x[diffRows,][featuresBool, cols]), na.rm = TRUE) == 0.0 #no NAs per feature allowed --> they can't be allowed in general!
          }
        }
        
        
        features <- diffRows[featuresBool]
        
      }
      features
    }
    
    start_time <- Sys.time()
    
    source("GABi_parallel.R")
    results <- GABi_parallel(
      t(dataSet),
      maxLoop = amountOfBicluserSearches,
      popsize = popsize,
      experiod = 10,
      nsubpops = nsubpops,
      maxNgens = 5000,
      convergenceGens = 100,
      mfreq = 1,
      xfreq = 0.9,
      identityThreshold = 1,
      useCores = useCores,
      clusteringName = clusteringName,
      libraryLocations = libFolder,
      verbose = FALSE,
      fitnessFun = weighted.fitness,
      featureSelFun = weighted.selection,
      diffThreshold = 1,
      fitnessArgs = fitnessArgs
    )
    
    # use this for GABi_fixed
    # source("GABi_fixedl.R")
    # results <- GABi_fixed(
    #   t(dataSet),
    #   maxLoop = amountOfBicluserSearches,
    #   nSols=Inf,
    #   popsize = popsize,
    #   experiod = 10,
    #   nsubpops = nsubpops,
    #   maxNgens = 5000, #if you run the testdata, you can set this to 500 to get faster results
    #   convergenceGens = 100,
    #   mfreq = 1,
    #   xfreq = 0.9,
    #   identityThreshold = 1,
    #   verbose = FALSE,
    #   fitnessFun = weighted.fitness,
    #   featureSelFun = weighted.selection,
    #   diffThreshold = 1,
    #   fitnessArgs = fitnessArgs
    # )
    
    print(paste(
      "Weighted GABi took",
      difftime(Sys.time(), start_time, units = "mins"),
      "mins for ",
      length(results),
      "biclusters"
    ))
    
    biclusterGroups <- list()
    
    while (length(biclusterGroups) != length(results)) {
      if (length(biclusterGroups) > 0) {
        results <- biclusterGroups
        biclusterGroups <- list()
      }
      alreadyIn <- c()
      for (bI1 in 1:length(results)) {
        if (sum(alreadyIn == bI1) == 0) {
          print(paste0(
            "Merging BiclusterGroup: ",
            length(biclusterGroups) + 1
          ))
          biclusterGroupSamples <- c(results[[bI1]]$samples)
          biclusterGroupFeatures <- c(results[[bI1]]$features)
          biclusterGroupScores <- c(results[[bI1]]$score)
          
          areSimilarTobI1 <- c(0)
          countSimilar <- 0
          while (length(areSimilarTobI1) > 0) {
            leftToCheck <- setdiff(1:length(results), alreadyIn)
            if (length(leftToCheck) > 0) {
              areSimilarTobI1 <- leftToCheck[(sapply(leftToCheck, function(bI2) {
                (length(
                  intersect(biclusterGroupSamples, results[[bI2]]$samples)
                ) / min(
                  length(biclusterGroupSamples),
                  length(results[[bI2]]$samples)
                )) > 0.5 &&
                  (length(
                    intersect(biclusterGroupFeatures, results[[bI2]]$features)
                  ) == length(biclusterGroupFeatures)) &&
                  (length(
                    intersect(biclusterGroupFeatures, results[[bI2]]$features)
                  ) == length(results[[bI2]]$features))
              }))]
              if (length(areSimilarTobI1) > 0) {
                biclusterGroupSamples <-
                  unique(c(biclusterGroupSamples, as.vector(unlist(
                    sapply(areSimilarTobI1, function(bI2) {
                      return(results[[bI2]]$samples)
                    })
                  ))))
                biclusterGroupFeatures <-
                  unique(c(biclusterGroupFeatures, as.vector(unlist(
                    sapply(areSimilarTobI1, function(bI2) {
                      return(results[[bI2]]$features)
                    })
                  ))))
                biclusterGroupScores <-
                  unique(c(biclusterGroupScores, as.vector(unlist(
                    sapply(areSimilarTobI1, function(bI2) {
                      return(results[[bI2]]$score)
                    })
                  ))))
                alreadyIn <- c(alreadyIn, areSimilarTobI1)
                countSimilar <-
                  countSimilar + length(areSimilarTobI1)
              }
            } else{
              areSimilarTobI1 <- c()
            }
          }
          #print(paste0(biclusterGroupFeatures))
          print(paste0(
            "Similar Clusters/Similar Scores: ",
            countSimilar,
            "/",
            sum(results[[bI1]]$score == unlist(sapply(1:length(results), function(bI2) {
              return(results[[bI2]]$score)
            })))
          ))
          
          
          biclusterGroups[[length(biclusterGroups) + 1]] <- list()
          biclusterGroups[[length(biclusterGroups)]]$samples <-
            sort(unique(biclusterGroupSamples))
          biclusterGroups[[length(biclusterGroups)]]$features <-
            sort(unique(biclusterGroupFeatures))
          biclusterGroups[[length(biclusterGroups)]]$score <-
            max(biclusterGroupScores)
          
        }
      }
    }
    
    
    print(paste("Returning ", length(biclusterGroups), "bicluster groups!"))
    return(biclusterGroups)
    
  }



################################################################
####################### Test data ##############################
################################################################


# Creates a testset with 150 samples and 50 variables (first 10 variables are 
# one class (simulated dnds) - values between 0 and 3, 11-50 the second class (simulated gene expression) - binary data)
set.seed(42)

dataSet<-matrix(rnorm(7500),150,50)/2.5
dataSet[,11:50]<-dataSet[,11:50]>0.5
dataSet[21:30,1:3]<-rnorm(30,0.3,0.03)
dataSet[21:30,4:6]<-rnorm(30,0.4,0.03)
dataSet[21:30,7:8]<-rnorm(20,0.7,0.03)
dataSet[21:30,9:10]<-rnorm(20,0.85,0.03)
dataSet[21:30,41:50]<-1

dataSet[101:120,6:8]<-rnorm(60,0.9,0.03)
dataSet[101:120,9:10]<-rnorm(20,0.8,0.03)
dataSet[101:120,21:30]<-runif(200,0,1)>0.01

dataSet[46:55,3:10]<-rnorm(80,0.65,0.03)
dataSet[46:55,11:25]<-runif(150,0,1)>0

dataSet[61:80,3:6]<-rnorm(80,0.5,0.02)
dataSet[61:80,31:35]<-runif(100,0,1)>0.01

dataSet[116:135,1:4]<-rnorm(80,0.9,0.03)
dataSet[116:135,36:45]<-runif(200,0,1)>0.01

dataSet[dataSet>1]<-1
dataSet[dataSet<0]<-0

#Add NAs for testing
# dataSet[11:40,1][sample(length(dataSet[11:40,1]),length(dataSet[11:40,1])*0.7)]<-NA
# dataSet[91:110,4:6][sample(length(dataSet[91:110,4:6]),length(dataSet[91:110,4:6])*0.6)]<-NA
# dataSet[106:125,3:4][sample(length(dataSet[106:125,3:4]),length(dataSet[106:125,3:4])*0.5)]<-NA
# dataSet[101:118,10]<-NA
# dataSet[150,1:10]<-NA
# dataSet[30,10]<-NA


#DNDS values are not normalized, but using the rank instead
dataSet[,1:10]<-(apply(as.matrix(dataSet[,1:10]),2,function(x){rank(x,na.last="keep")/sum(!is.na(x))}))

dataSet<-apply(dataSet,2,function(x){
  x[x==min(x)]<-0
  x
})


# Visualize data
heatmap.2(dataSet,scale="none",dendrogram="none",Rowv=FALSE,Colv=FALSE, trace="none",col=colorpanel(100, "white","black"))

# Gives the first 10 variables the class "1", 
# variables 11-50 class "2"
variableClasses<-c(rep(1,10),rep(2,40))

# Biclusters with columns of class "1" will have a maximum standard 
# deviation of 0.1, biclusters with columns of class "2" will 
# have a maximum standard deviation of 0.2. Since class 2 is binary, 
# this would be around 95% similar values (zeros or ones). 

variableClassConsistencyFactors<-c(0.05,0.9) #one needs to try out different consistency factors - FIND "OPTIMAL" VALUES
variableClassConsistencyFactorsTypes<-c(0,1) #0=standard deviation, 1=mean, 2=minimum or 3=maximum as consistency type

variableWeightSimilar<-c(TRUE,FALSE)

# We are not interested in clusters with zero dN/dS (class 1)
variableClassAllowZeroSum<-c(FALSE,TRUE)
# dN/dS transitions should be used for all biclustering loops, while gene-expression correlation (networks) not 
# this leads to findings of biclusters with more diverse features/networks, especially if there are a few dominating ones
variableClassNotOnTabuList<-c(TRUE,FALSE) 



#Execute biclustering
bicluster_results<-weighted_GABi(dataSet,
                                 variableClasses,
                                 variableClassConsistencyFactors,
                                 variableWeightSimilar,
                                 variableClassAllowZeroSum,
                                 variableClassNotOnTabuList,
                                 variableClassConsistencyFactorsTypes,
                                 minimumAmountOfVariableClasses=2,
                                 amountOfBicluserSearches=5,
                                 useCores = 1,
                                 popsize = 1024,
                                 nsubpops = 32)



print(paste0("Bicluster found: ",length(bicluster_results)))

# Visualize bicluster
dataSet_with_bicluster<-dataSet
colorSet<-colorpanel(30, "white","black")
for(biclusterIndex in 1:length(bicluster_results)){
  dataSet_with_bicluster[bicluster_results[[biclusterIndex]]$samples,bicluster_results[[biclusterIndex]]$features]<-3+biclusterIndex
  colorSet<-c(colorSet,rep(rainbow(length(bicluster_results))[biclusterIndex],10))
}

heatmap.2(dataSet_with_bicluster,scale="none",dendrogram="none",Rowv=FALSE,Colv=FALSE, trace="none",col=colorSet)

################################################################
####################### Real data ##############################
################################################################

prefix<-"NGT_SPLIT_TREE_9000_NAto0_"
dndsColsToUse<-c()
includeNAGenes<-TRUE

# Biclustering with real data. To achieve good results, one has to play with the variableClassConsistencyFactors variable and
# MAYBE (but not likely) with the variableClassAllowZeroSum and amountOfBicluserSearches
# Load data
print("Load DNDS")
dnds_per_gene<-read.csv2(paste0("storage//csvs//",prefix,"normalized_dnds_per_gene.csv"),row.names=1)

if(length(dndsColsToUse)==0){
  dndsColsToUse<-1:ncol(dnds_per_gene)
}

print(paste0("Used DNDS column indices: ",paste0(dndsColsToUse,collapse="_")))
print(paste0("Used DNDS column names: ",paste0(colnames(dnds_per_gene)[dndsColsToUse],collapse="_")))

print("Load gene wise correlation networks")
gene_wise_network_correlation<-read.csv2(paste0("storage//csvs//",prefix,"gene_task_network_correlation_brainwide.csv"),row.names=NULL)[,-1]

#Remove all genes that have in general a small correlation to the networks
filteredGenes<-apply(gene_wise_network_correlation,1,function(x){sum(x<0.1)})<dim(gene_wise_network_correlation)[2]
print(paste0("Filtered ",sum(filteredGenes==FALSE)," from ",nrow(gene_wise_network_correlation)," for small correlation"))  

if(includeNAGenes==FALSE){
  #Remove all genes with NA values
  filteredGenesByNA<-apply(dnds_per_gene,1,function(x){sum(is.na(x))})==0
  print(paste0("Filtered ",sum(filteredGenesByNA==FALSE)," from ",nrow(dnds_per_gene)," for NA"))
  
  filteredGenes<-filteredGenes & filteredGenesByNA
}

#filteredGenes[filteredGenes][sample(sum(filteredGenes),sum(filteredGenes)*0.75)]<-FALSE

print(paste0("Chromosome size: ",sum(filteredGenes))) 
gene_wise_network_correlationOrg<-gene_wise_network_correlation

#scale block wise, to treat fMRI activation maps and literature networks similar
gene_wise_network_correlation[filteredGenes,1:11]<-rank(gene_wise_network_correlation[filteredGenes,1:11],ties.method="max",na.last=FALSE)/length(unlist(gene_wise_network_correlation[filteredGenes,1:11]))
gene_wise_network_correlation[filteredGenes,12:22]<-rank(gene_wise_network_correlation[filteredGenes,12:22],ties.method="max",na.last=FALSE)/length(unlist(gene_wise_network_correlation[filteredGenes,12:22]))


#Rank normalization
gene_wise_network_correlation[filteredGenes,]<-t(apply(gene_wise_network_correlation[filteredGenes,],1,function(x){
  rank(x,ties.method="max",na.last=FALSE)/length(x)
}))

#Set smallest rank to zero (would be 1/#columns)
gene_wise_network_correlation[filteredGenes,]<-t(apply(gene_wise_network_correlation[filteredGenes,],1,function(x){
  x[x==min(x)]<-0
  x
}))

#set previously small correlations to 0
gene_wise_network_correlation[gene_wise_network_correlationOrg<=0.1]<-0
gene_wise_network_correlationOrg[gene_wise_network_correlationOrg<=0.1]<-0

gene_wise_network_correlation[!filteredGenes,]<-NA #set rest to NA, but it's not necessary since it is anyway excluded

fullColumnIndices<-c(dndsColsToUse,(1:ncol(gene_wise_network_correlation))+ncol(dnds_per_gene))

print("Define parameters")

#define parameters
variableClassConsistencyFactors<-c(0.9,0.75) # for contrast use c(0.1,0.75). #one needs to try out different consistency factors - FIND "OPTIMAL" VALUES
variableClasses<-c(rep(1,ncol(dnds_per_gene[,dndsColsToUse])),rep(2,ncol(gene_wise_network_correlation))) #the 2 variable classes - don't change this
variableClassConsistencyFactorsTypes<-c(2,2)  #0=standard deviation, 1=mean, 2=minimum or 3=maximum as consistency type. Here we use the minimum, so that the biclustedrs have at least the value of the consistency factor. For contrast use c(3,2), so that the gen expresssion maximum is 0.1 (given in the variableClassConsistencyFactors)
variableClassAllowZeroSum<-c(TRUE,FALSE) #GO terms is very sparse, therefore no "zero GO columns" in bicluster - no need to change this
variableClassNotOnTabuList<-c(TRUE,FALSE) #GO terms is very sparse, therefore no "zero GO columns" in bicluster - no need to change this

amountOfBicluserSearches<-5 #can be higher, but you should look at the bicluster scores - FIND A "SUITABLE" AMOUNT
minimumAmountOfVariableClasses<-2 #we want biclusters with at least 3 types of variables.
nsubpops<-80 #good size for this data set

#build dataSet
dataSet <-
  cbind(dnds_per_gene[, dndsColsToUse], gene_wise_network_correlation)
print(paste0("DNDS columns have ", paste0(sum(
  is.na(dnds_per_gene[, dndsColsToUse])
)), " NAs"))

#Filter genes that are anyway smaller than the minimum threshold (if consistency factor type is minimum (2))
for(parameterClass in 1:length(variableClassConsistencyFactors)){
  if(variableClassConsistencyFactorsTypes[parameterClass]==2){
    filteredGenes <- filteredGenes & apply(dataSet[,variableClasses==parameterClass], 1, function(x) {
      sum(x < variableClassConsistencyFactors[parameterClass], na.rm = TRUE) < sum(!is.na(x))
    })
    print(paste0("Filtered all genes that have all values paramter ",parameterClass," values below ",variableClassConsistencyFactors[parameterClass],", since minimum constraint applied, new chromosome size: ", sum(filteredGenes)))
  }
  if(variableClassConsistencyFactorsTypes[parameterClass]==3){
    filteredGenes <- filteredGenes & apply(dataSet[,variableClasses==parameterClass], 1, function(x) {
      sum(x > variableClassConsistencyFactors[parameterClass], na.rm = TRUE) < sum(!is.na(x))
    })
    print(paste0("Filtered all genes that have all values paramter ",parameterClass," values above ",variableClassConsistencyFactors[parameterClass],", since maximum constraint applied, new chromosome size: ", sum(filteredGenes)))
  }
}

#Execute biclustering

bicluster_results <-
  weighted_GABi(
    dataSet[filteredGenes, ],
    variableClasses,
    variableClassConsistencyFactors,
    variableWeightSimilar,
    variableClassAllowZeroSum,
    variableClassNotOnTabuList,
    variableClassConsistencyFactorsTypes,
    minimumAmountOfVariableClasses,
    amountOfBicluserSearches,
    useCores = 1,
    popsize = dim(dataSet[filteredGenes, ])[1] * 250,
    nsubpops = round(dim(dataSet[filteredGenes, ])[1] * 250 / 80)
  )


print(paste0("Bicluster found: ",length(bicluster_results)))
#adapt sample index for the filtered genes/samples/rows
if(length(bicluster_results)>0){ 
  for(bcIndex in 1:length(bicluster_results)){
    for(actFiltered in which(filteredGenes==FALSE)){
      bicluster_results[[bcIndex]]$samples[bicluster_results[[bcIndex]]$samples>=actFiltered]<-bicluster_results[[bcIndex]]$samples[bicluster_results[[bcIndex]]$samples>=actFiltered]+1
    }
  }
  print(paste0("Filtered ",sum(filteredGenes==FALSE)," genes of ",length(bicluster_results)," results!"))  
  for(bcIndex in 1:length(bicluster_results)){
    bicluster_results[[bcIndex]]$features<-fullColumnIndices[bicluster_results[[bcIndex]]$features]
  }
  
}

clusteringName <-
  paste0("/storage/biclustering_results/", "ALL_DNDS_ALL_NETWORKS_", paste(gsub("\\.", "_", c(
    variableClassConsistencyFactors
  )), collapse = "_"))
saveRDS(bicluster_results, paste0(clusteringName, ".rds"))

print(paste0("Generate heatmap for ", length(bicluster_results)," biclusters...."))

dataSet_with_biclusterRaw<-cbind(dnds_per_gene,gene_wise_network_correlation)
dataSet_with_biclusterOrg<-cbind(dnds_per_gene,gene_wise_network_correlationOrg)
dataSet_with_bicluster<-dataSet_with_biclusterRaw
colorSet<-colorpanel(10, "white","black")
for(biclusterIndex in 1:length(bicluster_results)){
  dataSet_with_bicluster[bicluster_results[[biclusterIndex]]$samples,bicluster_results[[biclusterIndex]]$features]<-1+biclusterIndex
  colorSet<-c(colorSet,rep(rainbow(length(bicluster_results))[biclusterIndex],10))
}

distfun<- function(x) dist(dataSet_with_bicluster[filteredGenes,])
hclustfunc<-function(x) hclust(dist(as.matrix(dataSet_with_bicluster[filteredGenes,])))
RowvHelp <- as.dendrogram(hclust(dist(as.matrix(dataSet_with_bicluster[filteredGenes,]))))


heatmap.2(as.matrix(dataSet_with_biclusterRaw[filteredGenes,]),scale="none",distfun=distfun,hclustfun=hclustfunc,dendrogram="row",Colv=NA,Rowv=RowvHelp, trace="none",col=colorpanel(30, "white","black"),margins = c(12, 12))

heatmap.2(as.matrix(dataSet_with_bicluster[filteredGenes,]),distfun=distfun,hclustfun=hclustfunc,Rowv=RowvHelp,scale="none",dendrogram="row",Colv=NA, trace="none",col=colorSet,margins = c(12, 12))



################################################################
####################### Real data Single Cell ##################
################################################################

###REMOVE THIS FOR PAPER###

prefix<-"NGT_SPLIT_TREE_SCR9000_NAto0_"
dndsColsToUse<-c(11, 12, 13, 14, 15, 18, 19, 16, 17 ,20, 21, 10)
singleCellColsToUse<-c()
includeNAGenes<-TRUE

# Biclustering with real data. To achieve good results, one has to play with the variableClassConsistencyFactors variable and
# MAYBE (but not likely) with the variableClassAllowZeroSum and amountOfBicluserSearches
# Load data

print("Load DNDS")
dnds_per_gene<-read.csv2(paste0("storage//csvs//",prefix,"normalized_dnds_per_gene.csv"),row.names=1)

if (length(dndsColsToUse) == 0) {
  dndsColsToUse <- 1:ncol(dnds_per_gene)
}


dndsColsToUseOrg <- dndsColsToUse

if(max(dndsColsToUse)>ncol(dnds_per_gene)){
  print(paste0(
    "Some dndsColsToUse where ommited since they did not point to a usable column: ",
    paste0(dndsColsToUse[dndsColsToUse>ncol(dnds_per_gene)], collapse = "_")
  ))
  dndsColsToUse<-dndsColsToUse[dndsColsToUse<=ncol(dnds_per_gene)]
}

dndsColsToUse <-
  dndsColsToUse[apply(dnds_per_gene[, dndsColsToUse], 2, function(x) {
    sum(!is.na(x))
  }) > 0]

if (length(dndsColsToUse) < length(dndsColsToUseOrg)) {
  print(paste0(
    "Removed ",
    length(setdiff(dndsColsToUseOrg, dndsColsToUse)),
    " DNDS columns for being only NA: ",
    paste0(setdiff(dndsColsToUseOrg, dndsColsToUse), collapse =
             " ")
  ))
}

print(paste0(
  "Used DNDS column indices: ",
  paste0(dndsColsToUse, collapse = "_")
))
print(paste0("Used DNDS column names: ", paste0(colnames(dnds_per_gene)[dndsColsToUse], collapse =
                                                  "_")))

print("Load gene wise correlation singleCells")
scDataPerRegionAndCelltype <-
  read.csv2(paste0("storage//csvs//",prefix, "scDataPerRegionAndCelltype.csv"),
            row.names = NULL)[,-1]

scDataPerRegionAndCelltypeOrg <- scDataPerRegionAndCelltype

#Rank normalization
scDataPerRegionAndCelltype[,]<-rank(scDataPerRegionAndCelltype,ties.method="max",na.last="keep")/sum(!is.na(scDataPerRegionAndCelltype))
scDataPerRegionAndCelltype[scDataPerRegionAndCelltypeOrg == 0] <- 0

#Remove all genes that have in general a small expression singleCellData
filteredGenes <-
  apply(scDataPerRegionAndCelltype, 1, function(x) {
    sum(x == 0, na.rm = TRUE) < sum(!is.na(x))
  })

print(paste0(
  "Filtered ",
  sum(filteredGenes == FALSE),
  " from ",
  nrow(scDataPerRegionAndCelltype),
  " for small expression"
))

if (includeNAGenes == FALSE) {
  #Remove all genes with NA values
  filteredGenesByNA <-
    apply(dnds_per_gene, 1, function(x) {
      sum(is.na(x))
    }) == 0
  print(paste0(
    "Filtered ",
    sum(filteredGenesByNA == FALSE),
    " from ",
    nrow(dnds_per_gene),
    " for NA"
  ))
  
  filteredGenes <- filteredGenes & filteredGenesByNA
}

print(paste0("Chromosome size: ", sum(filteredGenes)))

scDataPerRegionAndCelltype[!filteredGenes, ] <-
  NA #set rest to NA, but it's not necessary since it is anyway excluded


if (length(singleCellColsToUse) == 0) {
  singleCellColsToUse <- 1:ncol(scDataPerRegionAndCelltype)
}

if(max(singleCellColsToUse)>ncol(scDataPerRegionAndCelltype)){
  print(paste0(
    "Some singleCellColsToUse where ommited since they did not point to a usable column: ",
    paste0(singleCellColsToUse[singleCellColsToUse>ncol(scDataPerRegionAndCelltype)], collapse = "_")
  ))
  singleCellColsToUse<-singleCellColsToUse[singleCellColsToUse<=ncol(scDataPerRegionAndCelltype)]
}

singleCellColsToUseOrg <- singleCellColsToUse
singleCellColsToUse <-
  singleCellColsToUse[apply(scDataPerRegionAndCelltype[, singleCellColsToUse], 2, function(x) {
    sum(!is.na(x))
  }) > 0]

if (length(singleCellColsToUse) < length(singleCellColsToUseOrg)) {
  print(paste0(
    "Removed ",
    length(setdiff(
      singleCellColsToUseOrg, singleCellColsToUse
    )),
    " single cell columns for being only NA: ",
    paste0(
      setdiff(singleCellColsToUseOrg, singleCellColsToUse),
      collapse =
        " "
    )
  ))
}

print(paste0(
  "Used singleCell column indices: ",
  paste0(singleCellColsToUse, collapse = "_")
))
print(paste0(
  "Used singleCell column names: ",
  paste0(colnames(scDataPerRegionAndCelltype)[singleCellColsToUse], collapse =
           "_")
))


fullColumnIndices <-
  c(dndsColsToUse, (singleCellColsToUse) + ncol(dnds_per_gene))

print("Define parameters")

#define parameters
variableClassConsistencyFactors <-
  c(0.9,0.75) #one needs to try out different consistency factors - FIND "OPTIMAL" VALUES
variableClasses <-
  c(rep(1, ncol(dnds_per_gene[, dndsColsToUse])), rep(2, ncol(scDataPerRegionAndCelltype[, singleCellColsToUse]))) #the 2 variable classes - don't change this
variableClassConsistencyFactorsTypes <-
  c(2, 2) #one needs to try out different consistency factors - FIND "OPTIMAL" VALUES
variableClassAllowZeroSum <-
  c(FALSE, FALSE) #GO terms is very sparse, therefore no "zero GO columns" in bicluster - no need to change this
variableClassNotOnTabuList <-
  c(TRUE, FALSE) #GO terms is very sparse, therefore no "zero GO columns" in bicluster - no need to change this

amountOfBicluserSearches <-
  50  #can be higher, but you should look at the bicluster scores - FIND A "SUITABLE" AMOUNT
minimumAmountOfVariableClasses <-
  2 #we want biclusters with at least 3 types of variables.

#build dataSet
dataSet <-
  cbind(dnds_per_gene[, dndsColsToUse], scDataPerRegionAndCelltype[, singleCellColsToUse])
print(paste0("DNDS columns have ", paste0(sum(
  is.na(dnds_per_gene[, dndsColsToUse])
)), " NAs"))


#Execute biclustering

bicluster_results <-
  weighted_GABi(
    dataSet[filteredGenes, ],
    variableClasses,
    variableClassConsistencyFactors,
    variableWeightSimilar,
    variableClassAllowZeroSum,
    variableClassNotOnTabuList,
    variableClassConsistencyFactorsTypes,
    minimumAmountOfVariableClasses,
    amountOfBicluserSearches,
    useCores = 1,
    popsize = dim(dataSet[filteredGenes, ])[1] * 250,
    nsubpops = round(dim(dataSet[filteredGenes, ])[1] * 250 / 80)
  )

print(paste0("Bicluster found: ",length(bicluster_results)))
#adapt sample index for the filtered genes/samples/rows
if(length(bicluster_results)>0){
  for(bcIndex in 1:length(bicluster_results)){
    for(actFiltered in which(filteredGenes==FALSE)){
      bicluster_results[[bcIndex]]$samples[bicluster_results[[bcIndex]]$samples>=actFiltered]<-bicluster_results[[bcIndex]]$samples[bicluster_results[[bcIndex]]$samples>=actFiltered]+1
    }
  }
  print(paste0("Filtered ",sum(filteredGenes==FALSE)," genes of ",length(bicluster_results)," results!"))  
  for(bcIndex in 1:length(bicluster_results)){
    bicluster_results[[bcIndex]]$features<-fullColumnIndices[bicluster_results[[bcIndex]]$features]
  }
  
}


clusteringName <-
  paste0("/storage/biclustering_results/", "ALL_DNDS_ALL_CELLTYPES_", paste(gsub("\\.", "_", c(
    variableClassConsistencyFactors
  )), collapse = "_"))
saveRDS(bicluster_results, paste0(clusteringName, ".rds"))

print(paste0("Generate heatmap for ", length(bicluster_results)," biclusters...."))

dataSet_with_biclusterRaw<-cbind(dnds_per_gene,scDataPerRegionAndCelltype)
dataSet_with_biclusterOrg<-cbind(dnds_per_gene,scDataPerRegionAndCelltypeOrg)
dataSet_with_bicluster<-dataSet_with_biclusterRaw
colorSet<-colorpanel(10, "white","black")
for(biclusterIndex in 1:length(bicluster_results)){
  dataSet_with_bicluster[bicluster_results[[biclusterIndex]]$samples,bicluster_results[[biclusterIndex]]$features]<-1+biclusterIndex
  colorSet<-c(colorSet,rep(rainbow(length(bicluster_results))[biclusterIndex],10))
}

distfun<- function(x) dist(dataSet_with_bicluster[filteredGenes,])
hclustfunc<-function(x) hclust(dist(as.matrix(dataSet_with_bicluster[filteredGenes,])))
RowvHelp <- as.dendrogram(hclust(dist(as.matrix(dataSet_with_bicluster[filteredGenes,]))))


heatmap.2(as.matrix(dataSet_with_biclusterRaw[filteredGenes,]),scale="none",distfun=distfun,hclustfun=hclustfunc,dendrogram="row",Colv=NA,Rowv=RowvHelp, trace="none",col=colorpanel(30, "white","black"),margins = c(12, 12))

heatmap.2(as.matrix(dataSet_with_bicluster[filteredGenes,]),distfun=distfun,hclustfun=hclustfunc,Rowv=RowvHelp,scale="none",dendrogram="row",Colv=NA, trace="none",col=colorSet,margins = c(12, 12))

