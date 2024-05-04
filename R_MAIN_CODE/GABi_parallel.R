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
# Parallel version of "GABi_fixed", which are minor bugfixes for GABi (https://cran.r-project.org/web/packages/GABi). 
# Changes from "GABi" to "GABi_fixed" are documented in "GABi_fiexed". The optimization for parallel computing, which
# has been done in this file is not documented in detail, since it required to change a large portion of the original code.
GABi_parallel <- function (x,
                             convergenceGens = 40,
                             popsize = 256,
                             mfreq = 1,
                             xfreq = 0.5,
                             maxNgens = 200,
                             identityThreshold = 0.75,
                             nsubpops = 4,
                             experiod = 10,
                             diffThreshold = 0.9,
                             useCores = 1,
                             libraryLocations = NULL,
                             clusteringName = "",
                             minBiclusterSampleSize = 2, #minimum amount of samples that a bicluster should have
							 useTabuSamples = TRUE, #use tabu list for samples 
                             verbose = FALSE,
                             maxLoop = 1,
                             fitnessArgs = list(consistency = 0.8,
                                                featureWeights = rowMeans(x, na.rm = TRUE)),
                             fitnessFun = getFitnesses.entropy,
                             featureSelFun = featureSelection.basic)
  {
	#optimized version of the fps function by using raw values instead of bool (much better memory consumtion)
    fpsRaw <- function (population, fitnesses, elitism)
    {
      popsize <- nrow(population)
      goodSols <- which(fitnesses > 0)
      probShare <- fitnesses[goodSols] / sum(fitnesses[goodSols])
      cutoffs <- cumsum(probShare)
      intpop <- array(as.raw(0), dim = dim(population))
      if (elitism) {
        fittest <-
          sort(fitnesses,
               index.return = TRUE,
               decreasing = TRUE)$ix[1]
        intpop[1, ] <- population[fittest, ]
        selectionPoints <- runif(nrow(intpop) - 1)
        selectedSols <- goodSols[unlist(lapply(selectionPoints,
                                               function(x, cutoffs)
                                                 min(which(cutoffs > x)), cutoffs = cutoffs))]
        intpop[c(2:popsize), ] <- population[selectedSols, ]
      }
      else {
        selectionPoints <- runif(nrow(intpop))
        selectedSols <- goodSols[unlist(lapply(selectionPoints,
                                               function(x, cutoffs)
                                                 min(which(cutoffs > x)), cutoffs = cutoffs))]
        intpop <- population[selectedSols, ]
      }
      intpop
    }
    
	#optimized version of the reproduction function by using raw values instead of bool (much better memory consumtion)
    reproductionRaw <-
      function (population,
                xfreq,
                mfreq,
                xoverpoints,
                pinvert,
                elitism)
      {
        newpop <- array(as.raw(0), dim = dim(population))
        popsize <- dim(population)[1]
        if (elitism) {
          xindex <- c(2:popsize)[runif((popsize - 1)) > xfreq]
          xindex <- c(1, xindex)
        }
        else {
          xindex <- c(1:popsize)[runif(popsize) > xfreq]
        }
        newpop[xindex, ] <- population[xindex, ]
        if (length(xindex) < popsize) {
          xindex <- sample(setdiff(c(1:popsize), xindex))
          newpop[xindex, ] <-
            crossoverRaw(population[xindex, ], xoverpoints,
                         pinvert)
        }
        if (elitism) {
          newpop[c(2:popsize), ] <- mutationRaw(newpop[c(2:popsize), ], mfreq)
        }
        else {
          newpop <- mutationRaw(newpop, mfreq)
        }
        newpop
      }
    
	#optimized version of the crossover function by using raw values instead of bool (much better memory consumtion)
    crossoverRaw <- function (subpop, xoverpoints, pinvert = 0)
    {
      nsamples <- dim(subpop)[2]
      newpop <- array(as.raw(2), dim(subpop))
      if (identical(dim(subpop), NULL)) {
        newpop <- subpop
      }
      else {
        for (i in 1:floor(dim(subpop)[1] / 2)) {
          mother <- subpop[(2 * i), ]
          father <- subpop[(2 * i) - 1, ]
          xpoints <-
            round(runif(xoverpoints, min = 1, max = nsamples))
          xpoints <- c(xpoints, nsamples)
          new1 <- rep(as.raw(2), nsamples)
          new2 <- new1
          m <- 1
          parent <- "mum"
          
          for (j in xpoints) {
            if (parent == "mum") {
              if (runif(1) < pinvert) {
                new1[c(m:j)] <- mother[c(j:m)]
              }
              else {
                new1[c(m:j)] <- mother[c(m:j)]
              }
              if (runif(1) < pinvert) {
                new2[c(m:j)] <- father[c(m:j)]
              }
              else {
                new2[c(m:j)] <- father[c(m:j)]
              }
              parent <- "dad"
            }
            else {
              if (runif(1) < pinvert) {
                new1[c(m:j)] <- father[c(j:m)]
              }
              else {
                new1[c(m:j)] <- father[c(m:j)]
              }
              if (runif(1) < pinvert) {
                new2[c(m:j)] <- mother[c(m:j)]
              }
              else {
                new2[c(m:j)] <- mother[c(m:j)]
              }
              parent <- "mum"
            }
            m <- j + 1
          }
          
          newpop[(2 * i), ] <- new1[1:dim(subpop)[2]]
          newpop[(2 * i) - 1, ] <- new2[1:dim(subpop)[2]]
        }
        if ((newpop[dim(subpop)[1], 1]) == as.raw(2)) {
          newpop[dim(subpop)[1], ] <- subpop[dim(subpop)[1], ]
        }
      }
      newpop
    }
    
	#optimized version of the mutation function by using raw values instead of bool (much better memory consumtion)
    mutationRaw <- function (pop, mfreq)
    {
      changes <- array(as.vector(runif(pop) < mfreq, mode = "raw"),
                       dim = dim(pop))
      array(as.vector(xor(pop, changes), mode = "raw"), dim = dim(pop))
    }
    
	#optimized version of the exchangeSolsRaw function by using raw values instead of bool (much better memory consumtion)
    exchangeSolsRaw <-
      function (demes,
                fitnesses,
                fittestonly,
                proximity)
      {
        nsubpops <- length(demes)
        if (fittestonly) {
          chrs <- list(rep(as.raw(2), nsubpops))
          swapindices <- rep(NA, nsubpops)
          for (i in 1:nsubpops) {
            fittest <- sort(fitnesses[[i]],
                            index.return = TRUE,
                            decreasing = TRUE)$ix[1]
            chrs[[i]] <- demes[[i]][fittest, ]
            swapindices[i] <- fittest
          }
        }
        else {
          chrs <- list(rep(as.raw(2), nsubpops))
          swapindices <- rep(NA, nsubpops)
          for (i in 1:nsubpops) {
            swapchoice <- runif(1, min = 1, max = dim(demes[[i]])[1])
            chrs[[i]] <- demes[[i]][swapchoice, ]
            swapindices[i] <- swapchoice
          }
        }
        if (proximity) {
          swapindex <- swapindices[1]
          demes[[1]][swapindex, ] <- chrs[[nsubpops]]
          for (i in 2:nsubpops) {
            swapindex <- swapindices[i]
            demes[[i]][swapindex, ] <- chrs[[i - 1]]
          }
        }
        else {
          listrefs <- sample(c(1:nsubpops))
          for (i in 1:nsubpops) {
            swapindex <- swapindices[i]
            demes[[i]][swapindex, ] <- chrs[[listrefs[i]]]
          }
        }
        demes
      }
    fitnessArgs$tabu <- list(features = 0, samples = 0)
    diffRows <- 1:nrow(x)
    # row.mins <- apply(x, MARGIN = 1, min, na.rm = TRUE)
    # row.maxes <- apply(x, MARGIN = 1, max, na.rm = TRUE)
    # invariant.rows <-
    #   which(mapply(function(x, y)
    #     x == y, x = row.mins,
    #     y = row.maxes))
    # diffRows <- setdiff(1:nrow(x), invariant.rows)
    # low.rows <-
    #   diffRows[which(((rowMeans(x[diffRows, ], na.rm = TRUE) -
    #                      row.mins[diffRows]) / (row.maxes[diffRows] - row.mins[diffRows])) <
    #                    (1 - diffThreshold))]
    # hi.rows <-
    #   diffRows[which(((rowMeans(x[diffRows, ], na.rm = TRUE) -
    #                      row.mins[diffRows]) / (row.maxes[diffRows] - row.mins[diffRows])) >
    #                    diffThreshold)]
    # cat(paste(
    #   length(low.rows),
    #   "constantly low features and",
    #   length(hi.rows),
    #   "constantly high\n"
    # ))
    # diffRows <-
    #   setdiff(1:nrow(x), c(low.rows, hi.rows, invariant.rows))
    # cat(paste(
    #   length(diffRows),
    #   "sufficiently variable features of",
    #   nrow(x),
    #   "\n"
    # ))
    if (length(diffRows) < 2) {
      warning("no sufficiently variable features in dataset")
      return(list())
    }
    subpopsize <- ceiling(popsize / nsubpops)
    demes <- list(rep(NA, nsubpops))
    featureSelFun <- featureSelFun
    fitnessFun <- fitnessFun
    environment(featureSelFun) <- environment()
    environment(fitnessFun) <- environment()
    allSols <- list()
    loopCount <- 0
    
    #For loading/storing intermediate results
    if (file.exists(paste0(clusteringName, "_tmp.rds")) &&
        nchar(clusteringName) > 1) {
      allSols <- readRDS(paste0(clusteringName, "_tmp.rds"))
      print(paste0("Already biclusters found: ", length(allSols)))
      
      newTabuFeatures <-
        unlist(lapply(allSols, function(x)
          x$features))
      #newTabuSamples <-
      #  unlist(lapply(allSols, function(x)
       #   x$samples))
      fullTabuFeatures <- union(fitnessArgs$tabu$features,
                                newTabuFeatures)
      #fullTabuSamples <-
      #  union(fitnessArgs$tabu$samples, newTabuSamples)
      fitnessArgs$tabu$features <- fullTabuFeatures
      #fitnessArgs$tabu$samples <- fullTabuSamples
      
      loopCount <-
        length(unique(sapply(allSols , function(x) {
          x$score
        })))
    }
    
    `%op%` <- if (useCores > 1)
      `%dopar%`
    else
      `%do%`
    
    if (useCores > 1) {
      clComputing <- 0
      clEvalRes <- NULL
      while (is.null(clEvalRes)) {
        tryCatch({
          print(paste("Allocating ",useCores," for parallel processing."))
          clComputing <-
            makeCluster(
              useCores,
              type = "FORK",
              outfile = paste0(clusteringName, "_parallel_log.txt"),
              timeout = 36000
            )
          Sys.sleep(5)
          registerDoParallel(clComputing)
          m <- matrix(rnorm(useCores*useCores), useCores, useCores)
          clEvalRes <- foreach(i=1:nrow(m), .combine=rbind) %dopar%{
            (m[i,] / mean(m[i,]))
          } 
        }, error = function(e) {
          clEvalRes<-NULL
          print(paste(
            "Allocate Clustering cores error:",
            e,
            "...retrying"
          ))
        })
      }
    }
    
    
    while (loopCount < maxLoop) {
      start_time_loop <- Sys.time()
      cat(paste("initialising population for GA loop", loopCount +
                  1, "\n"))
      for (j in 1:nsubpops) {
        demes[[j]] <- array(as.raw(0), dim = c(subpopsize, ncol(x)))
        for (i in 1:subpopsize) {
          demes[[j]][i, sample(ncol(x), 2)] <- as.raw(1)
        }
      }
      bestSols <- list()
      genNo <- 0
      keptSols <- list()
      for (j in 1:nsubpops) {
        bestSols[[j]] <- rep(0, ncol(demes[[j]]))
      }
      convergenceCounter <- 0
      notConverged <- TRUE
      transfergen <- experiod
      fitnesshash <- hash()
      cat("evolving GA solutions\n")
      while (notConverged & (genNo < maxNgens)) {
        cat(paste("trying generation no.", genNo, "\n"))
        genNo <- genNo + experiod
        all_fitnesses <- list()
        stagnancy <- rep(FALSE, nsubpops)
        
        start_time <- Sys.time()
        #parallel computing for demes
        demes_and_fitnesses <-
          foreach(
            actDeme = demes,
            j = icount(),
            .combine = "c",
            .inorder = TRUE,
            .multicombine = TRUE
          ) %op% {
            for (k in 1:experiod) {
              fitnesses <- rep(0, nrow(actDeme))
              fitnesses <- apply(actDeme, MARGIN = 1, fitnessFun)
              fitnesses[is.infinite(fitnesses)] <- 0
              fitnesses[is.na(fitnesses)] <- 0
              new_bestSol <- actDeme[which(fitnesses ==
                                             max(fitnesses))[1], ]
              identical <- as.numeric(unlist(lapply(bestSols,
                                                    function(x) {
                                                      identical(x, new_bestSol)
                                                    })))
              if (sum(identical) > 0) {
                stagnancy[j] <- TRUE
              }
              else {
                convergenceCounter <- 0
              }
              
              if (convergenceCounter < convergenceGens) {
                if (sum(fitnesses) == 0) {
                  if (verbose)
                    cat("re-initialising population\n")
                  actDeme <- array(as.raw(0), dim = c(subpopsize,
                                                      ncol(x)))
                  for (i in 1:subpopsize) {
                    actDeme[i, sample(ncol(x), 2)] == as.raw(1)
                  }
                  
                  fitnesses <- rep(0, nrow(actDeme))
                  
                  fitnesses <-
                    apply(actDeme, MARGIN = 1, fitnessFun)
                  fitnesses[is.infinite(fitnesses)] <- 0
                  fitnesses[is.na(fitnesses)] <- 0
                  
                  
                }
                
                f <- fitnesses / mean(fitnesses[fitnesses > 0])
                
                intpop <- fpsRaw(
                  population = actDeme,
                  fitnesses = f,
                  elitism = TRUE
                )
                actDeme <-
                  reproductionRaw(
                    intpop,
                    xfreq,
                    mfreq / ncol(intpop),
                    xoverpoints = 1,
                    pinvert = 0,
                    elitism = TRUE
                  )
              }
              else {
                if (verbose)
                  cat(paste(
                    "population convergence after",
                    genNo,
                    "generations\n"
                  ))
                notConverged <- FALSE
              }
            }
            
            return(list(
              list(
                deme = actDeme,
                fitnesses = fitnesses,
                new_bestSol = new_bestSol,
                stagnancy = stagnancy[j],
                notConverged = notConverged
              )
            ))
          }
        
        
        if (length(demes_and_fitnesses) < 1) {
          if (useCores > 1) {
            stopCluster(clComputing)
          }
          return(allSols)
        }
        #Setting returns from parallel computing
        demes <- lapply(demes_and_fitnesses, function(x) {
          return(x$deme)
        })
        all_fitnesses <- lapply(demes_and_fitnesses, function(x) {
          return(x$fitnesses)
        })
        
        if (verbose) {
          cat("fitnesses:\n")
          cat(round(sort(
            sapply(demes_and_fitnesses, function(x) {
              return(max(x$fitnesses))
            }),
            decreasing = TRUE
          )))
          cat("\n")
        } else{
          cat(paste())
        }
        
        bestSols <- lapply(demes_and_fitnesses, function(x) {
          return(x$new_bestSol)
        })
        
        stagnancy <- sapply(demes_and_fitnesses, function(x) {
          return(x$stagnancy)
        })
        notConverged <-
          sum(sapply(demes_and_fitnesses, function(x) {
            return(x$notConverged)
          })) == nsubpops
        
        if (sum(as.numeric(stagnancy)) == nsubpops) {
          convergenceCounter <- convergenceCounter + experiod
        } else{
          convergenceCounter <- 0
        }
        if (genNo == transfergen) {
          demes <- exchangeSolsRaw(demes, all_fitnesses, TRUE,
                                   FALSE)
          transfergen <- transfergen + experiod
        }
        cat(paste(
          experiod,
          "gens took",
          difftime(Sys.time(), start_time, units = "secs"),
          "secs, best fitness:",
          (max(
            sapply(demes_and_fitnesses, function(x) {
              return(max(x$fitnesses))
            })
          )),
          "\n"
        ))
      }
      if (verbose)
        cat("out of convergence loop\n")
      
      start_time <- Sys.time()
      all_fitnesses <-
        foreach(
          actDeme = demes,
          j = icount(),
          .combine = "c",
          .inorder = TRUE,
          .multicombine = TRUE
        ) %op% {
          return(list(apply(actDeme, MARGIN = 1, fitnessFun)))
        }
      
      cat(paste(
        "final fitness calc took ",
        difftime(Sys.time(), start_time, units = "secs"),
        "secs\n"
      ))
      
      maxFitness <- max(sapply(all_fitnesses, function(x) {
        return(x)
      }))
      cat(paste("Best fitness: ", maxFitness, "\n"))
      
      #Computed sols with maxFitness that are no duplicates by inflating the all_fitnesses and demes from subpops
      start_time <- Sys.time()
      goodsolsCombined <-
        which(unlist(all_fitnesses) == maxFitness &
                duplicated(do.call(rbind, demes)) == FALSE)
      cat(
        paste(
          "Check for (",
          sum(unlist(all_fitnesses) == maxFitness),
          "-->",
          length(goodsolsCombined),
          ") unique goodsols ",
          difftime(Sys.time(), start_time, units = "mins"),
          "mins\n"
        )
      )
      ######
      
      start_time <- Sys.time()
      if (maxFitness > 0) {
        keptSols <- list()
        for (j in 1:nsubpops) {
          #goodsols<-which(all_fitnesses[[j]]==maxFitness)
          #goodsols<-goodsols[duplicated(demes[[j]][goodsols,])==FALSE] #removes duplicates
          
          #unique max fitnesses have been already computed, now they have to be extracted out of goodsolsCombined
          solsBeforeDemej <- 0
          if (j > 1) {
            solsBeforeDemej <- sum(sapply(demes[1:(j - 1)], nrow))
          }
          endOfSolsDemej <- sum(sapply(demes[1:j], nrow))
          goodsols <-
            goodsolsCombined[goodsolsCombined > solsBeforeDemej &
                               goodsolsCombined <= endOfSolsDemej] - solsBeforeDemej
          
          if (identityThreshold == 1) {
            newSols <- lapply(goodsols, function(x) {
              list(fitness = all_fitnesses[[j]][x],
                   chr = demes[[j]][x, ])
            })
            keptSols <- append(keptSols, newSols)
            goodsols <- c()
          }
          
          ######
          
          if (length(goodsols) > 0) {
            if (length(keptSols) == 0) {
              keptSols[[1]] <-
                list(chr = demes[[j]][goodsols[1],], fitness = all_fitnesses[[j]][goodsols[1]])
              goodsols <- goodsols[-1]
            }
            if (length(goodsols) > 0) {
              for (i in 1:length(goodsols)) {
                identity <- lapply(keptSols, function(x, comp) {
                  # This similarity will also be 1 if the smaller is contained in the larger
                  # Since the score is based on the size, this shouldn't be an issue
                  # Could influence stuff on tabu list!
                  sum(colMeans(rbind(
                    x$chr == as.raw(1), comp == as.raw(1)
                  )) > 0.5) / min(c(sum(x$chr == as.raw(
                    1
                  ))),
                  sum(comp ==
                        as.raw(1)))
                }, comp = demes[[j]][goodsols[i],])
                toRemove <- which(identity > identityThreshold)
                
                newSols <- list()
                if (length(toRemove) > 0) {
                  newSols <- keptSols[toRemove]
                  keptSols <- keptSols[-toRemove]
                }
                newSols[[length(newSols) + 1]] <-
                  list(fitness = all_fitnesses[[j]][goodsols[i]],
                       chr = demes[[j]][goodsols[i], ])
                removeScores <- unlist(lapply(newSols, function(x) {
                  x$fitness
                }))
                ReplaceSol <-
                  which(removeScores == max(removeScores))[1]
                keptSols[[length(keptSols) + 1]] <-
                  newSols[[ReplaceSol]]
                
              }
            }
          }
        }
      }
      cat(paste(
        "Check similar identity calc took ",
        difftime(Sys.time(), start_time, units = "mins"),
        "mins\n"
      ))
      
      if (length(keptSols) > 0 &&
          sum(sapply(keptSols, function(x) {
            sum(x$chr == as.raw(1))
          }) > minBiclusterSampleSize) > 0) {
        keptSols<-keptSols[sapply(keptSols, function(x) {
          sum(x$chr == as.raw(1))
        }) > minBiclusterSampleSize]
        keptSols <- lapply(keptSols, function(x) {
          list(
            samples = which(x$chr == as.raw(1)),
            features = featureSelFun(which(x$chr ==
                                             as.raw(1))),
            score = fitnessFun(x$chr)
          )
        })
      } else{
        cat(paste(
          "No new biclusters can be found after ",
          loopCount,
          " loops!"
        ))
        if (useCores > 1) {
          stopCluster(clComputing)
        }
        return(allSols)
      }
      
      
      allSols <- c(allSols, keptSols)
      cat(paste(
        "Loop calculation took ",
        difftime(Sys.time(), start_time_loop, units = "secs"),
        "secs \n"
      ))
      cat(paste(
        length(keptSols),
        "biclusters discovered in loop, total",
        length(allSols),
        "\n"
      ))
	  
      newTabuFeatures <-
        unlist(lapply(keptSols, function(x)
          x$features))
      newTabuSamples <-
        unlist(lapply(keptSols, function(x)
         x$samples))
      fullTabuFeatures <- union(fitnessArgs$tabu$features,
                                newTabuFeatures)
      fullTabuSamples <-
        union(fitnessArgs$tabu$samples, newTabuSamples)
		if(is.null(rownames(x))){
		  rownames(x)<-1:nrow(x)
		}
			  
	  for(actKeptSol in keptSols){
		  cat(paste(
			  "New bicluster (features&sampleAmount):",
			  paste(c(rownames(x)[sort(actKeptSol$features)],length(actKeptSol$samples)),collapse = ", "),
			  "\n"
			))
	  }
	  
		cat(paste(
		  "Tabu features:",
		  paste(rownames(x)[fullTabuFeatures],collapse = ", "),
		  "\n"
		))
      fitnessArgs$tabu$features <- fullTabuFeatures
	  if(useTabuSamples){
        fitnessArgs$tabu$samples <- fullTabuSamples
	  }
      cat("updating arguments to re-run biclustering loop\n")
      environment(featureSelFun) <- environment()
      environment(fitnessFun) <- environment()
      loopCount <- loopCount + 1
      
      #For loading/storing intermediate results
      if (nchar(clusteringName) > 1) {
        saveRDS(allSols, paste0(clusteringName, "_tmp.rds"))
      }
      
    }
    if (useCores > 1) {
      stopCluster(clComputing)
    }
    
    allSols
  }