# Copyright (C) 2016 VRVis.
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


setsToRunFolder<-"..\\storage\\setsToRunHuman"

while(length(list.files(setsToRunFolder))>0){
  
  runGeneticAlgorithm<-FALSE
  
  actGeneList<-list.files(setsToRunFolder)[length(list.files(setsToRunFolder))]
  if(substr(actGeneList,1,4)=="_GA_"){
    actGeneList<-list.files(setsToRunFolder)[1]
  }
  
  print(paste0("Run setlist: ",actGeneList))
  if(substr(actGeneList,1,4)=="_GA_"){
    runGeneticAlgorithm<-TRUE
    print("Geneset selection activated!")
  }
  
  useCores<-4                                #cores that will be used for p-value computation and clustering. 
  #If one chooses 1, no parallel computing will be used
  
  precision<-5000                                  #Amount of random gene-sets (synergies) that will be drawn for p-value computation
  
  splitInBatches<-1                                #Random gene-sets (synergies) will be stored in one matrix with "precision" cols.
  #This can be a problem regarding the size of the memory, therefore
  #you can split the computation into batches. Batchsize = precision/splitInBatches
  #The batchsize shouldn't be bigger than 1000
  
  
  libFolder<-NULL
  
  testSetName<-paste0(setsToRunFolder,"//",actGeneList)                               #Name of the testset (without file ending) 

  #(what will be stored by get_gene_expression_for_genesets)
  
  resultsFolder<-paste0("../storage//results_functionalmaps_human//",tools::file_path_sans_ext(actGeneList))
  
  allGenesData<-"../storage//all_genes_expression_Human.mat"                                 #Folder with random genes
  #(what will be stored by get_gene_expression_for_random_genes)
  
  trimExpression<-0.05                                                #Trim-factor for trimmed-mean gene-expression calculation of a gene set
  #0.05 --> The upper 5% and lower 5% percent gene-expression for every grid point
  #will be filtered. 
  
  minimumPValueToVisualize<-10^(-8)
  
  quantileScale<-TRUE
  getColorScale<-function(amountOfColors,threshold){
    #return(jet.colors(amountOfColors))
    greyColor<-colorpanel(amountOfColors,"#f9f9f9","#a8a8a8","#7a7a7a")
    heatColor<-colorpanel(amountOfColors,"#ffffb2","#feb24c","#bd0026")
    
    if(is.na(threshold))
      return(heatColor)
    if(threshold<=1 || threshold==Inf)
      return(heatColor)
    return(c(greyColor[1:floor(threshold)],heatColor[(floor(threshold)+1):amountOfColors]))
  }
  
  
  if(sum(!(c("XML","R.matlab","rjson","circlize","foreach","gplots","png","iterators","fields","igraph","fastcluster","clue","R.utils","GA") %in% installed.packages(lib.loc=libFolder)[,"Package"]))>0){
    install.packages("iterators", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("fields", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("png_0.1-7.tar.gz", repos = NULL, lib=libFolder)
    install.packages("foreach", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("R.matlab", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("rjson", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("doParallel", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("GA", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("gplots", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("fastcluster", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("clue", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("R.utils", repos = "http://cran.at.r-project.org/", lib=libFolder)
    install.packages("XML", repos = "http://cran.at.r-project.org/", lib=libFolder)

  }
  
  
  library(R.matlab, lib.loc=libFolder)
  library(rjson, lib.loc=libFolder)
  library(foreach, lib.loc=libFolder)
  library(doParallel, lib.loc=libFolder)
  library(gplots, lib.loc=libFolder)
  library(png, lib.loc=libFolder)
  library(iterators, lib.loc=libFolder)
  library(fields, lib.loc=libFolder)
  library(fastcluster, lib.loc=libFolder)
  library(clue, lib.loc=libFolder)
  library(GA, lib.loc=libFolder)
  library(R.utils, lib.loc=libFolder)
  library(rsvg, lib.loc=libFolder) 
  library(XML, lib.loc=libFolder)
  
  
  setFile<-read.csv2(paste0(setsToRunFolder,"//",actGeneList),header=FALSE,stringsAsFactors=FALSE)
  setFile[,1][is.na(setFile[,1])]<-""
  setnames<-setFile[is.na(setFile[,2])&nchar(setFile[,1])>0,1]
  startsOfSets<-(1:nrow(setFile))[is.na(setFile[,2])&nchar(setFile[,1])>0]+1
  endsOfSets<-c(startsOfSets[2:length(startsOfSets)]-3,max((1:nrow(setFile))[!is.na(setFile[,2])]))
  
  
  `%op%` <- if (useCores>1) `%dopar%` else `%do%`
  
  #ontology file (gives information for every region in the Allen brain atlas)
  document <- fromJSON(file="../storage//ontologyHuman.json", method='C')
  
  #gets acronym of a region given by a region id (from ontology file)
  getAcronymByID <- function(child,ID){
    if(!is.null(child)){
      if(child$id==ID){
        return(child)
      }
    }
    
    for(actchild in child$children){
      res<-getAcronymByID(actchild,ID)
      if(!is.null(res)){
        if(res$id==ID){
          return(res)
        }
      }
    }
    return(NULL)
  }
  
  #create results folder
  dir.create(resultsFolder)
  
  #gives a vector with the size of the atlas (3D volume rescaled to 1D vector) where every position
  #is zero, execpt for positions where the region (given by the ID) can be found
  getAtlasRegionsOfID <- function(atlasRegions,ID){
    newAtlasRegions<-atlasRegions==ID
    
    childrens<-getAcronymByID(document$msg[[1]],ID)$children
    while(length(childrens)>0){
      newChildren<-c()
      for(i in 1:length(childrens)){
        newAtlasRegions<-newAtlasRegions|(atlasRegions==childrens[[i]]$id)
        newChildren<-c(newChildren,getAcronymByID(document$msg[[1]],childrens[[i]]$id)$children)
      }
      
      childrens<-newChildren
    }
    return(newAtlasRegions)
  }
  
  
  ########################################################################################################################
  # Functions for plotting
  ########################################################################################################################
  
  recursiveColorChange<-function(actNode,structureIDs,colorsOfStructures){
    for(actChild in 1:xmlSize(actNode)){
      attribs<-xmlAttrs(actNode[[actChild]])
      
      if(!is.na(attribs["style"])){
        isStructure<-attribs["structure_id"]==structureIDs
        if(sum(isStructure)>0){
          #print(attribs["style"])
          if(!is.na(colorsOfStructures[isStructure]) && !(substr(attribs["style"],nchar(attribs["style"])-6,nchar(attribs["style"]))=="#f2f1f0")){
              if(!(substr(attribs["style"],nchar(attribs["style"])-6,nchar(attribs["style"]))=="#231f20")){
                if(!(substr(attribs["style"],nchar(attribs["style"])-6,nchar(attribs["style"]))=="#6d6e70")){
                  attribs["style"]<-paste0(substr(attribs["style"],0,nchar(attribs["style"])-7),colorsOfStructures[isStructure])
                }else{
                  attribs["style"]<-paste0(substr(attribs["style"],0,nchar(attribs["style"])-7),"#6d6e70")
                }
              }else{
                attribs["style"]<-paste0(substr(attribs["style"],0,nchar(attribs["style"])-7),"#231f20")
              }
              
           }else{
              attribs["style"]<-paste0(substr(attribs["style"],0,nchar(attribs["style"])-7),"#ffffff")
           }
          xmlAttrs(actNode[[actChild]])<-attribs
        }
      }
      if(xmlSize(actNode[[actChild]])>0){
        recursiveColorChange(actNode[[actChild]],structureIDs,colorsOfStructures)
      }
      
    }
  }
  
  recursiveGetAllStructures<-function(actNode){
    structureIDs<-c()
    for(actChild in 1:xmlSize(actNode)){
      attribs<-xmlAttrs(actNode[[actChild]])
      
      if(!is.na(attribs["style"])){
        structureIDs<-c(attribs["structure_id"],structureIDs)
      }
      if(xmlSize(actNode[[actChild]])>0){
        structureIDs<-unique(c(structureIDs,recursiveGetAllStructures(actNode[[actChild]])))
      }
      
    }
    return(structureIDs)
  }
  
  getSmallestParentAtlasRegions <- function(atlasRegionsHuman,ID){
    if(is.null(ID)){
      return(rep(FALSE,length(atlasRegionsHuman)))
    }
    newAtlasRegions<-getAtlasRegionsOfID(atlasRegionsHuman,ID)
    
    if(sum(newAtlasRegions>0)==0){
      newAtlasRegions<-getSmallestParentAtlasRegions(atlasRegionsHuman,getAcronymByID(document$msg[[1]],ID)$parent_structure_id)
    }
    
    return(newAtlasRegions)
  }
  
  #plots (stores a png to path_to_file) of spatial p-values (plotPvals) into several slices using the custom palette
  #A color bar will be added, that shows significance at fdr005 and fdr01.
  plot_mri_style_image<-function(path_to_file,plotPvals,fdr005,fdr01){
    
    #volume that will be plotted 
    m<-array(0, length(uniqueRegions))
    plotPvals[is.na(plotPvals)]<-1
  
    for(actRegion in 1:length(uniqueRegions)){
      m[actRegion]<-quantile(-log10(plotPvals[uniqueRegionAtlas[[actRegion]]]),0.9,na.rm=TRUE)
    }
    
    m[is.nan(m)]<-0
    m[is.na(m)]<-0
    m[m<0]<-0
  
    minimumPval<--log10(minimumPValueToVisualize)
    
    if(quantileScale){
      if(quantile(m[m>0],probs=0.99)<Inf){
        minimumPval<-quantile(m[m>0],probs=0.99)
      }
    }
    
    m[m>minimumPval]<-minimumPval
    
    png(paste0(path_to_file),width=2000,height=2666)  
    par(mfrow=c(4,3),mar = rep(0.1, 4)) 
    
    
    m<-round(m*1000)
    rb<-getColorScale(minimumPval*1000+1,-log10(fdr01)*1000+1)
    
    for (x in 1:11){

      fileName <- paste0("../storage/slice_svg/",x*9,".svg")
      
      xmlfile=xmlParse(fileName)
      xmltop = xmlRoot(xmlfile)
      
      recursiveColorChange(xmltop[['g']],uniqueRegions,rb[m+1])
      
      brainSVG <- rsvg(charToRaw(as(xmlfile, "character")),width=1000)
      plot(c(0,0), type="n", axes=F, xlab="", ylab="")
      rasterImage(brainSVG,1,-1,2,1)
    }
    image.scale <- function(z, zlim, col = heat.colors(12),
                            breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
      if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
      }
      if(missing(breaks) & !missing(zlim)){
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
      }
      if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
      }
      poly <- vector(mode="list", length(col))
      for(i in seq(poly)){
        poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
      }
      xaxt <- ifelse(horiz, "s", "n")
      yaxt <- ifelse(horiz, "n", "s")
      if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
      if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
      if(missing(xlim)) xlim=XLIM
      if(missing(ylim)) ylim=YLIM
      plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i",cex.axis=8, ...)  
      for(i in seq(poly)){
        if(horiz){
          polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
        }
        if(!horiz){
          polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
        }
      }
      if(!(fdr005==0)){
        lines(c(-log10(fdr005),-log10(fdr005)),c(0,1),lty=3,lwd=5)
      }
      if(!(fdr01==0)){
        lines(c(-log10(fdr01),-log10(fdr01)),c(0,1),lwd=5)
      }
    }
    
    par(mgp=c(3,6,0))
    par(mar=c(15,5,15,5))
    
    image.scale(c(0,minimumPval),col=rb)
    dev.off()
    
  }
  
	# This script is for maximizing contrast of gene expression synergy for gene sets.m)
	# The code for the rbga.int function is from the R package genalg (general code) and gramEvol (genalg extended for integer chromosome) adapted for parallel computing
	# genalg: https://github.com/egonw/genalg by Egon Willighagen and Michel Ballings
	# gramEvol: https://github.com/fnoorian/gramEvol/ by Farzad Noorian and Anthony Mihirana
  rbga.int <- function (ints=c(0,1), size = 10, suggestions = NULL, popSize = 200, iters = 100,
                      mutationChance = NA, elitism = NA, monitorFunc = NULL, ints.probs=rep(1/length(ints),length(ints)),
                      evalFunc = NULL, showSettings = FALSE, verbose = FALSE,useCores=1,dataMatrix=NULL,equalsIterationsForConvergence=10) {
  `%op%` <- if (useCores>1) `%dopar%` else `%do%`
  newIters<-0
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.")
  }
  vars = size
  if (is.na(mutationChance)) {
    mutationChance = 1/(vars + 1)
  }
  if (is.na(elitism)) {
    elitism = floor(popSize/5)
  }
  if (verbose)
    cat("Testing the sanity of parameters...\n")
  if (popSize < 5) {
    stop("The population size must be at least 5.")
  }
  if (iters < 1) {
    stop("The number of iterations must be at least 1.")
  }
  if (!(elitism < popSize)) {
    stop("The population size must be greater than the elitism.")
  }
  if (length(ints)<2) {
    stop("The number of integers must be at least 2.")
  }
  if (showSettings) {
    if (verbose)
      cat("The start conditions:\n")
    result = list(size = size, suggestions = suggestions,
                  popSize = popSize, iters = iters, elitism = elitism,
                  mutationChance = mutationChance)
    class(result) = "rbga"
    cat(summary(result))
  } else {
    if (verbose)
      cat("Not showing GA settings...\n")
  }
  if (vars > 0) {
    if (!is.null(suggestions)) {
      if (verbose)
        cat("Adding suggestions to first population...\n")
      population = matrix(nrow = popSize, ncol = vars)
      suggestionCount = dim(suggestions)[1]
      for (i in 1:suggestionCount) {
        population[i, ] = suggestions[i, ]
      }
      if (verbose)
        cat("Filling others with random values in the given domains...\n")
      for (child in (suggestionCount + 1):popSize) {
        population[child, ] = sample(ints, vars, rep = TRUE, prob=ints.probs)
      }
    } else {
      if (verbose)
        cat("Starting with random values in the given domains...\n")
      population = matrix(nrow = popSize, ncol = vars)
      for (child in 1:popSize) {
        population[child, ] = sample(ints, vars, rep = TRUE, prob=ints.probs)
      }
    }
    bestEvals = rep(NA, iters)
    meanEvals = rep(NA, iters)
    evalVals = rep(NA, popSize)
    for (iter in 1:iters) {
      if (verbose)
        cat(paste("Starting iteration", iter, "\n"))
      if (verbose)
        cat("Calucating evaluation values... ")
      
      clComputing <- c()
      if(useCores>1){
        clComputing <- makeCluster(useCores,type="PSOCK")
        registerDoParallel(clComputing)
      }
      evalVals<-foreach(allPop=1:popSize, i=icount(),.combine="c", .inorder=TRUE,.multicombine=TRUE,.export=c("geneExpressionSynergyOrg","expressionIndex","expressionOfGenes","trimExpression")) %op%{
        return(evalFunc(population[i,]))
      }
      if(useCores>1){
        stopCluster(clComputing)
      }
      
      bestEvals[iter] = min(evalVals)
      meanEvals[iter] = mean(evalVals)
      if (verbose)
        cat(" done.\n")
      if (!is.null(monitorFunc)) {
        if (verbose)
          cat("Sending current state to rgba.monitor()...\n")
        result = list(type = "integer chromosome", size = size,
                      popSize = popSize, iter = iter, iters = iters,
                      population = population, elitism = elitism,
                      mutationChance = mutationChance, evaluations = evalVals,
                      best = bestEvals, mean = meanEvals)
        class(result) = "rbga"
        monitorFunc(result)
      }
      
      newIters<-newIters+1
      
      if(iter>equalsIterationsForConvergence){
        if(length(unique(bestEvals[(iter-equalsIterationsForConvergence+1):iter]))==1){
          print("Converged!")
          break;
        }
      }
      
      if (iter < iters) {
        if (verbose)
          cat("Creating next generation...\n")
        newPopulation = matrix(nrow = popSize, ncol = vars)
        newEvalVals = rep(NA, popSize)
        if (verbose)
          cat("  sorting results...\n")
        sortedEvaluations = sort(evalVals, index = TRUE)
        sortedPopulation = matrix(population[sortedEvaluations$ix,], ncol = vars)
        if (elitism > 0) {
          if (verbose)
            cat("  applying elitism...\n")
          newPopulation[1:elitism, ] = sortedPopulation[1:elitism,]
          newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
        }
        if (vars > 1) {
          if (verbose)
            cat("  applying crossover...\n")
          for (child in (elitism + 1):popSize) {
            parentProb = dnorm(1:popSize, mean = 0, sd = (popSize/3))
            parentIDs = sample(1:popSize, 2, prob = parentProb)
            parents = sortedPopulation[parentIDs, ]
            # Crossover probability???
            crossOverPoint = sample(0:vars, 1)
            if (crossOverPoint == 0) {
              newPopulation[child, ] = parents[2, ]
              newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            } else if (crossOverPoint == vars) {
              newPopulation[child, ] = parents[1, ]
              newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            } else {
              newPopulation[child, ] = c(parents[1, ][1:crossOverPoint],
                                         parents[2, ][(crossOverPoint + 1):vars])
            }
          }
        }
        else {
          if (verbose)
            cat("  cannot crossover (#vars=1), using new randoms...\n")
          newPopulation[(elitism + 1):popSize, ] = sortedPopulation[sample(1:popSize,
                                                                           popSize - elitism), ]
        }
        population = newPopulation
        evalVals = newEvalVals
        
        if (mutationChance > 0) {
          if (verbose)
            cat("  applying mutations... ")
          mutatedCells=which(runif((popSize-elitism)*vars)<mutationChance)
          m_rows=elitism+ceiling(mutatedCells/vars)
          m_cols=mutatedCells%%vars
          m_cols=ifelse(m_cols==0,vars,m_cols)
          mutationCount=length(mutatedCells)
          for (muts in 1:mutationCount) {
            if (length(ints)==2) {
              population[m_rows[muts], m_cols[muts]]=ints[ints!=population[m_rows[muts], m_cols[muts]]]
            } else {
              population[m_rows[muts], m_cols[muts]] = sample(ints[ints!=population[m_rows[muts], m_cols[muts]]], 1, prob=ints.probs[ints!=population[m_rows[muts], m_cols[muts]]])
            }
          }

          if (verbose)
            cat(paste(mutationCount, "mutations applied\n"))
        }
      }
      
    }
  }
  result = list(type = "integer chromosome", size = size, popSize = popSize,
                iters = newIters, suggestions = suggestions, population = population,
                elitism = elitism, mutationChance = mutationChance, evaluations = evalVals,
                best = bestEvals, mean = meanEvals)
  class(result) = "rbga"
  return(result)
}
  
  
  
  #############################START OF THE METHOD#########################################
  
  
  print("Load atlas")
  atlasRegions<-readMat("../storage//atlasRegionsHuman.mat")
  
  #atlasRegions is a 3D volume reformated to a 1D vector
  #this is also done for gene expression volumes. In this script
  #we do not work with 3D volumes, since it's easier in R to work 
  #with vectors. The variable index maps [x,y,z] coordinates to the
  #1D vector 
  indexX<-atlasRegions$indexX
  indexY<-atlasRegions$indexY
  indexZ<-atlasRegions$indexZ
  atlasRegions<-as.vector(atlasRegions$atlasRegions)
  
  #Basically all atlas regions (everything that is within the brain)
  atlasRegionsBiggerZero<-atlasRegions[atlasRegions>0]
  
  print("Load humanBrain")
  humanBrain<-readMat("../storage//humanBrain_MNI_talairach.mat")
  humanBrainOffset<-humanBrain$offset
  humanBrain<-humanBrain$atlas
  
  
  
  #Create region-index for plotting
  uniqueRegions<-c()
  
  for (x in 1:11){
    
    fileName <- paste0("../storage/slice_svg/",x*9,".svg")
    
    xmlfile=xmlParse(fileName)
    xmltop = xmlRoot(xmlfile)
    
    uniqueRegions<-unique(c(uniqueRegions,recursiveGetAllStructures(xmltop[['g']])))
  }
  
  uniqueRegionAtlas<-list()
  if(file.exists("../storage/uniqueRegionAtlas.rds")){
    uniqueRegionAtlas<-readRDS("../storage/uniqueRegionAtlas.rds")
  }else{
    for(actRegion in 1:length(uniqueRegions)){
      uniqueRegionAtlas[[actRegion]]<-getSmallestParentAtlasRegions(atlasRegions,uniqueRegions[actRegion])
    }
    
    saveRDS(uniqueRegionAtlas,"../storage/uniqueRegionAtlas.rds")
  }
  
  
  
  
  all_genes<-readMat(allGenesData)
  expressionOfRandomGenes<-all_genes$expressionMatrix
  
  #REMOVE THIS FOR NORMAL DATA, ONLY IMPORTANT FOR EVO PAPER!
  ########################################################################
   expressionOfRandomGenes<-expressionOfRandomGenes[,match(read.csv2("../storage/humanEntrezIDs.csv")[,2],unlist(all_genes$entrezIDsOfRows))]
   all_genes$entrezIDsOfRows<-all_genes$entrezIDsOfRows[match(read.csv2("../storage/humanEntrezIDs.csv")[,2],unlist(all_genes$entrezIDsOfRows))]
  ########################################################################
  
   
   meanExprsForInjSiteNorm<-apply(expressionOfRandomGenes,1,function(x){median(x,na.rm=TRUE)})
   sdExprsForInjSiteNorm<-apply(expressionOfRandomGenes,1,function(x){mad(x,na.rm=TRUE)})
   
   expressionOfRandomGenes<-(expressionOfRandomGenes-meanExprsForInjSiteNorm)/sdExprsForInjSiteNorm
   expressionOfRandomGenes<-expressionOfRandomGenes[atlasRegions>0,]
   
   print(paste0("Mean mean of injection sites: ",mean(apply(expressionOfRandomGenes,1,function(x){mean(x,na.rm=TRUE)}))))
   print(paste0("Mean sd of injection sites: ",mean(apply(expressionOfRandomGenes,1,function(x){sd(x,na.rm=TRUE)}))))
   
   
   
   #Stanardize every random gene by its mean and standard deviation
   expressionOfRandomGenes<-apply(expressionOfRandomGenes,2,function(x){(x-median(x,na.rm=TRUE))/mad(x,na.rm=TRUE)})
   
   
   meanExprs<-apply(expressionOfRandomGenes,1,function(x){median(x,na.rm=TRUE)}) #mean gene expression for every grid point
   sdExprs<-apply(expressionOfRandomGenes,1,function(x){mad(x,na.rm=TRUE)})     #standard deviation of gene expression for every grid point
   
   
   amountOfSets<-length(setnames)   #Amount of gene-sets in testSetName.
   
   geneExpressionSynergyMatrix<-matrix(0,nrow=sum(atlasRegions>0),ncol = amountOfSets) #stores gene-expression synergy (grid-point-wise trimmed-mean gene-expression)
   ratiosMatrix<-list() #every gene in the gene-expression set can have a weighting ratio that gives them more weight in the gene-expression synergy , usually 1 for all
   amountOfGenesVector<-c()
   genenamesList<-list()
   geneExressionMatrix<-list()
   geneIndizesList<-list()
   
   #load all relevant data from all gene-expression sets and computes gene-expression synergies
   print("Load gene-expression sets and computes gene-expression synergies")
   start <- Sys.time ()
   dir.create(paste0(resultsFolder,'//',"GA_sets"))
   for(actSet in 1:amountOfSets){
     fname = setnames[actSet]
     print(paste0("",Sys.time()," ","Load: ",fname))
     
     #conData<-readMat(paste0(testSetName,"/",fname,".mat"))
     conData<-c()
     conData$genenames<-setFile[startsOfSets[actSet]:endsOfSets[actSet],1]
     conData$entrez<-setFile[startsOfSets[actSet]:endsOfSets[actSet],2]
     if(ncol(setFile)==3){
       conData$ratios<-setFile[startsOfSets[actSet]:endsOfSets[actSet],3] 
     }
     conData$expressionIndex<-(1:length(conData$genenames))[is.element(conData$entrez,unlist(all_genes$entrezIDsOfRows))]
     conData$expressionMatrix<-matrix(0,sum(atlasRegions>0),length(conData$genenames))
     conData$expressionMatrix[,conData$expressionIndex]<-expressionOfRandomGenes[,match(conData$entrez[conData$expressionIndex],unlist(all_genes$entrezIDsOfRows))]
     
     geneIndizesList[[actSet]]<-match(conData$entrez[conData$expressionIndex],unlist(all_genes$entrezIDsOfRows))
     
     expressionIndex<-as.vector(conData$expressionIndex)
     
     
     
     expressionOfGenes<-conData$expressionMatrix
     expressionOfGenes<-as.matrix(expressionOfGenes)
     
     originalSize<-length(expressionIndex)
     
     expressionOfGenes[expressionOfGenes==Inf]<-NA
     expressionOfGenes[expressionOfGenes==-Inf]<-NA

     
     if(length(expressionIndex)>1){
       #Perform gene set selection with genetic algorithm
       if(runGeneticAlgorithm){
         if(!file.exists(paste0(resultsFolder,'//','GA_sets//',fname,'.csv'))){
           expressionOfGenesRegions <-matrix(0,nrow=length(uniqueRegions),ncol=ncol(expressionOfGenes))
           
           for(actRegion in 1:length(uniqueRegions)){
             if(sum(atlasRegions==uniqueRegions[actRegion])==1){
               expressionOfGenesRegions[actRegion,]<-expressionOfGenes[atlasRegions==uniqueRegions[actRegion],]
             }else{
               expressionOfGenesRegions[actRegion,]<-apply(expressionOfGenes[atlasRegions==uniqueRegions[actRegion],],2,mean)
             }
           }
           
           
           repCount<-min((1:ceiling(1/trimExpression))[(length(expressionIndex)*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
           geneExpressionSynergyOrgRegions<-apply(expressionOfGenesRegions[,expressionIndex],1,function(x){mean(rep(x, repCount),trim=trimExpression,na.rm=TRUE)})
           
           
           GAmodel = rbga.int(size=length(expressionIndex), ints=expressionIndex, popSize=800, iters=500, mutationChance = 1/length(expressionIndex), evalFunc = function(y) {
             if(length(unique(y))<=1){
               return(Inf)
             }
             
             repCount<-min((1:ceiling(1/trimExpression))[(length(unique(y))*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
             aaa<-apply(expressionOfGenesRegions[,unique(y)],1,function(x){mean(rep(x, repCount),trim=trimExpression,na.rm=TRUE)})
             
             if(cor(geneExpressionSynergyOrgRegions,aaa,method="spearman",use="pairwise.complete.obs")<0.99){
               return(Inf)
             }
             
             return(-mad(aaa,na.rm=TRUE)*sqrt(length(unique(y))))
           },monitorFunc=function(x){
             if(x$iter==1 || (x$iter%%10)==0){
               lengthOfBestSet<-0
               if(sum(x$evaluations==min(x$evaluations,na.rm=TRUE))==1){
                 lengthOfBestSet<-length(unique(x$population[x$evaluations==min(x$evaluations),]))
                 print(paste0(x$iter,': ',-x$best[x$iter]," length: ",lengthOfBestSet))
               }else{
                 lengthOfBestSet<-min(apply(x$population[x$evaluations==min(x$evaluations),],1,function(y){length(unique(y))}))
                 print(paste0(x$iter,': ',-x$best[x$iter]," length: ",lengthOfBestSet))
               }
               
             }
           },useCores=useCores,equalsIterationsForConvergence=60)
           
           lengthOfBestSet<-0
           newExpressionIndex<-c()
           if(sum(GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE))==1){
             newExpressionIndex<-sort(unique(GAmodel$population[!is.na(GAmodel$evaluations)][GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE)]))
             lengthOfBestSet<-length(newExpressionIndex)
           }else{
             bestSets<-GAmodel$population[!is.na(GAmodel$evaluations),][GAmodel$evaluations[!is.na(GAmodel$evaluations)]==min(GAmodel$evaluations,na.rm=TRUE),]
             
             lengthOfBestSet<-min(apply(bestSets,1,function(y){length(unique(y))}))
             newExpressionIndex<-sort(unique(as.matrix(t(bestSets)[,lengthOfBestSet==apply(bestSets,1,function(y){length(unique(y))})])[,1]))
           }
           if(min(GAmodel$evaluations,na.rm=TRUE)<(-mad(geneExpressionSynergyOrgRegions,na.rm=TRUE)*sqrt(length(expressionIndex)))){
             print(paste0("Finished after ",GAmodel$iters,' iterations: ',-GAmodel$best[GAmodel$iters]," length: ",lengthOfBestSet))
             expressionIndex<-newExpressionIndex
           }else{
             print(paste0("Finished after ",GAmodel$iters,' iterations: ',-GAmodel$best[GAmodel$iters]," length: ",length(expressionIndex)))
           }
           
           if(length(expressionIndex)>1){
             write.table(data.frame(expressionindex=expressionIndex,genenames=unlist(conData$genenames[expressionIndex]),entrez=unlist(conData$entrez[expressionIndex])),file=paste0(resultsFolder,'//',"GA_sets//",fname,".csv"),sep = ",",row.names = FALSE)
           }
         }else{
           expressionIndex<-read.csv2(paste0(resultsFolder,'//','GA_sets//',fname,'.csv'),sep=",")$expressionindex
         }
       }
       
       print(paste0("Total: ",originalSize," left: ",length(expressionIndex)))
       expressionOfGenes<-expressionOfGenes[,expressionIndex]
     }
     
     if(length(expressionIndex)>0){
       genenames<-as.vector(conData$genenames[expressionIndex])
       genenamesList[[actSet]]<-genenames
       print(unlist(genenames))
       
       ratios<-conData$ratios
       if(length(ratios)==0){
         ratios<-rep(1,length(expressionIndex))
       }else{
         ratios<-ratios[expressionIndex]
       }
       
       ratiosMatrix[[actSet]]<-ratios
       
       expressionOfGenes<-t(t(expressionOfGenes)*ratios)
       amountOfGenes<-dim(expressionOfGenes)[2]
       
       geneExressionMatrix[[actSet]]<-expressionOfGenes
       amountOfGenesVector[actSet]<-amountOfGenes
       
       #calculate gene expression synergy (trimmed mean)
       if(amountOfGenes>1){
         repCount<-min((1:ceiling(1/trimExpression))[(amountOfGenes*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
         geneExpressionSynergy<-apply(expressionOfGenes,1,function(x){mean(rep(x, repCount),trim=trimExpression,na.rm=TRUE)})
         #the repetition (rep) of genes is necessary, since a 5% lower and upper cutoff would only work for at least 20 genes.
         #example: one could not compute the 10% trimmed mean of 5 numbers [0,1,2,3,10], but taking those 5 numbers 2 times [0,1,2,3,10,0,1,2,3,10]
         #its possible (it doesn't trim the "10" perfectly, but still better than without)
       }else{
         geneExpressionSynergy<-expressionOfGenes
       }
       
       geneExpressionSynergy[is.na(geneExpressionSynergy)]<-0
       
       geneExpressionSynergy<-geneExpressionSynergy
       geneExpressionSynergyMatrix[,actSet]<-geneExpressionSynergy
       
     }else{
       geneExressionMatrix[[actSet]]<-0
       amountOfGenesVector[actSet]<-0
       ratiosMatrix[[actSet]]<-0
       genenamesList[[actSet]]<-0
       geneIndizesList[[actSet]]<-0
     }
     
     
     
   }
   print(paste0("Finished after: ",difftime(Sys.time(),start,units="mins"),"min"))
   
   
   for(actSet in c(1:amountOfSets)[amountOfGenesVector>0]){
     fname = setnames[actSet]
     
     print(paste0("",Sys.time()," ","Start pvalue calculation for: ",fname," (",actSet,"/",amountOfSets,")"))
     
     
     dir.create(paste0(resultsFolder,'//',paste0(fname)))
     dir.create(paste0(resultsFolder,'//',"pvalues"))
    
     
     start <- Sys.time ()
     
     #get current gene set (actSet) specific data
     geneExpressionSynergy<-geneExpressionSynergyMatrix[,actSet]
     expressionOfGenes<-geneExressionMatrix[[actSet]]
     amountOfGenes<-amountOfGenesVector[actSet]
     genenames<-genenamesList[[actSet]]
     ratios<-ratiosMatrix[[actSet]]
     geneIndizes<-geneIndizesList[[actSet]]
     
     #if p-values havent been computed yet, compute them now
     if(!file.exists(paste0(resultsFolder,'//','pvalues/',fname,'_ttest_geneExpr.rds'))){
       
       
       #those are the mean/standard deviation of gene-expression, computed
       #from random-gene set synergies
       geneExpressionMeans<-rep(0,sum(atlasRegions>0))
       geneExpressionSDs<-rep(0,sum(atlasRegions>0))
       
       
       computationStart<-Sys.time()
       for(actBatch in 1:splitInBatches){
         
         inBetweenStart<-Sys.time()
         
         
         chunkSize<-min(100,ceiling(dim(connectivityMatrix)[2]/useCores)) #splits the connectivityMatrix into chunks column wise, to fit the computation into the memory
         
         if(useCores>1){
           print("Setup parallel computing...")
           clComputing <- makeCluster(useCores,outfile="parallel_log.txt",type="PSOCK")
           registerDoParallel(clComputing)
           print("Parallel computing setup done")
         }
         
         
         print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5)),": Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Generate ",ceiling(precision/splitInBatches)," random gene-sets..."))
         
         #Generates an amount (precision/splitInBatches) of random gene-sets and computes the gene expression synergy (trimmed-mean)
         expressionOfRandomGenesRowSumsMatrix<-foreach(o=1:ceiling(precision/splitInBatches),.combine="cbind",.multicombine=TRUE, .inorder=FALSE) %op%{
           randnum<-sample((1:dim(expressionOfRandomGenes)[2])[-geneIndizes],amountOfGenes)
           expressionOfRandomGenesRowSums<-0
           if(amountOfGenes>1){
             repCount<-min((1:ceiling(1/trimExpression))[(amountOfGenes*(1:ceiling(1/trimExpression)))%%ceiling(1/trimExpression)==0])
             expressionOfRandomGenesRowSums<-apply(t(expressionOfRandomGenes[,randnum])*ratios,2,function(x){mean(rep(x,repCount),trim=trimExpression,na.rm=TRUE)})
             #the repetition (rep) of genes is necessary, since a 5% lower and upper cutoff would only work for at least 20 genes.
             #example: one could not compute the 10% trimmed mean of 5 numbers [0,1,2,3,10], but taking those 5 numbers 2 times [0,1,2,3,10,0,1,2,3,10]
             #its possible (it doesn't trim the "10" perfectly, but still better than without)
           }else{
             expressionOfRandomGenesRowSums<-expressionOfRandomGenes[,randnum]
           }
           
           return(expressionOfRandomGenesRowSums)
         }
         print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5))," Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Random gene-sets generated after ",difftime(Sys.time(),inBetweenStart,units="mins"),"min"))
         
         geneExpressionMeans<-geneExpressionMeans+apply(expressionOfRandomGenesRowSumsMatrix,1,function(x){mean(x,na.rm=TRUE)})
         geneExpressionSDs<-geneExpressionSDs+apply(expressionOfRandomGenesRowSumsMatrix,1,function(x){sd(x,na.rm=TRUE)})
         
         expressionOfRandomGenesRowSumsMatrix[is.na(expressionOfRandomGenesRowSumsMatrix)]<-0
         
         
         
         if(useCores>1){
           print("Parallel computing finished")
           stopCluster(clComputing)
         }
         
         rm(expressionOfRandomGenesRowSumsMatrix)
         gc()
         
         print(paste0("",Sys.time()," ",actSet,":",substr(fname,1,min(nchar(fname),5))," Batch ",(((actBatch-1)*100)/splitInBatches),"%-",((actBatch*100)/splitInBatches),"% - Node-strengths computed after ",difftime(Sys.time(),inBetweenStart,units="mins"),"min"))
         
         
       }
       
       print("Calculate p-values...")
       pvalsExpr_ttest<-(1-pnorm((geneExpressionSynergy),mean=geneExpressionMeans/splitInBatches,sd=geneExpressionSDs/splitInBatches))
       pvalsExpr_ttest[geneExpressionSDs==0]<-NA
       
       print("p-values calculated!")
       saveRDS(pvalsExpr_ttest, file = paste0(resultsFolder,'//','pvalues/',fname,'_ttest_geneExpr.rds'))
       
     }
     
     #Load calculated p-values from hard-drive
     subfolder<-""
     print("TTest: ")
     pvalsExpr<-readRDS(paste0(resultsFolder,'//','pvalues//',fname,'_ttest_geneExpr.rds'))
     subfolder<-"//TTest"
     
     #Calculate FDR (q-values)
     qvalsExpr<-rep(1,length(pvalsExpr))
     qvalsExpr[!is.na(pvalsExpr)]<-p.adjust(pvalsExpr[!is.na(pvalsExpr)], method="BH")
     
     dir.create(paste0(resultsFolder,'//',fname,subfolder))
    
     pvalsExpr[is.na(pvalsExpr)]<-1
     qvalsExpr[is.na(qvalsExpr)]<-1
     
     print("Percentage of significant regions in the brain: ")
     print(paste0("First order FDR(0.05): ",sum(qvalsExpr<=0.05)/length(qvalsExpr)))
     print(paste0("First order FDR(0.1): ",sum(qvalsExpr<=0.1)/length(qvalsExpr)))
     
     print("Plotting results (MRI style plot)...")
     plot_mri_style_image(paste0(resultsFolder,'//',fname,subfolder,'/mri_style_pvals_max10_firstOrder.png'),pvalsExpr,max(pvalsExpr[qvalsExpr<=0.05]),max(pvalsExpr[qvalsExpr<=0.1]))
     
     
   }
}