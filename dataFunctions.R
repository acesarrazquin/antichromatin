################################################################################
# createModList creates a list of modifications with repetition 
# (if several modifications in same peptide) from the peptide list.
#
# The peptide list should contain at least the fields:
#   - protein name (H2, H3, p53...) Here 'Protein'
#   - protein modification: whole mod.whole e.g. K12me1S41ph Here 'Modification'
# 
################################################################################

createModList <- function(peptide.list){
  mod.list <- peptide.list[, c('Protein', 'Modification', 'Fragment', 'ID')]
  result.mod.list <- data.frame()
  
  uniq.prots <- sort(as.character(unique(mod.list$Protein)))
  
  # iterate on each protein
  for (prot in uniq.prots){
    prot.mod <- mod.list[mod.list$Protein == prot, ]
    new.prot.mod <- data.frame()
    
    # iterate on each whole modification 
    for (i in (1:nrow(prot.mod))){ 
      mod.whole <- as.character(prot.mod[i, 'Modification'])
      mod.regex <- '[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*'
      pept_id <- as.numeric(prot.mod[i, 'ID'])
      print(mod.whole)
      
      if (!grepl(mod.regex, mod.whole)){
        name <- paste(prot, mod.whole, sep=' ')
        rest <- prot.mod[i, "Fragment"]
        new.prot.mod <- rbind(new.prot.mod, data.frame(name=name, pept_id=pept_id, protein=prot, modification="",
                                                       nb=0, residue.first="", ptm.first="",
                                                       mod.rest=rest))
      }else{
        name <- paste(prot, mod.whole, sep=' ')
        # separate if several modifications
        mods <- regmatches(mod.whole, gregexpr(mod.regex, mod.whole))[[1]]
        
        # iterate on each simple modification
        for (mod in mods){ 
          residue <- regmatches(mod, gregexpr('[A-Z]{1}[0-9]+', mod))[[1]]
          nb <- as.integer(substring(residue, 2)) # residue number
          
          ptm <- regmatches(mod, gregexpr('[a-z]{2}[0-9]*[s|u]*', mod))[[1]]
          rest <- sub(mod, '', mod.whole)
          
          new.prot.mod <- rbind(new.prot.mod, data.frame(name = name, pept_id=pept_id, protein=prot, modification=mod.whole,
                                                         nb=nb, residue.first=residue, ptm.first=ptm,
                                                         mod.rest=rest))
        }                                               
      }
    }

    #order by residue number
  
    order.index <- order(new.prot.mod$nb)
    new.prot.mod <- new.prot.mod[order.index,]

    #order by modification
    unique.res <- as.character(unique(new.prot.mod$residue.first))

    new.prot.mod.ord <- data.frame()
    for (res in unique.res){
      res.mod <- new.prot.mod[new.prot.mod$residue.first==res, ]
       
      order.index2 <- order(as.character(paste(res.mod$ptm.first, res.mod$mod.rest)))
      res.mod <- res.mod[order.index2, ]
      
      
      new.prot.mod.ord <- rbind(new.prot.mod.ord, res.mod)
    }
    
    
    result.mod.list <- rbind(result.mod.list, new.prot.mod.ord)
  }

  rownames(result.mod.list) <- 1:nrow(result.mod.list)
  result.mod.list$nr <- factor(rownames(result.mod.list), 
                               levels=rownames(result.mod.list), ordered=TRUE) # so that the level order is correct; important for plots
  return(result.mod.list[, c(9,2,1,3,5,6,7,8)]) # exclude modification
}

#######################################################################################
# annotatePrimaryTargets annotates, in the delfia results table, the peptides that are
# supposed to be target of the antibody they are tested with.
#
# Input: ab table, peptide table, and delfia results. Necessary fields:
#   - Peptide_ID and Antibody_ID in delfia, and corresponding ID in ab and pepts table
#   - ab table: Modification and Protein
#   - pept table: Modification and Protein
# Output: reduced delfia table, with "Primary_Target" column
#
# Mode option:
#   - "simple" mode (default): in the case of antibodies against complex modifications 
#     (e.g.K9acS10ph) is considered target any peptide with any of the individual ptm
#   - "complex" mode: targets are only the peptides that bear the complex ptm
#
# Some other considerations:
#   - names of antibodies with exclusive expecificities are separated by / 
#   - general antibodies (against H3, or K-) are considered as special cases
#        - for abs against any modification of a histone, I substitute the
#         "-" symbol by a regexp "[0-9]" -> DOESN'T WORK???
#   - the output is a reduced table from the merged three tables - choose fields
#
# ISSUES: general antibodies, K4unmod antibody
#######################################################################################
annotatePrimaryTargets <- function(delfia, pepts, abs, mode="simple"){

  # merge the three tables to get the modifications
  delf <- merge(delfia, pepts, by.x="Peptide_ID", by.y="ID")
  ##delf <- merge(delf, abs, by.x="Antibody_ID", by.y="ID")
  delf <- merge(delf, abs, by="Barcode")
 
  #regex <- '[A-Z]{1}[0-9-]+[a-z-]{2}[0-9]*[u|s]*' #put "-" for case where no fixed residue or modification
  regex <- '([A-Z]{1}[0-9-]+)?([a-z-]{2}[0-9]*[u|s]*)?'
  prim.tars <- vector() 
  for (i in 1:nrow(delf)){
    
    # select modification in antibody
    ab.modifs <- as.character(delf$Modification.y[i]) #mods separated if /
    ab.prot <- as.character(delf$Protein.y[i])
    ab.frag <- delf$Fragment.y[i]
      
    # select modification in peptide

    pep.frag <- substr(delf$Fragment.x[i], 2, nchar(as.character(delf$Fragment.x[i]))-1)
    pep.mod <- as.character(delf$Modification.x[i])
    pep.mods <- regmatches(pep.mod, gregexpr(regex, pep.mod))[[1]] # mods in peptide
    pep.prot <- as.character(delf$Protein.x[i])
    
    # include primary_target (yes/no) column
    prim.tar <- "no"
    if (ab.prot==pep.prot | ab.prot=="-" | grepl(ab.prot, pep.prot)){
      if (ab.modifs == ""){ # or change by "-"
        if(ab.frag == ""){
          prim.tar <- "yes"
        }else{     
          ab.frag.split <- as.numeric(strsplit(ab.frag, split="-")[[1]])
          pep.frag.split <- as.numeric(strsplit(pep.frag, split="-")[[1]])
          if (pep.frag.split[1] >= ab.frag.split[1] & 
                pep.frag.split[1] <= ab.frag.split[2]){
            prim.tar <- "yes"
          }
        }
      }else{
#         if (mode=="complex"){#NOT UPDATED
#           ab.ex.mods <- strsplit(ab.modifs, '/')[[1]]
#           for (ab.ex.mod in ab.ex.mods){ # exclusive modifications: e.g. K4me1/K4ac
#             prim.tar <- 'yes'
#             ab.mods <- regmatches(ab.ex.mod, gregexpr(regex, ab.ex.mod))[[1]] # mods in antibody
#             
#             for (ab.mod in ab.mods){ # each modification in the antibody
#               # correction
#               ab.mod <- gsub("-", "[0-9]*", ab.mod)
#               
#               if (sum(grepl(ab.mod, pep.mods))==0){ # if the modifications not equal, then its not a primary target
#                 prim.tar <- 'no'
#               }
#             } 
#             
#             if (prim.tar=='yes'){
#               break
#             }
#           }
        if (mode=="simple"){ # could remove it if finally it is never "complex"
      
          ab.modifs <- gsub("/", "", ab.modifs) # we don't consider whole modifications
          ab.mods <- regmatches(ab.modifs, gregexpr(regex, ab.modifs))[[1]]
          
          for (ab.mod in ab.mods){ # each modification in the antibody   
            #different cases
            if (grepl("^[A-Z][0-9-]*$", ab.mod)){ #only residue, with or without number:K-, K4
              
              nb <- as.numeric(regmatches(ab.mod, gregexpr("[0-9-]+", ab.mod))[[1]]) # residue number, can be ""
              if (length(nb)==0){
                nb==0
              }
              
              frag.split <- strsplit(pep.frag, split="-")[[1]] # peptide fragment
              
              if(sum(grepl(paste(ab.mod, "[a-z-]{2}[0-9]*[u|s]*", sep=""), 
                           pep.mods))==0){ # check peptide is not modified on that residue
                if(nb == "-" | nb == 0 | 
                   nb >= as.numeric(frag.split[1]) & nb <= as.numeric(frag.split[2])){
                  prim.tar <- "yes"
                  break
                }
              }
            }else{ #case of only modification: ac, me2 me2u, me-, or complete: K4ac, K-ac...
              # correction
              ab.mod <- gsub("-", "[0-9]*", ab.mod)
              
              if (sum(grepl(ab.mod, pep.mods))!=0){
                prim.tar <- "yes"
                break
              }
            }
          }
        }
      }
    }
      
    prim.tars <- c(prim.tars, prim.tar) 
  }

  delf$Primary_Target <- factor(prim.tars) 
  delf <- delf[,c('Signal', 'Peptide_ID', 'Antibody_ID', 'Barcode', 'Signal_Mean', 'Signal_StdDev', 'Signal_Mean_noOutliers', 
                  'Signal_StdDev_noOutliers', 'OutlierWell', 'Primary_Target')]# select only some columnms
  return(delf)
}

#####################################################################################
# normalizeDelfia calculates summary statistics for replicates and normalizes signal
# values.
#
# The normalization can be done by:
#   - by="max": dividing by the maximum signal value for each antibody.
#     - issue if that value is an oulier?
#   - by="maxmedian": dividing by the maximum median signal value for each antibody
#     - this is the default
#
# Input: - delfia table
# Output: - modified delfia table, with statistical descriptors, raw and normalized
#####################################################################################
normalizeDelfia <- function(delfia, by="maxmedian"){
  # Identify antibodies by Barcode, and ID (just in case...)
  ab.unique <- unique(delfia[,c("Antibody_ID", "Barcode")])

  delfia.norm <- data.frame()
  for (i in 1:nrow(ab.unique)){   
    id <- as.numeric(ab.unique[i, 1])
    barcode <- as.character(ab.unique[i,2])
    
    print(paste(i, barcode, id))
    
    dataset <- delfia[delfia$Antibody_ID==id & delfia$Barcode==barcode,]
    
    delfia.ab <- data.frame()
    
    # calculate statistics for replicates
    pep.unique <- unique(as.numeric(dataset$Peptide_ID))
    for (pep in pep.unique){ # values for antibody and peptide
      dataset.pep <- dataset[dataset$Peptide_ID == pep,]
      signal <- as.numeric(dataset.pep$Signal)
      
      summary.signal <- summary(signal)
      
      min.signal <- summary.signal[[1]]
      q1.signal <- summary.signal[[2]]
      median.signal <- summary.signal[[3]]
      mean.signal <- summary.signal[[4]]
      q3.signal <- summary.signal[[5]]
      max.signal <- summary.signal[[5]]
      sd.signal <- sd(signal)
      
      outliers.val <- boxplot.stats(signal)$out
      outliers.index <- which(signal %in% outliers.val)
      dataset.pep$outlier <- "no"
      dataset.pep$outlier[outliers.index] <- "yes"
     
      
      delfia.ab.new <- cbind(dataset.pep, signal, min.signal, q1.signal, median.signal,
                         mean.signal, q3.signal, max.signal, sd.signal)
      delfia.ab <- rbind(delfia.ab, delfia.ab.new)
    }
    
    # normalization
    if (by=="max"){
      factor <- max(delfia.ab$signal)
    }else if (by=="maxmedian"){
      factor <- max(delfia.ab$median.signal)
    }
     
    delfia.ab$signal.norm <- delfia.ab$signal/factor
    #delfia.ab[delfia.ab$signal.norm > 1, "signal.norm"] <- 1 #in case values are higher than 1
    delfia.ab$q1.signal.norm <- delfia.ab$q1.signal/factor
    delfia.ab$median.signal.norm <- delfia.ab$median.signal/factor
    delfia.ab$mean.signal.norm <- delfia.ab$mean.signal/factor
    delfia.ab$q3.signal.norm <- delfia.ab$q3.signal/factor
    delfia.ab$max.signal.norm <- delfia.ab$max.signal/factor
    delfia.ab$sd.signal.norm <- delfia.ab$sd.signal/factor
    
    delfia.norm <- rbind(delfia.norm, delfia.ab)
  }

  return(delfia.norm)
}

#####################################################################################
# calculateThresholds calculates a threshold value for each antibody in order to
# distinguish between signal and background values.
#
# The calculation is done by applying kmeans with 2 clusters on the intensity values 
# each antibody, and considering as threshold the maximum value of the background cluster.
# K-means can be applied to:
#   -by="median": the median normalized signal values
#       - default
#   -by="values": the individual normalized signal values
# 
# Input:  - delfia data
#         - antibody list
# Output: - modified antibody list with "Threshold" column
#         - the list is extended to include different replicates (barcodes) for same ID
# Note: other functions for calculating thresholds at the end of this file.
#####################################################################################
calculateThresholds <- function(delfia, abs, by="median", method="kmeans"){
  ab.bars <- unique(abs$Barcode)

  ab.threshold <- data.frame()
  for (barcode in ab.bars){
    # perform k-means clustering with 2 centers (noise & detection)
 
    id = abs[abs$Barcode==barcode, "ID"]

    print(paste(barcode, id))
    
    if(by=="median"){
      dataset <- unique(delfia[delfia$Barcode==barcode, 
                               c("Barcode", "Antibody_ID", "Peptide_ID", "median.signal.norm")])#CONSIDER JUST MEDIANS
      signal <- dataset$median.signal.norm
    }else if(by=="values"){
      dataset <- delfia[delfia$Barcode==barcode,]
      signal <- dataset$signal.norm
    }

    km <- kmeans(signal, centers=2)
    noise.cluster <- names(sort(km$centers[,1])[1])
    noise.signals <- signal[which(km$cluster==noise.cluster)]
    if (method=="kmeans"){
      threshold <- max(noise.signals)
      #threshold <- min(0.25, threshold) # probably better to remove this, not good in many cases
    }else if(method=="quantiles"){
      threshold <- quantile(noise.signals, 0.99) # other options?
    }

    ab.threshold <- rbind(ab.threshold, data.frame(Barcode=barcode,
                                                   Threshold=threshold))
  }
  
  abs <- merge(abs, ab.threshold, by="Barcode")
    
  return(abs)
}
###############################################################################
# classifyAntibodies
###############################################################################
classifyAntibodies <- function(delfia, abs, correction=FALSE){
  ab.bars <- abs$Barcode
  ab.sps <- data.frame()
  for (barcode in ab.bars){
    id = abs[abs$Barcode==barcode, "ID"]
    
    print(paste(barcode, id))
    
    dataset <- unique(delfia[delfia$Barcode==barcode, 
                             c("Barcode", "Antibody_ID", "Peptide_ID", "median.signal.norm",
                               "Primary_Target")])
    
    # separate signal from background
    threshold <- abs[abs$Barcode==barcode, 'Threshold']

    signal.values <- dataset[dataset$median.signal.norm > threshold, ]
    background.values <- dataset[dataset$median.signal.norm <= threshold, ]

    # calculate FP, FN, TP, TN

    #TRY: 
    tp.signs <- signal.values[signal.values$Primary_Target=="yes", "median.signal.norm"]
    fp.signs <- signal.values[signal.values$Primary_Target=="no", "median.signal.norm"]
    fn.signs <- background.values[background.values$Primary_Target=="yes", "median.signal.norm"]
    tn.signs <- background.values[background.values$Primary_Target=="no", "median.signal.norm"]

    tp <- nrow(signal.values[signal.values$Primary_Target=="yes", ])
    fp <- nrow(signal.values)-tp
    fn <- nrow(background.values[background.values$Primary_Target=="yes", ])
    tn <- nrow(background.values)-fn
    
    # calculate sp and s
    sp <- tn/(tn+fp)
    s <- tp/(tp+fn)
    ppv <- tp/(tp+fp)
    npv <- tn/(tn+fn)
    
    ppv.adj <- ppv
    if (correction==TRUE){
      if(length(fp.signs) > 0 && length(tp.signs)>0){
        ppv.adj <- ppv + 0.5 * (median(tp.signs) - median(fp.signs))
      }
    }

    # classify in categories
    if (is.nan(ppv.adj) | is.nan(s)){
      category=0
    }else{
      if(ppv.adj>0.75){
        if(s>0.75){
          category=1
        }else if(s<=0.75 & s>0.5){
          category=2
        }else if(s<=0.5){
          category=3
        }
      }else if(ppv.adj>0.5 & ppv.adj<=0.75){
        if(s>0.75){
          category=4
        }else if(s<=0.75 & s>0.5){
          category=5
        }else if(s<=0.5){
          category=6
        }
      }else{
        if(s>0.75){
          category=7
        }else if(s<=0.75 & s>0.5){
          category=8
        }else if(s<=0.5){
          category=9
        }
      }
    }
    
    # fill data frame
    ab.sps <- rbind(ab.sps, data.frame(Barcode=barcode, 
                                       SP=sp,SE=s,PPV=ppv, NPV=npv, PPV.adj=ppv.adj,
                                       Class=category))
  }
  abs <- merge(abs, ab.sps, by="Barcode")
  return(abs)
}


####################################################################################
# Annotate with a replicate number the antibodies that have the same Name, Company
# and CatNr.
#
# In some cases these are tested with different dilutions, or just tested two or more
# times. Some have different ID (i.e. for different dilutions) and some don't (true)
# replicate.
#
# Input: - antibody list (modified with barcodes!!)
# Ouput: - antibody list with an extra "Replicate" column
####################################################################################
annotateRepeatedAb <- function(dataset){
  # Antibodies are the same if they have same ID
  # Barcode is what distinguishes replicates
  # Replicates if they have same ID, but also Dilution and LotNr
  
  antibody.uniqs <- unique(dataset[, c("ID", "LotNr", "Dilution")])
  
  rep.table <- data.frame() 
  for (i in 1:nrow(antibody.uniqs)){   
    id <- antibody.uniqs[i, "ID"]
    lotnr <- antibody.uniqs[i, "LotNr"]
    dil <- antibody.uniqs[i, "Dilution"]
    
    barcodes <- dataset[dataset$ID==id & dataset$LotNr==lotnr & dataset$Dilution==dil,
                             "Barcode"]
    
    print(id)
    print(paste("barcodes:", barcodes))
    
    if (length(barcodes)==1){
      rep.table <- rbind(rep.table, data.frame(Barcode=barcodes, 
                                               Replicate=""))
    }else{
      for (j in 1:length(barcodes)){
        rep.table <- rbind(rep.table, data.frame(Barcode=barcodes[j], 
                                                 Replicate=as.character(j)))
      }
    }
  }
  
  #print(rep.table)
  
  dataset <- merge(dataset, rep.table, by="Barcode")
  return(dataset)
}

############################################################################################################
# orderMods orders the lists of antibodies/peptides according to their modifications.
#
# Needed: column "Modification"
#
###########################################################################################################

orderMods <- function(data, type="abs", level=1){
  mod.regex <- "[A-Z]{1}[0-9]+([a-z]{2}[0-9]*[s|u]*)*" # parenthesis to catch also cases with no modification (abs e.g. K4, K9)
  sorted.ids <- vector()
  
  mods <- data$Modification
  ids <- row.names(data)
  
  level.mods <- sapply(mods, function(x){regmatches(x, gregexpr(mod.regex, x))[[1]][level]}) #level modification
  names(level.mods) <- ids
  
  # NAs (no modification - with that pattern): order by fragment (peptides) or no order (antibodies)
  na.mods <- level.mods[is.na(level.mods)]
  if (length(na.mods) != 0){
    frags <- as.character(data[row.names(data) %in% names(na.mods), "Fragment"])
    
    if (type=="pepts" & length(frags) !=0){
      first.aa <- sapply(frags, function(x){as.numeric(sub("\\(", "", strsplit(x, split="-")[[1]][1]))})
      names(first.aa) <- names(na.mods) # because the operation removes the names...
      
      sorted.ids <- c(sorted.ids, names(sort(first.aa)))
      
      if (length(na.mods) == length(level.mods)){ # no more levels
        return(sorted.ids)
      }
      
    }else{
      sorted.ids <- c(sorted.ids, names(na.mods))
    }
  }
  
  # modifications (no NA): order by residue nb, residue, ptm
  full.mods <- level.mods[!is.na(level.mods)]
  if (length(full.mods != 0)){
    res <- sapply(full.mods, function(x){regmatches(x, gregexpr("[A-Z]{1}[0-9]+", x))[[1]][1]})
    aas <- sapply(res, function(x){substring(x, 1, 1)})
    nbs <- sapply(res, function(x){as.numeric(substring(x, 2, nchar(x)))}) #when i convert to numeric the names disappear
    names(nbs) <- names(aas)
    
    ptms <- sapply(full.mods, function(x){regmatches(x, gregexpr("[a-z]{2}[0-9]*[s|u]*", x))[[1]][1]})
    
    for (nb in sort(unique(nbs))){
      ids <- names(nbs[nbs==nb])
      aas.nb <- aas[names(aas) %in% ids]
      
      for (aa in sort(unique(aas.nb))){
        ids <- names(aas.nb[aas.nb==aa])
        ptms.aa <- ptms[names(ptms) %in% ids]
        
        for (ptm in sort(unique(ptms.aa), na.last=FALSE)){ # na.last=FALSE to put NAs at beggining of list
          if (is.na(ptm)){ # if NA: don't sort
            ids <- names(ptms.aa[is.na(ptms.aa)])
            sorted.ids <- c(sorted.ids, ids)
          }else{
            ids <- names(ptms.aa[ptms.aa==ptm])
            sorted.ids <- c(sorted.ids, orderMods(data[row.names(data) %in% ids, ], level=level+1))
          }
        }
      }
    } 
  } 
  return(sorted.ids)
}

orderElms <- function(elms, type="abs"){  
  uniq.prots <- sort(as.character(unique(elms$Protein)))
  sorted.ids <- vector()
  
  for (prot in uniq.prots){
    print(prot)
    
    prot.data <- elms[elms$Protein == prot, ]
    
    sorted.ids <-  c(sorted.ids, orderMods(prot.data, type=type))
  } 
  elms.ord <- elms[sorted.ids,]
  
  elms.ord$ID <- as.character(elms.ord$ID)
  
  elms.ord[sapply(elms.ord, is.character)] <- lapply(elms.ord[sapply(elms.ord, is.character)], 
                                                     function(x){factor(x, levels=unique(x))})
  
  return(elms.ord)   
}







####OLD#####################################
annotateRepeatedPep <- function(dataset){
  rep.table <- data.frame()
  pep.catnr.name <- unique(paste(dataset$CatNr, dataset$Name, sep="$$"))
  for (i in pep.catnr.name){
    pep.catnr <- strsplit(i, split="\\$\\$" )[[1]][1]
    pep.name <- strsplit(i, split="\\$\\$" )[[1]][2]
    nb.pep <- unique(dataset[dataset$CatNr==pep.catnr & dataset$Name==pep.name,
                             'ID'])
    if (length(nb.pep)==1){
      rep.table <- rbind(rep.table, data.frame(ID=nb.pep, Replicate=""))
    }else{
      for (j in 1:length(nb.pep)){
        rep.table <- rbind(rep.table, data.frame(ID=nb.pep[j], Replicate=as.character(j)))
      }
    }
  }
  print(rep.table)
  
  dataset.new <- merge(dataset, rep.table, by.x="ID", by.y="ID")
}


#### More options for threshold calculation####
calculateThresholds.ckmeans <- function(delfia, abs){ # gets the same as k-means
  ab.unique <- unique(as.numeric(delfia$Antibody_ID))
  thresholds <- vector()
  for (ab in ab.unique){
    # perform k-means clustering with 2ls centers (noise & detection)
    print(ab)
    dataset <- delfia[delfia$Antibody_ID==ab,]
    ckm <- Ckmeans.1d.dp(dataset$signal.norm, k=2)
    
    threshold <- max(dataset$signal.norm[which(ckm$cluster==1)])
    
    thresholds <- c(thresholds, threshold)
  }
  ab.threshold <- data.frame(ab.id=ab.unique, thres.ckm=thresholds)
  
  abs <- merge(abs, ab.threshold, by.x="ID", by.y="ab.id")
  
  return(abs)
}

calculateThresholds.outs <- function(delfia, abs){ 
  ab.unique <- unique(as.numeric(delfia$Antibody_ID))
  thresholds <- vector()
  for (ab in ab.unique){
    print(ab)
    dataset <- delfia[delfia$Antibody_ID==ab,]
    outliers <- boxplot.stats(dataset$signal.norm)$out
    
    noout <- dataset[!(dataset$signal.norm %in% outliers),]
    threshold <- quantile(noout$signal.norm, 0.99)
    
    thresholds <- c(thresholds, threshold)
  }
  ab.threshold <- data.frame(ab.id=ab.unique, thres.out=thresholds)
  
  abs <- merge(abs, ab.threshold, by.x="ID", by.y="ab.id")
  
  return(abs)
}

calculateThresholds.tan <- function(delfia, abs){
  ab.unique <- unique(as.numeric(delfia$Antibody_ID))
  thresholds <- vector()
  for (ab in ab.unique){
    
    print(ab)
    dataset <- delfia[delfia$Antibody_ID==ab,]
    
    quantiles <- quantile(dataset$signal.norm, 1:100*0.01)
    angles <- vector()
    for (i in 1:length(quantiles)-1){
      num <- quantiles[i+1]-quantiles[i]

      angle <- asin(num/1)*180/pi
      angles <- c(angles, angle)
    }

    angles.rest <- 45-angles
    min.rest <- min(angles.rest)
    min.rest.index <- which(angles.rest==min.rest)
    threshold <- quantiles[min.rest.index]
    
    thresholds <- c(thresholds, threshold)
    print(threshold)
  }
  ab.threshold <- data.frame(ab.id=ab.unique, thres.tan=thresholds)
    
  abs <- merge(abs, ab.threshold, by.x="ID", by.y="ab.id")
    
  return(abs)
}
    

calculateThresholds.mclust <- function(delfia, abs){ # mclust - it doesn't perform a good separation...
  ab.unique <- unique(as.numeric(delfia$Antibody_ID))
  thresholds <- vector()
  for (ab in ab.unique){
    # perform k-means clustering with 2 centers (noise & detection)
    print(ab)
    dataset <- delfia[delfia$Antibody_ID==ab,]
    
    #noise.size <- nrow(dataset)/4
    #noise.interval <- which(dataset$signal.norm %in% dataset[dataset$signal.norm>median(dataset$signal.norm) & 
    #                 dataset$signal.norm<quantile(dataset$signal.norm)[4],
    #                 "signal.norm"])
    #noise <- which(dataset$signal.norm %in% sample(dataset$signal.norm, noise.size))
    data.log <- log(dataset$signal.norm)               
    clusters <- mclustBIC(data=data.log, G=2)# 
                       #initialization=list(noise=noise))
    clusters.summary <- summary(data=data.log, clusters)
    classif <- clusters.summary$classification
    
    #noise.cluster <- which(clusters$parameters$mean==min(clusters$parameters$mean))

    threshold <- max(dataset$signal.norm[which(classif=="1")])
  
    thresholds <- c(thresholds, threshold)
  }
  ab.threshold <- data.frame(ab.id=ab.unique, thres.mclust2=thresholds)
  
  abs <- merge(abs, ab.threshold, by.x="ID", by.y="ab.id")
  
  return(abs)
}



# normalizeDelfia <- function(delfia){
#   ab.unique <- unique(as.numeric(delfia$Antibody_ID))
#   max.signal <- max(delfia$Signal) # OR saturation of the machine
#   delfia.norm <- data.frame()
#   for (ab in ab.unique){
#     # perform k-means clustering with 2 centers (noise & detection)
#     print(ab)
#     dataset <- delfia[delfia$Antibody_ID==ab,]
#     km <- kmeans(dataset$Signal, centers=2)
#     noise.avg <- min(km$centers)
#     #max.signal <- max(dataset$Signal) - noise.avg #doesn't make sense to compare different abs...
#     
#     dataset$Signal.norm <- (dataset$Signal - noise.avg)/max.signal
#     #dataset[dataset$Signal.norm < 0,]$Signal.norm  <- 0 #make 0 the negative values?
#     
#     pep.unique <- unique(as.numeric(dataset$Peptide_ID))
#     for (pep in pep.unique){
#       dataset.pep <- dataset[dataset$Peptide_ID == pep,]
#       mean.signal <- mean(as.numeric(dataset.pep$Signal.norm))
#       dataset.pep$Signal.norm_allReplicates_Mean <- mean.signal
#       
#       delfia.norm <- rbind(delfia.norm, dataset.pep)
#     }
#     
#   }
#   return(delfia.norm)
# }
# antibody specificity? penalize strongly the nonspecific hits up to half of the others
# score <- function(tp,fn,fp){
#   +     sp.found <- tp/(tp+fn)
#   +     offtarg <- fp/sp.found
#   +     formula <- sp.found -offtarg
#   +     return(formula)
# OR just based on sensitivity and specificity 
# score <- function(t,p,tp,fn){
#   +     n <- t-p
#   +     tn <- n-fn
#   +     fp <- p-tp
#   +     sens <- tp/(tp+fn)
#   +     spe <- tn/(tn+fp)
#   +     return(list(sens,spe))