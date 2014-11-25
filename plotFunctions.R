require(ggplot2)
require(grid)
require(gridExtra)

################################################################################################
# function to plot as a heatmap the antibodies that bind the selected histone modification(s)
# input:
#   - datasets: delfia results (with Peptide_Name field)
#               peptides (with peptide.name field)
#   - other parameters:
#       - prot: only a single possibility e.g. H3, is a character vector of length 1. Can be " "
#       - res: can be a vector of residues or only one e.g. c(K4, K8...). Can be NULL.
#       - ptm: a vector of ptm c(me1, me2...). Can be NULL.
#                - if more than one residue, then they are e.g. me1 (K4), me2 (K9)...
#       - clonality: Any, Monoclonal, Polyclonal
#       - host:
#       - application:
################################################################################################
abHeatmap <- function(delfia, peptides, abs, 
                      prot, res, ptm, thres=0.25,
                      clonality="Any", host=c("Mouse", "Rat", "Rabbit"), application=NULL,
                      specificity="Low", target=FALSE){
  
  # select peptides with the modification
  pep <- selectByMod(peptides, prot, res, ptm)
  
  # merge those peptides with the delfia 
  delf <- merge(delfia, pep, by.x="Peptide_ID", by.y="ID")
  
  # apply filters to antibodies and merge with delfia+peptides
  if(target==TRUE){
    abs <- selectByMod(abs, prot, res, ptm)
  }
  
  abs <- abs[abs$Host %in% host,]
  if(!is.null(application)){
    abs <- abs[abs$Application %in% application,]
  }
  if(clonality!="Any"){
    abs <- abs[abs$Clonality==clonality,]
  }
  
  if(specificity=="High"){
    abs <- abs[abs$Class %in% 1:3,]
  }else if(specificity=="Moderate"){
    abs <- abs[abs$Class %in% 1:6,]
  }
  
  delf <- merge(delf, abs, by="Barcode")
  

  # apply threshold OR FILTER OUT INSIGNIFICANT COLUMNS/ROWS
  delf.filt <- delf[delf$median.signal.norm >= thres, ]
  
  abs.filt <- unique(delf.filt$Barcode)
  pepts.filt <- unique(delf.filt$Peptide_ID)
  
  delf <- delf[delf$Barcode %in% abs.filt &
                 delf$Peptide_ID %in% pepts.filt, ]
  print(delf)
  # if no data, return NULL
  if (nrow(delf) == 0){
    return(NULL)
  }

  # generate peptides table
  pep.table <- pep[, c('Name', 'Vendor', 'CatNr', 'LotNr', 
                       'Sequence', 'Purity', 'Aminoacids', 'Molecular_weight_g.mol' )]
  colnames(pep.table) <- c('Name', 'Vendor', 'CatNr', 'LotNr',  
                           'Sequence', 'Purity', 'Aminoacids', 'Molecular weight (g/mol)' )
  rownames(pep.table) <- 1:nrow(pep.table)
  
  # generate antibodies table
  ab.table <- abs[abs$Barcode %in% unique(delf$Barcode), 
                  c('Name', "Company", "CatNr", "LotNr", "Replicate", "Dilution", 
                    "Host", "Clonality", "Isotype",  "Clone",  "Application")]
  ab.table <- unique(ab.table)
  rownames(ab.table) <- 1:nrow(ab.table)
  colnames(ab.table) <- c("Name", "Company", "CatNr", "LotNr", "Replicate", "Dilution",
                          "Host", "Clonality",  "Isotype",	"Clone",	"Application")
  
  # generate data frame to plot
  delf <- delf[, c('median.signal.norm', 'median.signal',
                   'Barcode', 'Name.y', 'CatNr.y', 'Replicate.y',
                   'Peptide_ID', 'Name.x', 'CatNr.x')]
  delf <- unique(delf)
  
  plot.delf <- data.frame(ab=factor(delf$Barcode, levels=unique(abs$Barcode)), 
                          pep=factor(as.character(delf$Peptide_ID),
                                     levels=unique(peptides$ID)), 
                          signal=delf$median.signal, 
                          signal.norm=delf$median.signal.norm)
  
  # create results table
  res.table <- data.frame(Antibody = paste(delf$Name.y, " (", delf$CatNr.y, ") ", delf$Replicate.y, sep=""),
                           Peptide = delf$Name.x,
                           Signal = plot.delf$signal,
                           Normalized_Signal = plot.delf$signal.norm)
  rownames(res.table) <- 1:nrow(res.table)
  colnames(res.table) <- c("Antibody", "Peptide", "Med_signal", "Med_signal_norm")

  # create plot

  # data
  p <- ggplot(plot.delf, aes(ab, pep)) + geom_tile(aes(fill = signal.norm), colour = "grey") 
  #p <- p + coord_fixed(ratio = nrow(ab.table)/nrow(pep.table))
  p <- p + scale_fill_gradientn(colours=c("grey97","grey97","steelblue"), values=c(0,thres,1),
                                limits=c(0,1),
                                name="Normalized signal")#set_limits, same scale all the time
  
  # layout
  p <- p + theme(panel.grid=element_blank(), 
                 panel.background=element_blank())
  p <- p + theme(plot.title=element_text(size=16, vjust=1, face="bold"),
                 axis.title.y=element_text(size=14, face="bold", vjust=0.2), 
                 axis.text.y=element_text(size=10),
                 axis.title.x=element_text(size=14, face="bold", vjust=0.2),
                 axis.text.x=element_text(angle=45, hjust=1, size=10))
  p <- p + theme(legend.position="top", legend.direction="horizontal",
                 legend.key.size=unit(18, 'points'), legend.key=element_blank(),
                 legend.text=element_text(size=8), legend.title=element_text(size=10)) 
  
  # axis ticks labels

   print(levels(plot.delf$ab))
   print(abs[abs$Barcode %in% levels(plot.delf$ab), c("Barcode", "Name", "CatNr")])
   
  p <- p + scale_y_discrete(breaks=levels(plot.delf$pep),
                            labels=peptides[peptides$ID %in% levels(plot.delf$pep), "Name"])
  p <- p + scale_x_discrete(breaks=levels(plot.delf$ab),
                            labels=unname(apply(abs[abs$Barcode %in% levels(plot.delf$ab), c("Name", "CatNr")],
                                         1,
                                         function(x){paste(x, collapse=" #")})))

  #p <- p + theme(aspect.ratio=nrow(pep.table)/nrow(ab.table))
  #p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))

  # create title 
  mods <- vector()
  if (length(res) == 0){
    mods <- paste(ptm, collapse="/")
  }else{ 
    for(ri in res){
      if(length(ptm) != 0){
        if (length(res)>1){
          ptms <- ptm[grepl(ri, ptm)]
        }else{
          ptms <- ptm
        }
      }else{
        ptms <- ""
      }
     
      ptms <- paste(ptms, collapse="/")
      ptms <- gsub(" \\([A-Z][0-9]*\\)", "", ptms)
      mod <- paste(ri, ptms, sep="") 
      mods <- c(mods, mod)
    }
  }
  mods <- paste(mods, collapse=", ")
  
  title <- paste("Abs vs. ", prot, " ", mods,
                 " [signal > ", thres, " ]", sep="")
  
  p <- p + labs(title=title, 
                y='Peptides', x='Antibodies')
  

  return(list(p, pep.table, ab.table, res.table))
}


###########################################################################################
# Function to plot the results of the DELFIA assay PER ANTIBODY, plotting peptides in the 
# same order (as in mod.list) with repetition for multi-modifications
# 
###########################################################################################

delfiaplot <- function(ab.list, pept.list, delfia.data, mod.list, antibody.name, antibody.catnr){
  
  # Select antibody number and the corresponding delfia data
  nb <- selectAntibodies(ab.list, antibody.name, antibody.catnr)
  print(nb)
  if (antibody.name == ""){
    antibody.name <- unique(as.character(ab.list[ab.list$ID %in% nb, "Name"]))
    antibody.name <- paste(antibody.name, collapse="/")
  }
  if (antibody.catnr == ""){
    antibody.catnr <- unique(as.character(ab.list[ab.list$ID %in% nb, "CatNr"]))
  }

  
  # Select delfia data
  delfia <- data.frame()
  for (i in nb){ # in case there are several numbers
    delfia.nb <- delfia.data[delfia.data$Antibody_ID == i, ]
    delfia <- rbind(delfia, delfia.nb)
  }
  
  # Prepare data for plot -if empty, then NULL
  if (nrow(delfia) == 0){
    return(NULL)
  }
  
  # Map delfia results to the modification list
  
  delfia.new <- merge(delfia, ab.list, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  delfia.new <- merge(delfia.new, mod.list, by.x="Peptide_ID", by.y="pept_id") #look at dplyr package - could be faster (inner_join doesnt work?)
  
  
  # Order the dataset according to mod-list
  ord <- order(as.numeric(delfia.new$nr, decreasing=FALSE))
  delfia.new.ord <- delfia.new[ord,]
  
  # PLOT
  
  # general
  p <- ggplot()
  p <- p + theme_bw() 
  p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())
  p <- p + theme(plot.margin = unit(c(0.5, 1, 0.5, 2), "lines"))
  p <- p + theme(plot.title=element_text(size=20, vjust=1, face="bold"),
                 axis.title.y=element_text(size=14, vjust=0.2),
                 axis.text.y=element_text(size=12),
                 panel.grid.major.x=element_blank())
  p <- p+scale_x_discrete() + scale_y_continuous(limits=c(0,1.25), 
                                                 labels = function(x){format(x, nsmall=1)})
  
  # legend and title
  if (length(nb)>1){# && length(antibody.catnr)>1){   
    antibodyname <- antibody.name
    p <- p + theme(legend.position="top", legend.direction="horizontal",
                   legend.key=element_blank(), legend.key.size=unit(15, 'points'), 
                   legend.text=element_text(size=14), legend.title=element_text(size=14))

  }else{
    antibodyname <- paste(antibody.name, " (CatNr:", 
                          ab.list[ab.list$ID==nb[1], 'CatNr'], ")", sep="")
    p <- p + theme(legend.position="none")
  }
  p <- p + labs(title=antibodyname, y= "Normalized signal")
  
  # change grids
  cols <- length(unique(delfia.new.ord$nr)) 
  vlines <- seq(1.5,cols,1)
  
  p <- p + geom_vline(xintercept=vlines, color="#E0E0E0")
  
  # create rectangles
  uniq.ptar <- unique(delfia.new.ord[,c('nr', 'Primary_Target')])
  indexes <- which(uniq.ptar$Primary_Target=="yes")
  if(length(indexes)!=0){
    p <- p + annotate("rect", xmin=indexes-0.5, xmax=indexes+0.5, ymin=-Inf, ymax=Inf, 
                      fill="grey", alpha=0.2)##00ACC6
  }
  
  
  # plot data
  #colors <- c("#4A789C","#999999", "#7FCA9F", "#E96D63", "#F4BA70", "#9C8AA5", "#F75B68")#add more colors
  
  if (length(nb)>1){
    delfia.means <- unique(delfia.new.ord[,c("nr", "median.signal.norm", "q1.signal.norm", "q3.signal.norm", 
                                             "CatNr", "Replicate")])
    delfia.means$Cat_Nr <- paste(delfia.means$CatNr, delfia.means$Replicate, sep=" ")
    
    p <- p + geom_pointrange(data=delfia.means,
                             aes(y=median.signal.norm, x=nr, 
                                 ymax=q3.signal.norm, ymin=q1.signal.norm,
                                 colour=Cat_Nr), 
                             alpha=0.5, position="identity")
    p <- p + scale_color_brewer(palette="Set1")
    #p <- p + scale_color_manual(values=colors)
    
#     # plot threshold
#     y = ab.list[ab.list$ID %in% nb, 'Threshold']
#     print(y)
#     for (i in y){
#       p <- p + geom_hline(data=data.frame(i),
#                         aes(yintercept=i), alpha=0.5, linetype="dashed")
#     }
#     
#     # plot threshold2
#     y = ab.list[ab.list$ID %in% nb, 'Threshold2']
#     print(y)
#     for (i in y){
#       p <- p + geom_hline(data=data.frame(i),
#                           aes(yintercept=i), alpha=0.5, linetype="dashed", colour="red")
#     }
    
  }else{
    delfia.means <- unique(delfia.new.ord[,c("Antibody_ID", "Replicate", "Peptide_ID", "nr",
                                             "median.signal.norm", "q1.signal.norm", "q3.signal.norm")])
    p <- p + geom_pointrange(data=delfia.means,
                        aes(y=median.signal.norm, x=nr, 
                            ymax=q3.signal.norm, ymin=q1.signal.norm),
                        colour="#4A789C", #outlier.colour=NULL,
                        alpha=0.5, position="identity")
    
#     # plot threshold
#     y = ab.list[ab.list$ID==nb, 'Threshold']
#     print(y)
#     p <- p + geom_hline(data=data.frame(y),
#                         aes(yintercept=y), alpha=0.5, linetype="dashed")
#     
  }

  # create xlabelplot
  dataset <- unique(delfia.new.ord[, c("protein", "residue.first",
                                       "ptm.first", "mod.rest", "nr")])  
  dataset <- dataset[,c("protein", "residue.first", "ptm.first", "mod.rest")]  
  
  xplot <- xLabelplot(dataset)
  
  # merge plots
  delfia.plot <- arrangeGrob(p, xplot, ncol=1, heights=2:1)
  
  # table 
  delfia.table <- merge(delfia.new.ord, pept.list, by.x="Peptide_ID", by.y="ID")
  delfia.table <- delfia.table[,c('Name.x', 'CatNr.x', "Replicate.x", 
                                  'Name.y', 'CatNr.y',  
                                  'signal', 'signal.norm')]
  colnames(delfia.table) <- c('Ab_Name', 'Ab_CatNr', 'Ab_Replicate', 'P_Name', 
                              'P_CatNr', 'Signal', 'nSignal')
  rownames(delfia.table) <- 1:nrow(delfia.table)
  
  # list to return
  result <- list(delfia.plot, delfia.table)
  
  
  return(result)
  
}

comparePlot <- function(delfia, pepts, abs, 
                        prot1, res1, ptm1,
                        prot2, res2, ptm2,
                        fold=5,
                        clonality="Any", host=c("Rat", "Rabbit", "Mouse"), application=NULL,
                        specificity="Low"){
  
  # select peptides with the modification1
  pep1 <- selectByMod(pepts, prot1, res1, ptm1)
  
  # select peptides with the modification2 (to compare with)
  pep2 <- selectByMod(pepts, prot2, res2, ptm2)
  
  # combine both peptide sets
  pep <- rbind(pep1, pep2)
  
  # merge those peptides with the delfia results and antibodies list 
  delf <- merge(delfia, pep, by.x="Peptide_ID", by.y="ID")
  delf <- merge(delf, abs, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  
  # apply filters
  delf <- delf[delf$Host %in% host,]
  if(!is.null(application)){
    delf <- delf[delf$Application %in% application,]
  }
  if(clonality!="Any"){
    delf <- delf[delf$Clonality==clonality,]
  }
  if(specificity=="High"){
    delf <- delf[delf$Class %in% 1:3,]
  }else if(specificity=="Moderate"){
    delf <- delf[delf$Class %in% 1:6,]
  }
  
  # apply threshold for antibodies against pep1
  abs.bar <- unique(delf[delf$Peptide_ID %in% pep1$ID &
                 delf$median.signal.norm > delf$Threshold, "Barcode"])

  # for peptides in pep2 select the antibodies where the signal cocient between 
  # pep1 and pep2 is more than X fold
  results <- data.frame()
  ps <- list()
  poss <- vector()
  for (ab in abs.bar){
    print(ab)
    
    signals1 <- unique(delf[delf$Barcode==ab & delf$Peptide_ID %in% pep1$ID, 
                            c("Barcode", "CatNr.y", "Name.y", "Replicate.y",
                              "Peptide_ID", "CatNr.x", "Name.x", 
                              "median.signal.norm", "Threshold")])

    signals2 <- unique(delf[delf$Barcode==ab & delf$Peptide_ID %in% pep2$ID, 
                             c("Barcode", "CatNr.y", "Name.y", "Replicate.y",
                               "Peptide_ID", "CatNr.x", "Name.x", 
                               "median.signal.norm")])
    signals <- merge(signals1, signals2, by.x="Barcode", by.y="Barcode")
    
    signals$Fold <- signals$median.signal.norm.x/signals$median.signal.norm.y

    signals[signals$median.signal.norm.x <= signals$Threshold, 'Fold'] <- 0

   
    # porportion of positive signals/total
    pos <- sum(signals$Fold>=fold)/length(signals$Fold)
    print(pos)
    
    #signals <- signals[signals$Fold > fold, ]

    # plot
    
    if (pos != 0){
      poss <- c(poss, pos)
      results <- rbind(results, signals)

      
      plot.signals <- data.frame(Pep1=factor(signals$Peptide_ID.x), 
                              Pep2=factor(signals$Peptide_ID.y),
                              Fold=signals$Fold)

      if (unique(signals$CatNr.y.x) != ""){
        abname <- paste("Antibody: ", unique(signals$Name.y.x), " (", 
                       unique(signals$CatNr.y.x), ") ", unique(as.character(signals$Replicate.y.x)),
                       sep="")
      }else{
        abname <- paste("Antibody: ", unique(signals$Name.y.x),
                        unique(as.character(signals$Replicate.y.x)),
                        sep="")
      }
    
      p <- ggplot(plot.signals, aes(y=Pep1, x=Pep2)) + geom_tile(aes(fill = Fold), colour = "grey") 
      p <- p + scale_fill_gradient(low="grey97", high="steelblue", limits=c(fold,max(signals$Fold)),
                                   na.value = "white",
                                   name="Normalized signal fold")

      # Axis ticks labels
      p <- p + scale_y_discrete(breaks=levels(plot.signals$Pep1),
                          labels=pepts[pepts$ID %in% levels(plot.signals$Pep1), "Name"])
      p <- p + scale_x_discrete(breaks=levels(plot.signals$Pep2),
                          labels=pepts[pepts$ID %in% levels(plot.signals$Pep2), "Name"])
          
     # layout
      p <- p + theme(panel.grid=element_blank(), panel.background=element_blank())
      p <- p + theme(plot.title=element_text(size=16, vjust=1, hjust=0,face="bold"),
                   axis.title.y=element_text(size=12, face="bold", vjust=0.2), 
                   axis.text.y=element_text(size=10),
                   axis.title.x=element_text(size=12, face="bold", vjust=0.2),
                   axis.text.x=element_text(angle=45, hjust=1, size=10),
                   plot.margin = unit(c(0.5, 0, 4.5, 4), "lines"))
      p <- p + theme(legend.position="top", legend.direction="horizontal",
                    legend.key.size=unit(16, 'points'), legend.key=element_blank()) 
     
      p <- p + labs(title=abname, 
                    y=paste('Peptides w/ ', prot1, res1, ptm1, sep=""),
                    x=paste('Peptides w/ ', prot2, res2, ptm2, sep=""))
      

      
      ps <- c(ps, list(p))
     
    }
  }
  if (length(ps) == 0){
    return(NULL)
  }
  
  # results table
 
  results <- results[, c("Name.y.x", "CatNr.y.x", "Replicate.y.x",
                         "Name.x.x", "CatNr.x.x", "median.signal.norm.x",
                         "Name.x.y", "CatNr.x.y", "median.signal.norm.y",
                         "Fold")]
  colnames(results) <- c("Ab_Name", "Ab_CatNr", "Ab_Replicate",
                         "Pept1_Name", "Pept1_CatNr", "Median_nSignal1",
                         "Pept2_Name", "Pept2_CatNr", "Median_nSignal2",
                         "Fold")
  
  # order the plots according to porportion of positive signals
  ord <- order(poss, decreasing=TRUE)
  ps <- ps[ord]
  
  # arrange plots in one
  plot <- do.call(arrangeGrob, c(ps,list(ncol=1)))
  
  return(list(plot, results))

}

############################################################################
# selectByMod selects antibodies or peptides based on a given modification
# (protein, residue, ptm).
#
# So far, residues can only appear IF there is a protein specified (to be changed
# in future, probably)
#
# Special symbols (in abs list) to take into account --> patches for the moment...
# "-": "any"
# "H": "all histones"
#
# Needed:
#   - column 'Protein' e.g. H2A, H3, H4...
#   - column 'Modification' e.g. K4me1K9me2...
#############################################################################
selectByMod <- function(elements, prot, res, ptm){ 
  if (prot == " "){
    if (is.null(ptm)){ # case with nothing (show all)
      els <- elements
    }else{ # case with only ptm
      els <- data.frame()
      for (i in 1:length(ptm)){
        els <- rbind(els, elements[grepl(ptm[i], elements$Modification), ])                
      }
    }
    
  }else{ # cases with protein
    
    els.prot <- elements[elements$Protein==prot, ] 
    
    # special cases for antibodies
    if(grepl("H", prot)){
      els.prot <- rbind(els.prot, elements[elements$Protein=="H",])# case of anti-histone ab
    }
    els.prot <- rbind(els.prot, elements[elements$Protein=="-",]) # case of antibody "-"
    
    els <- data.frame()
    
    if (is.null(res)){
      if (!is.null(ptm)){ # case NO residue(s), but ptm(s)
        els <- data.frame()
        for (i in 1:length(ptm)){
          els <- rbind(els, els.prot [grepl(ptm[i], els.prot$Modification), ])
        }
      }else{ # case NO residue, NO ptm
        els <- els.prot
      }
    }else{ 
      if (is.null(ptm)){ # case resiudue(s), NO ptm
        for (i in 1:length(res)){  
          els <- rbind(els, els.prot[grepl(res[i], els.prot$Modification), ])
          
          # patch for "K-ac"-like cases
          r <- paste(substr(res[i],1,1),"-",sep="")
          els <- rbind(els, els.prot[grepl(r, els.prot$Modification), ])
        }
      }else{ # case residue(s) and ptm(s)
        if (length(res) == 1){ # ptm names are simple
          for (i in 1:length(ptm)){
            mod <- paste(res, ptm[i], sep='')
            els <- rbind(els, els.prot[grepl(mod, els.prot$Modification), ])
            
            # patch for "K-ac"-like cases
            r <- paste(substr(res,1,1),"-",sep="")
            mod <- paste(r,ptm[i], sep='')
            els <- rbind(els, els.prot[grepl(mod, els.prot$Modification), ])
          }   
        }else{ # ptm names contain residue names - problem if I have residues without modification
          for (i in 1:length(res)){
            ptm.res <- ptm[grepl(res[i], ptm)]
            
            if (length(ptm.res) == 0){
              els <- rbind(els, els.prot[grepl(res[i], els.prot$Modification), ]) 
              
              # patch for "K-ac"-like cases
              r <- paste(substr(res[i],1,1),"-",sep="")
              els <- rbind(els, els.prot[grepl(r, els.prot$Modification), ])
            }else{
              for (j in 1:length(ptm.res)){
                mod <- paste(res[i], strsplit(ptm.res[j], " ")[[1]][1], sep="")
                els <- rbind(els, els.prot[grepl(mod,els.prot$Modification), ]) 
                
                # patch for "K-ac"-like cases
                r <- paste(substr(res[i],1,1),"-",sep="")
                mod <- paste(r,strsplit(ptm.res[j], " ")[[1]][1], sep='')
                els <- rbind(els, els.prot[grepl(mod, els.prot$Modification), ])
              }
            }
          }
        }
      }
    }
  }
  return(unique(els))   
}


selectAntibodies <- function(ab.list, ab.name, ab.catnr){
  print(ab.name)
  # Select antibody number (same name can have several numbers, but catnr is unique)
  print(ab.catnr)
  if (ab.catnr != ""){
    if (ab.name == ""){
      nb <- ab.list[ab.list$CatNr == ab.catnr, "ID" ] # if several barcodes, there will be repeated elements
    }else{
      nb <- ab.list[ab.list$Name == ab.name & ab.list$CatNr == ab.catnr, "ID" ]
    }
  }else{
    if (ab.name != ""){
      nb <- ab.list[ab.list$Name == ab.name, "ID" ] # could be several
    }else{
      nb = NULL
    }
  }
  print(nb)
  return(nb)
}


xLabelplot <- function(dataset){
  data.size <- nrow(dataset)
  lab <- calculateCoords(dataset, data.size)
  
  p <- ggplot(lab) + geom_text(aes(x=x, y=y, label=label, fontface=fontface,
                                   vjust=vjust, angle=angle, hjust=hjust, size=size))
  p <- p + scale_size_identity()
  p <- p + geom_segment(aes(x=xini, xend=xfin, y=yini, yend=yini))
  p <- p + theme_bw()
  p <- p + theme(panel.grid=element_blank(), panel.border=element_blank(),
                 legend.position="none", 
                 plot.margin = unit(c(-1, 1, 0, 2), "lines"))
  #axis.text=element_blank(), axis.ticks=element_blank())#, axis.title=element_blank())
  p <- p + scale_x_continuous(labels="", breaks=NULL, limits=c(0.5, nrow(dataset)+0.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(min(lab$y), 
                                max(lab$y)+(max(nchar(as.character(dataset$mod.rest)))*data.size*0.02/110)),
                       labels = function(x){format(x, nsmall=1)}) # to align
  p <- p + labs(x="", y="")
  p <- p + theme(axis.title.y=element_text(size=14, vjust=0.2),
                 axis.text.y=element_text(size=12, color="white"),
                 axis.ticks.y=element_line(color="white"))
  #p <- p + theme(panel.margin=unit(c(0,0,0,10), "npc"))

  return(p)
  
}

calculateCoords <- function(dataset, data.size, i=1, index=0){
  labels <- data.frame()
  
  # text sizes
  base.text.size <- 5.5*110/data.size
  names <- as.character(dataset[, i])  
  uniq.names <- unique(names)
  counts <- table(names)
  
  for (el in uniq.names){
    #print(el)
    dataset.el <- dataset[dataset[, i] == el, ]
    
    span <- unname(counts[names(counts) == el])
    
    x <- index + (1 + span)/2
    y <- 0.1*i
    
    # some of these paremeters are interrelated - fine tunning
    color="#FFFFFF"
    if (i==1){
      angle = 0
      size=base.text.size
      fontface="bold"
      hjust=0.5
      vjust=0
      offset=0.07
    }else if(i==2){
      angle = 0
      size=base.text.size/1.2
      fontface="bold"
      hjust=0.5
      vjust=0
      offset=0.07
    }else if(i==3){
      angle=90
      size=base.text.size/1.4
      fontface="bold"
      hjust=0.2
      vjust=0.3
      offset=0.09
    }else if(i==4){
      angle=90
      size=base.text.size/1.6
      fontface="plain"
      hjust=0
      vjust=0.25
      offset=0.09
      
    }
    
    xini <- 0
    xfin <- 0
    yini <- 0
    if (i < ncol(dataset)){
      if (el != "" & (span > 1 | dataset.el[, i+1][1] != "")){
        xini <- index + 0.75
        xfin <- index + span + 0.25
        yini <- y + offset
        
      }
    }
    
    # correction for non-modified peptides
    if (i==4 & grepl("\\(.*\\)", el)){
      y <- 0.2
    }
    
    labels.el <- data.frame(type=i, label=el, x=x, y=y, 
                            angle=angle, size=size, fontface=fontface, 
                            hjust=hjust, vjust=vjust,
                            xini=xini, xfin=xfin, yini=yini)
    labels <- rbind(labels, labels.el)
    
    j <- i+1
    if (j <= ncol(dataset)){
      labels.n <- calculateCoords(dataset.el, data.size, i=j, index=index)
      labels <- rbind(labels, labels.n)
    }
    
    index <- index + span
  }

  return(labels)
}

delfiaplot2 <- function(ab.list, pept.list, delfia.data, mod.list, antibody.name, antibody.catnr){
  
  # Select antibody number and the corresponding delfia data
  nb <- selectAntibodies(ab.list, antibody.name, antibody.catnr)
  print(nb)
  if (antibody.name == ""){
    antibody.name <- unique(as.character(ab.list[ab.list$ID %in% nb, "Name"]))
    antibody.name <- paste(antibody.name, collapse="/")
  }
  if (antibody.catnr == ""){
    antibody.catnr <- unique(as.character(ab.list[ab.list$ID %in% nb, "CatNr"]))
  }
  
  
  # Select delfia data
  delfia <- data.frame()
  for (i in nb){ # in case there are several numbers
    delfia.nb <- delfia.data[delfia.data$Antibody_ID == i, ]
    delfia <- rbind(delfia, delfia.nb)
  }
  
  # Prepare data for plot -if empty, then NULL
  if (nrow(delfia) == 0){
    return(NULL)
  }
  
  # Map delfia results to the modification list
  
  delfia.new <- merge(delfia, ab.list, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  delfia.new <- merge(delfia.new, mod.list, by.x="Peptide_ID", by.y="pept_id") #look at dplyr package - could be faster (inner_join doesnt work?)
  
  
  # Order the dataset according to mod-list
  ord <- order(as.numeric(delfia.new$nr, decreasing=FALSE))
  delfia.new.ord <- delfia.new[ord,]
  
  # PLOT
  
  # plot data
  #colors <- c("#4A789C","#999999", "#7FCA9F", "#E96D63", "#F4BA70", "#9C8AA5", "#F75B68")#add more colors
  
  if (length(nb)>1){
    delfia.means <- unique(delfia.new.ord[,c("nr", "median.signal.norm", "q1.signal.norm", "q3.signal.norm", 
                                             "CatNr", "Replicate")])
    delfia.means$Cat_Nr <- paste(delfia.means$CatNr, delfia.means$Replicate, sep=" ")
    
    p <- ggplot(data=delfia.means)
    p <- p + geom_errorbarh(aes(x=median.signal.norm, y=nr, 
                                 xmax=q3.signal.norm, xmin=q1.signal.norm,
                                 colour=Cat_Nr), 
                             alpha=0.5, position="identity")
    p <- p + geom_point(aes(x=median.signal.norm, y=nr, 
                            color=Cat_Nr),
                        alpha=0.5)
    p <- p + scale_color_brewer(palette="Set1")
    
  }else{
    delfia.means <- unique(delfia.new.ord[,c("Antibody_ID", "Replicate", "Peptide_ID", "nr",
                                             "median.signal.norm", "q1.signal.norm", "q3.signal.norm")])
    
    p <- ggplot(data=delfia.means)
    p <- p + geom_errorbarh(aes(x=median.signal.norm, y=nr, 
                                xmax=q3.signal.norm, xmin=q1.signal.norm), 
                            color="#4A789C",
                            alpha=0.5, position="identity")
    p <- p + geom_point(aes(x=median.signal.norm, y=nr), 
                        color="#4A789C",
                        alpha=0.5)   
  }
  
  # layout
  p <- p + theme_bw() 
  p <- p + theme(axis.ticks.y=element_blank(), axis.title.y=element_blank(),
                 axis.text.y=element_blank())
  p <- p + theme(plot.margin = unit(c(0.5, 2, 0.5, -1), "lines"))
  p <- p + theme(plot.title=element_blank(),#text(size=20, vjust=1, face="bold"),
                 axis.title.x=element_text(size=14, vjust=0.2),
                 axis.text.x=element_text(size=12),
                 panel.grid.major.x=element_blank())
  maxval <- min(1.2, max(delfia.means$q3.signal.norm))
  p <- p+scale_y_discrete() + scale_x_continuous(limits=c(0,maxval),
                                                 breaks=c(0, 0.25, 0.5, 0.75, 1),
                                               labels = function(x){format(x, nsmall=1)})

  # legend and title
  if (length(nb)>1){# && length(antibody.catnr)>1){   
    antibodyname <- antibody.name
    p <- p + theme(legend.position="right", legend.direction="vertical",
                   legend.key=element_blank(), legend.key.size=unit(15, 'points'), 
                   legend.text=element_text(size=14), legend.title=element_text(size=14))
    
  }else{
    antibodyname <- paste(antibody.name, " (CatNr:", 
                          ab.list[ab.list$ID==nb[1], 'CatNr'], ")", sep="")
    p <- p + theme(legend.position="none")
  }
  p <- p + labs(title=antibodyname, x= "Normalized signal")


  # create rectangles
  uniq.ptar <- unique(delfia.new.ord[,c('nr', 'Primary_Target')])
  indexes <- which(uniq.ptar$Primary_Target=="yes")
  if(length(indexes)!=0){
    p <- p + annotate("rect", ymin=indexes-0.5, ymax=indexes+0.5, xmin=-Inf, xmax=Inf, 
                      fill="grey", alpha=0.2)##00ACC6
  }
  
  # create xlabelplot
  dataset <- unique(delfia.new.ord[, c("protein", "residue.first",
                                       "ptm.first", "mod.rest", "nr")])  
  dataset <- dataset[,c("protein", "residue.first", "ptm.first", "mod.rest")]  
  
  xplot <- xLabelplot2(dataset)
  
  # merge plots
  delfia.plot <- arrangeGrob(xplot, p, ncol=2, widths=1:2)
  
  # table 
#   delfia.table <- merge(delfia.new.ord, pept.list, by.x="Peptide_ID", by.y="ID")
#   delfia.table <- delfia.table[,c('Name.x', 'CatNr.x', "Replicate.x", 
#                                   'Name.y', 'CatNr.y',  
#                                   'signal', 'signal.norm')]
#   colnames(delfia.table) <- c('Ab_Name', 'Ab_CatNr', 'Ab_Replicate', 'P_Name', 
#                               'P_CatNr', 'Signal', 'nSignal')
#   rownames(delfia.table) <- 1:nrow(delfia.table)
#   
#   # list to return
#   delfia.plot=p
#   result <- list(delfia.plot, delfia.table)
  

  return(delfia.plot)
  
}

xLabelplot2 <- function(dataset){
  data.size <- nrow(dataset)
  lab <- calculateCoords2(dataset, data.size)

  p <- ggplot(lab) + geom_text(aes(x=x, y=y, label=label, fontface=fontface,
                                   vjust=vjust, angle=angle, hjust=hjust, size=size))
  p <- p + scale_size_identity()
  p <- p + geom_segment(aes(y=yini, yend=yfin, x=xini, xend=xini))
  p <- p + theme_bw()
  p <- p + theme(panel.grid=element_blank(), #panel.border=element_blank(),
                 legend.position="none", 
                 plot.margin = unit(c(1, 0, 0.5, 2), "lines"))
  p <- p + theme(axis.title.x=element_text(size=14, vjust=0.2),
                 axis.text.x=element_text(size=12, color="black"),
                 axis.ticks.x=element_line(color="black"))
  
  #axis.text=element_blank(), axis.ticks=element_blank())#, axis.title=element_blank())
  p <- p + scale_y_continuous(labels="", breaks=NULL, 
                              limits=c(0.5, nrow(dataset)+0.5), expand=c(0,0)) +
           scale_x_continuous(limits=c(min(lab$x), 
                                max(lab$xmax)),
                                labels = function(x){format(x, nsmall=1)}) # to align
  p <- p + labs(x="Normalized signal", y="")

  return(p)
}


calculateCoords2 <- function(dataset, data.size, i=1, index=0, 
                             hztls=NULL,
                             font.size.max=5, font.sizes=NULL){
  labels <- data.frame()
  
  # get horizontal positions, if not yet calculated
  if (i==1){
    hztls <- c(0)
    font.sizes <- c()
    for (j in 1:ncol(dataset)){
      print(paste("Level:", j))
      # corrected font sizes
      font.size <- font.size.max - font.size.max*0.1*(j-1)
      print(font.size)
      font.sizes <- c(font.sizes, font.size)
      
      # horizontal positions
      nchar.max <- max(nchar(as.character(dataset[, j])))
      print(nchar.max)
      hztl <- hztls[j] + (nchar.max+0.5)*font.size
      print(hztl)
      hztls <- c(hztls, hztl)
    }
    print(hztls)
    print(font.sizes)
  }

  # set common parameters:
  color <- "#FFFFFF"
  angle <- 0
  hjust <- 0
  vjust <- 0.5
  if (i %in% c(1,2)){
    fontface <- "bold"
  }else{
    fontface <- "plain"
  }
  offset <- hztls[i+1] - 0.5*font.sizes[i]
  x <- hztls[i] 
  xmax <- hztls[length(hztls)] 
  size <- font.sizes[i]
  
  # iterate
  names <- as.character(dataset[, i])  
  uniq.names <- unique(names)
  counts <- table(names)
  for (el in uniq.names){
    #print(el)
    dataset.el <- dataset[dataset[, i] == el, ]
  
    span <- unname(counts[names(counts) == el])

    y <- index + (1 + span)/2   
    yini <- 0
    yfin <- 0
    xini <- offset
    if (i < ncol(dataset)){
      if (el != "" & (span > 1 | dataset.el[, i+1][1] != "")){
        yini <- index + 0.75
        yfin <- index + span + 0.25
        xini <- offset  
      }
    }
    
    # correction for non-modified peptides
    if (i==4 & grepl("\\(.*\\)", el)){
      x <- hztls[2] 
    }
    
    labels.el <- data.frame(type=i, label=el, x=x, y=y, 
                            angle=angle, size=size, fontface=fontface, 
                            hjust=hjust, vjust=vjust,
                            yini=yini, yfin=yfin, xini=xini, xmax=xmax)
    
    labels <- rbind(labels, labels.el)
    
    j <- i+1
    if (j <= ncol(dataset)){
      labels.n <- calculateCoords2(dataset.el, data.size, i=j, index=index, 
                                   hztls=hztls, font.sizes=font.sizes)
      labels <- rbind(labels, labels.n)
    }   
    index <- index + span
  }  
  return(labels)
}

comparePlot <- function(delfia, pepts, abs, 
                        prot1, res1, ptm1,
                        prot2, res2, ptm2,
                        fold=5,
                        clonality="Any", host=c("Rat", "Rabbit", "Mouse"), application=NULL,
                        specificity="Low"){
  
  # select peptides with the modification1
  pep1 <- selectByMod(pepts, prot1, res1, ptm1)
  
  # select peptides with the modification2 (to compare with)
  pep2 <- selectByMod(pepts, prot2, res2, ptm2)
  
  # combine both peptide sets
  pep <- rbind(pep1, pep2)
  
  # merge those peptides with the delfia results and antibodies list 
  delf <- merge(delfia, pep, by.x="Peptide_ID", by.y="ID")
  delf <- merge(delf, abs, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  
  # apply filters
  delf <- delf[delf$Host %in% host,]
  if(!is.null(application)){
    delf <- delf[delf$Application %in% application,]
  }
  if(clonality!="Any"){
    delf <- delf[delf$Clonality==clonality,]
  }
  if(specificity=="High"){
    delf <- delf[delf$Class %in% 1:3,]
  }else if(specificity=="Moderate"){
    delf <- delf[delf$Class %in% 1:6,]
  }
  
  # apply threshold for antibodies against pep1
  abs.bar <- unique(delf[delf$Peptide_ID %in% pep1$ID &
                           delf$median.signal.norm > delf$Threshold, "Barcode"])
  
  # for peptides in pep2 select the antibodies where the signal cocient between 
  # pep1 and pep2 is more than X fold
  results <- data.frame()
  ps <- list()
  poss <- vector()
  for (ab in abs.bar){
    print(ab)
    
    signals1 <- unique(delf[delf$Barcode==ab & delf$Peptide_ID %in% pep1$ID, 
                            c("Barcode", "CatNr.y", "Name.y", "Replicate.y",
                              "Peptide_ID", "CatNr.x", "Name.x", 
                              "median.signal.norm", "Threshold")])
    
    signals2 <- unique(delf[delf$Barcode==ab & delf$Peptide_ID %in% pep2$ID, 
                            c("Barcode", "CatNr.y", "Name.y", "Replicate.y",
                              "Peptide_ID", "CatNr.x", "Name.x", 
                              "median.signal.norm")])
    signals <- merge(signals1, signals2, by.x="Barcode", by.y="Barcode")
    
    signals$Fold <- signals$median.signal.norm.x/signals$median.signal.norm.y
    
    signals[signals$median.signal.norm.x <= signals$Threshold, 'Fold'] <- 0
    
    
    # porportion of positive signals/total
    pos <- sum(signals$Fold>=fold)/length(signals$Fold)
    print(pos)
    
    #signals <- signals[signals$Fold > fold, ]
    
    # plot
    
    if (pos != 0){
      poss <- c(poss, pos)
      results <- rbind(results, signals)
      
      
      plot.signals <- data.frame(Pep1=factor(signals$Peptide_ID.x), 
                                 Pep2=factor(signals$Peptide_ID.y),
                                 Fold=signals$Fold)
      
      if (unique(signals$CatNr.y.x) != ""){
        abname <- paste("Antibody: ", unique(signals$Name.y.x), " (", 
                        unique(signals$CatNr.y.x), ") ", unique(as.character(signals$Replicate.y.x)),
                        sep="")
      }else{
        abname <- paste("Antibody: ", unique(signals$Name.y.x),
                        unique(as.character(signals$Replicate.y.x)),
                        sep="")
      }
      
      p <- ggplot(plot.signals, aes(y=Pep1, x=Pep2)) + geom_tile(aes(fill = Fold), colour = "grey") 
      p <- p + scale_fill_gradient(low="grey97", high="steelblue", limits=c(fold,max(signals$Fold)),
                                   na.value = "white",
                                   name="Normalized signal fold")
      
      # Axis ticks labels
      p <- p + scale_y_discrete(breaks=levels(plot.signals$Pep1),
                                labels=pepts[pepts$ID %in% levels(plot.signals$Pep1), "Name"])
      p <- p + scale_x_discrete(breaks=levels(plot.signals$Pep2),
                                labels=pepts[pepts$ID %in% levels(plot.signals$Pep2), "Name"])
      
      # layout
      p <- p + theme(panel.grid=element_blank(), panel.background=element_blank())
      p <- p + theme(plot.title=element_text(size=16, vjust=1, hjust=0,face="bold"),
                     axis.title.y=element_text(size=12, face="bold", vjust=0.2), 
                     axis.text.y=element_text(size=10),
                     axis.title.x=element_text(size=12, face="bold", vjust=0.2),
                     axis.text.x=element_text(angle=45, hjust=1, size=10),
                     plot.margin = unit(c(0.5, 0, 4.5, 4), "lines"))
      p <- p + theme(legend.position="top", legend.direction="horizontal",
                     legend.key.size=unit(16, 'points'), legend.key=element_blank()) 
      
      p <- p + labs(title=abname, 
                    y=paste('Peptides w/ ', prot1, res1, ptm1, sep=""),
                    x=paste('Peptides w/ ', prot2, res2, ptm2, sep=""))
      
      
      
      ps <- c(ps, list(p))
      
    }
  }
  if (length(ps) == 0){
    return(NULL)
  }
  
  # results table
  
  results <- results[, c("Name.y.x", "CatNr.y.x", "Replicate.y.x",
                         "Name.x.x", "CatNr.x.x", "median.signal.norm.x",
                         "Name.x.y", "CatNr.x.y", "median.signal.norm.y",
                         "Fold")]
  colnames(results) <- c("Ab_Name", "Ab_CatNr", "Ab_Replicate",
                         "Pept1_Name", "Pept1_CatNr", "Median_nSignal1",
                         "Pept2_Name", "Pept2_CatNr", "Median_nSignal2",
                         "Fold")
  
  # order the plots according to porportion of positive signals
  ord <- order(poss, decreasing=TRUE)
  ps <- ps[ord]
  
  # arrange plots in one
  plot <- do.call(arrangeGrob, c(ps,list(ncol=1)))
  
  return(list(plot, results))
  
}
#OLD STUFF

# exampleXtext <- function(){
#   ex <- data.frame(
#     name=c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6),
#     type=c(0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3),
#     value=c("", "", "","H3","","","","K4","","","K9","","ac","me1", "me2", "ac", "me1", "me3","","","S10ph","","",""),
#     size = c(14,14,14,14,14,14,13,13,13,13,13,13,12,12,12,12,12,12,11,11,11,11,11,11),
#     angle= c(0,0,0,0,0,0,0,0,0,0,0,0,90,90,90,90,90,90,90,90,90,90,90,90),
#     count=c(0,0,0,6,0,0,0,3,0,0,3,0,0,0,1,0,0,0,0,0,0,0,0,0))
# 
#   p <- ggplot() + geom_text(data=ex, aes(x=name, y=type, label=value, size=size, angle=angle))
# 
#   width <- nrow(ex)/length(unique(ex$type))
#   segments <- unique(ex[!ex$type==3 & !ex$value=="" &!ex$count==0, c('type', 'value', 'count')])
#   
#   xini=1
#   print(segments)
#   segments.draw <- data.frame()
#   for (i in 1:nrow(segments)){
#  
#     print(xini)
#     segment <- segments[i, ]
#     seg.type <- segment$type   
# 
#     xfin <- segment$count + xini - 1
#     print(xfin)
#     if(seg.type==0){
#       yini=0.8
#     }else if (seg.type == 1){
#       yini=1.8
#     }else if (seg.type == 2){
#       yini=2.8
#     }
#     #print(yini)
#     #p <- p + geom_segment(aes(x=xini , xend=xfin , y=yini, yend=yini, size=1))
#     # it doesnt work like this...
#     segments.draw <- rbind(segments.draw, data.frame(x=xini-0.25 , xend=xfin+0.25 , y=yini, yend=yini))
#     xini <- min(xfin+1, abs(width-xfin+1))
#     
#   }
#   
#   p <- p + geom_segment(aes(x=x, xend=xend, y=y, yend=yend), data=segments.draw)
#   p <- p + theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank(), legend.position="none",
#                               axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
#   print(p)
# }