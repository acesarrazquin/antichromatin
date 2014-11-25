categoryPiechart <- function(abs, vendor="all"){
  #remove vendor MFPL
  abs <- abs[abs$Company != "MFPL", ]
  
  uniqs <- unique(abs$ID)

  antibodies <- data.frame()

  class.names <- c("Highly specific/Highly sensitive",
                   "Highly specific/Moderately sensitive",
                   "Highly specific/Poorly sensitive",
                   "Moderately specific/Highly sensitive",
                   "Moderately specific/Moderately sensitive",
                   "Moderately specific/Poorly sensitive",
                   "Poorly specific/Highly sensitive",
                   "Poorly specific/Moderately sensitive",
                   "Poorly specific/Poorly sensitive",
                   "Not classified")
  for (uniq in uniqs){
    company <- unique(abs[abs$ID==uniq, "Company"])
    
    if (vendor != "all" && vendor != company){
      next
    }
    
    class <- abs[abs$ID == uniq, "Class"]
    
    if (length(class) > 1){
      count <- table(class)
      count <- count[names(count) != 0]
      class <- min(as.numeric(names(count[count==max(count)])))  
    }
    
    #class <- min(class[class!=0])
    if (class==Inf || class==0){
      class=Inf
      class.name <- class.names[10]
    }else{
      class.name <- class.names[class]
    }
    print(uniq)
    
    antibodies <- rbind(antibodies, data.frame(ab=uniq, company=company, class=class, class.name=class.name))
  }
  

  # calculate text positions
  counts <- unname(table(antibodies$class))
  sumcounts <- c(0, unname(cumsum(counts)))
  pos <- vector()
  for (i in 1:length(sumcounts)-1){
    pos <- c(pos, sumcounts[i]+(sumcounts[i+1]-sumcounts[i])/2)
  }
  
  perc <- paste(as.character(format(unname(counts/sum(counts))*100, digits=2, nsmall=1)),
                "%\n (", counts,")", sep="")
  
  # plot
  p <- ggplot(antibodies) + geom_bar(aes(x=factor(1), fill=factor(class)), 
                              size=1, width=1)
  p <- p + scale_y_continuous(breaks=pos, labels=perc)
  p <- p + coord_polar(theta="y")
  
  
  colors=c("Inf"="grey", "1"="#0F3B5F", "2"="#336699", "3"="#96B3D3", 
           "4"="#8E324C", "5"="#99495C", "6"="#C46279", "7"="#607848", "8"="#789048", "9"="#C0D860")
  
  
  p <- p + theme(axis.title=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text.y=element_blank(),
                 axis.text.x=element_text(color="black", vjust=-1, hjust=-1,size=20),
                 panel.background=element_rect(fill="white"),
                 plot.title=element_text(size=20, face="bold"),
#                  legend.position="none",
                 legend.key=element_blank(),
                 legend.text=element_text(size=14),
                 panel.margin=unit(c(0,20,0,0), "lines"))
  if (vendor !="all"){
    p <- p + theme(legend.position="none")
    nb.abs <- nrow(antibodies[antibodies$company==vendor,])
  }else{
    nb.abs <- nrow(antibodies)
  }
  print(nb.abs)
  p <- p + scale_fill_manual(values=colors,
                             labels=c("Inf"="Not classified", 
                                      "1"="Highly specific/Highly sensitive",
                                      "2"="Highly specific/Moderately sensitive",
                                      "3"="Highly specific/Poorly sensitive",
                                      "4"="Moderately specific/Highly sensitive",
                                      "5"="Moderately specific/Moderately sensitive",
                                      "6"="Moderately specific/Poorly sensitive",
                                      "7"="Poorly specific/Highly sensitive",
                                      "8"="Poorly specific/Moderately sensitive",
                                      "9"="Poorly specific/Poorly sensitive"),
                             name="")
  
  p <- p + labs(title=paste(vendor, " (", nb.abs, " abs)", sep=""))
  
  filename <- paste('images/abclass_', vendor, "_1109.pdf", sep="")
  pdf(filename, width=12, height=8)
  print(p)
  dev.off()
  
  return(p)
}

categoryPiechartVis <- function(abs){
  uniqs <- unique(paste(abs$Name, abs$CatNr, abs$Company, sep="$$"))
  antibodies <- data.frame()

  for (uniq in uniqs){
    split_name <- strsplit(uniq, split="\\$\\$")[[1]]
    name <- split_name[1]
    catnr <- split_name[2]
    company <- split_name[3]
    
    class <- abs[abs$Name==name & abs$CatNr==catnr & abs$Company==company,
                 "Class"]
    class <- min(class[class!=0])
#     if (class==Inf){
#       class.name <- class.names[10]
#     }else{
#       class.name <- class.names[class]
#     }
    antibodies <- rbind(antibodies, data.frame(ab=uniq, class=class))
  }
  
  counts <- table(antibodies$class)
  angles <- (counts*2*pi)/(sum(counts))

  ab.angle <- data.frame()

  end <- 0
  for (i in 1:length(angles)){
    type <- names(angles[i])
    end <- end + angles[i]
    start <- end - angles[i]
    
    ab.angle <- rbind(ab.angle, data.frame(type=type, start=start, end=end))
  }
  
  colors <- c("Inf"="grey", "1"="#0F3B5F", "2"="#336699", "3"="#96B3D3", 
                     "4"="#8E324C", "5"="#99495C", "6"="#C46279", "7"="#607848", "8"="#789048", "9"="#C0D860")
  
  classnames <- c("Not classified",
                  "Highly specific/Highly sensitive",
                   "Highly specific/Moderately sensitive",
                   "Highly specific/Poorly sensitive",
                   "Moderately specific/Highly sensitive",
                   "Moderately specific/Moderately sensitive",
                   "Moderately specific/Poorly sensitive",
                   "Poorly specific/Highly sensitive",
                   "Poorly specific/Moderately sensitive",
                   "Poorly specific/Poorly sensitive")
  p <- ab.angle %>% 
  ggvis() %>% 
  layer_arcs(x:=0, y:=0,
         innerRadius:=0, outerRadius:=100,
         startAngle:=~start, endAngle:=~end,
         fill = ~type) %>%
  scale_nominal("fill", domain=names(colors), range =unname(colors)) %>%
  scale_numeric("angle", domain=c(0, 3*pi), nice=FALSE) %>%
  #scale_numeric("x", range=c(-100,100), domain=c(-100,100)) %>%
  #scale_numeric("", range=c(0,2), domain=c(0, 2)) %>%
  add_guide_legend(fill="fill", stroke="fill", orient="right",
                   title="Ab_category", 
                   values=classnames,
                   properties=legend_props(
                     legend=list(
                       x=150,
                       y=-75))) 
  
  return(p)
}

niceHeatmap <- function(classes="all",
                        classzero=FALSE, mfpl=FALSE, 
                        unique.abs=TRUE, unique.pepts=TRUE, 
                        modinabs=TRUE, complexmod=FALSE, pep.vendor="New England Peptide Inc."){
  # select antibodies
  abs.sel <- abs
  if (!complexmod){
    abs.sel <- abs[grepl("(^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$)|/", abs$Modification), ]
  }
  if (!classzero){
    abs.sel <- abs.sel[abs.sel$Class != 0, ]
  }
  if (!mfpl){
    abs.sel <- abs.sel[abs.sel$Company != "MFPL", ]
  }
  if (classes != "all"){
    abs.sel <- abs.sel[abs.sel$Class %in% classes,]
  }
    # select only one CatNr if several
  if (unique.abs){
    abs.sel.uniq <- data.frame()
    for (id in unique(abs.sel$ID)){
      abs.id <- abs.sel[abs.sel$ID==id,]
      
      if (nrow(abs.id)>1){
        class <- abs.id$Class
        count <- table(class)
        class <- min(as.numeric(names(count[count==max(count)]))) 
        abs.id <- abs.id[abs.id$Class==class,][1,]
      }
      abs.sel.uniq <- rbind(abs.sel.uniq, abs.id)
    }
    abs.sel <- abs.sel.uniq
  }
  #abs.sel.no <- abs[!grepl("^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$", abs$Modification), ] 
  
  # select peptides
  pepts.sel <- pepts
  if (!complexmod){
  pepts.sel <- pepts[grepl("^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$", pepts$Modification), ]
  }
  if (modinabs){
    # modifications in antibodies targetting multiple
    abs.mods <- strsplit(paste(abs.sel[grepl("/", abs.sel$Modification), "Modification"], collapse="/"), split="/")[[1]]
    pepts.sel <- pepts[pepts$Modification %in% abs.sel$Modification |
                         pepts$Modification %in% abs.mods , ]
  }
  if (pep.vendor != "all"){
    pepts.sel <- pepts.sel[pepts.sel$Vendor==pep.vendor,]
  }
    # only unique modifications (even from same vendor sometimes there are several) # RANDOM choice
  if (unique.pepts){
    pepts.sel.uniq <- data.frame()
    for (name in unique(pepts.sel$Name)){
      pepts.name <- pepts.sel[pepts.sel$Name==name,]
      pepts.name <- pepts.name[pepts.name$ID %in% delfia$Peptide_ID,]
      
      if (nrow(pepts.name)>1){
        pepts.name <- pepts.name[1,]
      }
      pepts.sel.uniq <- rbind(pepts.sel.uniq, pepts.name)
    }
    pepts.sel <- pepts.sel.uniq
  }
  #pepts.sel.no <- pepts[!grepl("^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$", pepts$Modification), ]
  
  # select values
  delf <- delfia[delfia$Antibody_ID %in% abs.sel$ID, ]
  delf <- delf[delf$Peptide_ID %in% pepts.sel$ID,]
  delf <- merge(delf, pepts.sel, by.x="Peptide_ID", by.y="ID")
  delf <- merge(delf, abs.sel, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  
  delf <- delf[, c('median.signal.norm', 'median.signal',
               'Barcode', 'Name.y', 'CatNr.y', 'Replicate.y',
               'Peptide_ID', 'Name.x', 'CatNr.x')]
  delf <- unique(delf)
  
  plot.delf <- data.frame(ab=factor(delf$Barcode, levels=unique(abs.sel$Barcode)), 
                          pep=factor(as.character(delf$Peptide_ID),
                                     levels=unique(pepts.sel$ID)), 
                          signal=delf$median.signal, 
                          signal.norm=delf$median.signal.norm)
  #plot
  
  #data
  p <- ggplot(plot.delf, aes(ab, pep)) + geom_tile(aes(fill = signal.norm), colour = "grey")   #p <- p + coord_fixed(ratio = nrow(ab.table)/nrow(pep.table))
  p <- p + scale_fill_gradient(high="red", low="blue",
                               name="Normalized signal")#set_limits, same scale all the time
  # layout
  p <- p + theme(panel.grid=element_blank(), 
                 panel.background=element_blank())
  p <- p + theme(plot.title=element_text(size=16, vjust=1, face="bold"),
                 axis.title.y=element_text(size=14, face="bold", vjust=0.2), 
                 axis.text.y=element_text(size=10),
                 axis.title.x=element_text(size=14, face="bold", vjust=0.2),
                 axis.text.x=element_text(angle=90, hjust=1, size=10, vjust=0.5))
  p <- p + theme(legend.position="top", legend.direction="horizontal",
                 legend.key.size=unit(18, 'points'), legend.key=element_blank(),
                 legend.text=element_text(size=8), legend.title=element_text(size=10)) 
  p <- p + xlab("antibodies") + ylab("modifications")
  # axis ticks labels
  p <- p + scale_y_discrete(breaks=levels(plot.delf$pep),
                            labels=pepts.sel[pepts.sel$ID %in% levels(plot.delf$pep), "Name"])
  p <- p + scale_x_discrete(breaks=abs.sel$Barcode,
                            labels=unname(apply(abs.sel[abs.sel$Barcode %in% levels(plot.delf$ab), c("Name", "CatNr")],
                                                1,
                                                function(x){paste(x, collapse=" #")})))
  #print
  pdf(paste("images/whole_heatmap_allfilt_new_", paste(classes, collapse=""), ".pdf", sep=""),
      width=15, height=10)
  print(p)
  dev.off()
  return(p)
}



targetLists <- function(abs, pepts, delfia,
                        classes=c(1,2,3), complexmod.pep=FALSE, complexmod.ab=FALSE,
                        mfpl=FALSE){
  
  # select antibodies
  abs.sel <- abs[abs$Class %in% classes, ]
  if (!mfpl){
    abs.sel <- abs.sel[abs.sel$Company != "MFPL", ]
  }
  if (!complexmod.ab){
      abs.sel <- abs.sel[grepl("(^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$)|/", abs.sel$Modification),]
  }
  # select peptides
  pepts.sel <- pepts
  if (!complexmod.pep){
    pepts.sel <- pepts[grepl("^[A-Z]{1}[0-9]+[a-z]{2}[0-9]*[s|u]*$", pepts$Modification), ]
  }

  # get modifications above threshold
  prots <- unique(pepts.sel$Protein)

  mods.yes <- vector()
  abs.yes <- list()
  mods.all <- vector()
  
  for (prot in prots){
    print(paste("PROTEIN:", prot))
    
    mods <- unique(pepts.sel[pepts.sel$Protein==prot, "Modification"])
    for (mod in mods){
      mod.whole <- paste(prot, mod, sep="")
      mods.all <- c(mods.all, mod.whole)
      print(paste("Modification:", mod))
      pepts.ids <- pepts.sel[pepts.sel$Modification==mod &
                               pepts.sel$Protein==prot,
                             "ID"]
      
      abs.barcode <- unique(abs.sel[abs.sel$Protein==prot &
                                      grepl(mod, abs.sel$Modification), 
                                    "Barcode"])
      
      abs.yes.mod <- vector()
      for (ab.bar in abs.barcode){
        print(paste("ab:", as.character(abs[abs$Barcode==ab.bar, "Name"])))
        threshold <- abs[abs$Barcode==ab.bar, "Threshold"]
        print(paste("threshold:", threshold))
        values <- unique(delfia[delfia$Barcode==ab.bar & 
                                  delfia$Peptide_ID %in% pepts.ids, "median.signal.norm"])
        print(paste("delfia values:", values))
        values <- values[values >= threshold]
        print(paste("delfia values (above):", values))
        if (length(values) != 0){
          mods.yes <- c(mods.yes, mod.whole)
          abs.yes.mod <- c(abs.yes.mod, as.character(abs[abs$Barcode==ab.bar, "CatNr"]))
          print(as.character(abs[abs$Barcode==ab.bar, "CatNr"]))
        }
      }
      print(abs.yes.mod)
      if (length(abs.yes.mod) !=0){
        abs.yes[mod.whole] <- paste(unique(abs.yes.mod), collapse="|")
      }
    }
  }
  mods.yes <- unique(mods.yes)
  mods.no <- unique(setdiff(mods.all, mods.yes))
  
  #write
  write.table(data.frame(mod=mods.yes), paste("targeted", paste(classes, collapse=""),"_nocomplex.txt"), 
              quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(data.frame(mod=mods.no), paste("nontargeted", paste(classes, collapse=""),"_nocomplex.txt.txt"), 
              quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(data.frame(mod=abs.yes), paste("whotargeted", paste(classes, collapse=""),"_nocomplex.txt.txt"), 
              quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  return(list(unique(mods.yes), mods.no, abs.yes))
}