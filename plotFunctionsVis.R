################################################################################################
# function to plot as a heatmap the antibodies that bind the selected histone modification(s)
# Using GGVIS
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
#       - specificity
#       - target (only abs targetting those modifications?)
#
# Current ISSUES with GGVIS:
#   - binding with shiny  - doesn't work within renderUi....it shows only for some seconds
#   - map axis labels (because of names that are equal)
#   - legend position (using relative scales) -> it moves around when resizing...
#       - possibility to orient it on "top" (it's weird it's not implemented)
#   - plot title? I put it as a title for a top x axis
#   - tooltip: is it so slow? possibility to write hyperlinks? - doesn't seem to work
################################################################################################
abHeatmapV <- function(delfia, peptides, abs, 
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
  
  delf <- merge(delf, abs, by.x=c("Antibody_ID", "Barcode"), by.y=c("ID", "Barcode"))
  
  
  # apply threshold OR FILTER OUT INSIGNIFICANT COLUMNS/ROWS
  delf.filt <- delf[delf$median.signal.norm >= thres, ]
  abs.filt <- unique(delf.filt$Barcode)
  pepts.filt <- unique(delf.filt$Peptide_ID)
  
  delf <- delf[delf$Barcode %in% abs.filt &
                 delf$Peptide_ID %in% pepts.filt, ]
  
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
                          "Host", "Clonality",  "Isotype",  "Clone",	"Application")
  
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
      print(ptms)
      ptms <- paste(ptms, collapse="/")
      ptms <- gsub(" \\([A-Z][0-9]*\\)", "", ptms)
      mod <- paste(ri, ptms, sep="") 
      mods <- c(mods, mod)
    }
  }
  mods <- paste(mods, collapse=", ")
  title <- paste("Abs vs. ", prot, " ", mods,
                 " [signal > ", thres, " ]", sep="")
  
  # create plot
  abnames.labs <- unname(apply(abs[abs$Barcode %in% unique(plot.delf$ab), 
                              c("Name", "CatNr")],
                    1,
                    function(x){paste(x, collapse=" #")}))

  abnames <- factor(unique(plot.delf$ab), labels=abnames.labs)
  pepnames <- peptides[peptides$ID %in% levels(plot.delf$pep), "Name"]

  # tooltip
  show_abpep <- function(x) {
    if(is.null(x)) return(NULL)
    ab.name <- abs[abs$Barcode==x[,'factor(ab)'], 'Name']
    ab.catnr <- abs[abs$Barcode==x[,'factor(ab)'], 'CatNr']
    ab.wholename <- paste(ab.name, " (CatNr:", ab.catnr, ")", sep="")

    pep.name <- peptides[peptides$ID==x[,'factor(pep)'], 'Name']
    
    sig <- as.numeric(unique(delfia[delfia$Barcode==x[,'factor(ab)'] &
                         delfia$Peptide_ID==x[,'factor(pep)'],
                         'median.signal']))
    df <- data.frame(ab=ab.wholename,            
                     pep=pep.name ,
                     signal.norm=x$signal.norm,
                     signal=sig)
    paste0("<b>", c("Antibody", "Peptide", "Norm.signal", "Raw signal"), "</b>", 
           ": ", format(df, digits=3), collapse = "<br />")
  }
  
  # data
  p <- plot.delf %>%
        ggvis() %>%
          layer_rects(x=~factor(ab,labels=abnames.labs),y=~factor(pep), 
                   fill=~signal.norm, 
                   stroke:="#999999",
                   width=band(), height=band())
  
#   p <- p %>% layer_text(text:=title,
#                         prop("x", 0.35, "x_rel"),
#                         prop("y", 1.05, "y_rel"))
  p <- p %>%
        scale_nominal("x", padding = 0, points=FALSE) %>%
        scale_nominal("y", padding = 0, points =FALSE) %>%
        scale_numeric("fill", domain = c(thres,1), range = c("#f2f2f2","steelblue"), 
                      clamp=TRUE)
  
  p <- p %>% add_tooltip(show_abpep, "hover")

  p <- p %>%
    #add_axis("x_rel", orient="bottom", offset=-10)%>%
    add_axis("x", grid=FALSE, 
             title = "Antibodies", title_offset = 75,
             offset = 10,
             
             properties=axis_props(
               axis=list(stroke="white"),
               labels = list(angle = -45, align = "right"),
               title=list(fontSize=16))) %>%
    add_axis("y", grid=FALSE, title = "Peptides", title_offset = 50,
             offset=-10,
             properties=axis_props(
               axis=list(stroke="white"),
               title=list(fontSize=16))) %>%
    add_axis("x", grid=FALSE, orient="top", title=title, offset = 10,
             tick_size_major=0, tick_size_minor=0,
             properties=axis_props(
               axis=list(stroke="white"),
               labels=list(fontSize=0),
               title=list(fontSize=18)))
  #p <- p %>% set_options(padding=padding(150,150,200,50))
#   p <- p %>% set_options(height=40+nrow(pep.table)*20, 
#                          width=45+nrow(ab.table)*25,
#                          keep_aspect=TRUE)
#     width = 50*nrow(ab.table), height = 50*nrow(pep.table), keep_aspect = TRUE)
  p <- p %>% add_relative_scales()
  p <- p %>% add_legend("fill", title = "Normalized signal")
#                          properties=legend_props(
#                            legend=list(x=scaled_value("x_rel", 0.4), 
#                                        y=scaled_value("y_rel", 0.5))))

  # result list
  return(list(p, pep.table, ab.table, res.table))
}