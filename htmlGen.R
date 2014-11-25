createAbInfoHTML <- function(abs, ab.name, ab.catnr){
  head  <- "<html><body>"
  tail <- "</body></html>"
  
  #select antibody - use function
  ab.id <- unique(selectAntibodies(abs, ab.name, ab.catnr))
  
  # properties
  property.list <- list(c("Name", "CatNr", "LotNr", "Company", "URL", "Class", "Comment"),
                        c("Host", "Clonality", "Isotype", "Clone", "Purity",
                          "Dilution", "Applications", "Recommended_dilutions"))
  body <- ""
  ncol <- length(property.list)
  for (ab in ab.id){

    if (length(ab.id)>1){
      body <- paste(body, '<h4 style="float:center;margin:0;width:35%;">Antibody-', 
                    which(ab.id==ab),':</h4>', sep="")
    }
    for (i in 1:ncol){

      body <- paste(body, '<div id="col', i, '" style="float:left;margin:0;width:35%;">', 
                    sep="")   
      for (property in property.list[[i]]){

        prop.value <- paste(unique(as.character(abs[abs$ID==ab, property])), 
                            collapse="/")
        

        str <- paste("<br><b>", property, "</b>: ", prop.value, "<br>",
                     sep="")
        body <- paste(body, str, sep="")
      }
      body <- paste(body, '<br><br></div>', sep="")
    }
  }
  
  html <- paste(head, body, tail, sep="")

  return(body)
}

