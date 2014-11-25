library(shiny)
require(ggplot2)
require(grid)
require(gridExtra)

## PREPARE:

# load functions
source("plotFunctions.R")
source("dataFunctions.R")
source("htmlGen.R")


# read files

# delfia <- read.csv('data/delfia_0909.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
# pepts <- read.csv('data/peptides.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
# abs <- read.csv('data/antibodies_barcodes_20140922.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

# # modify data
# mod.list <- createModList(pepts)

# delfia <- annotatePrimaryTargets(delfia, pepts, abs) 
# delfia <- normalizeDelfia(delfia)


# abs <- calculateThresholds(delfia, abs)
# abs <- classifyAntibodies(delfia, abs)
# abs <- annotateRepeatedAb(abs)
# abs <- orderElms(abs, type="abs")

# pepts <- orderElms(pepts, type="pepts")

## save.image("antiChromatin.RData")

# # TO LOAD
load("antiChromatin.RData")

## SERVER:

# Define server logic required 
shinyServer(function(input, output, session) {
  
  
  # INPUTS
  
  proteins <- reactive({
    x <- sort(unique(as.character(mod.list$protein)))
    x <- c(" ", x[x!='p53'])  
  }) 
  
  observe({
    updateSelectInput(session, "prot", choices=proteins(), selected=" ")
  })
  
  
  residues <- reactive({
    x <- unique(as.character(mod.list[mod.list$protein==input$prot, ]$residue.first))
    x <- c("", x) # otherwise when there is no protein still shows residues...why?
    
  })
  observe({
    updateSelectInput(session, "res", choices=residues())    
  })
  
  ptm <- reactive({
    
    if (input$prot == " "){
      x <- sort(unique(as.character(mod.list$ptm.first)))
      
    }else if (input$prot != " " & is.null(input$res)){
      x <- sort(unique(as.character(mod.list[mod.list$protein==input$prot, ]$ptm.first)))
      
    }else if (input$prot != " " & !is.null(input$res)){
      
      if (length(input$res) == 1){
        x <- sort(unique(as.character(mod.list[mod.list$protein==input$prot &
                                                 mod.list$residue.first==input$res, ]$ptm.first)))
      }else{
        x <- unique(mod.list[mod.list$protein==input$prot &
                               mod.list$residue.first %in% input$res, c('ptm.first', 'residue.first')])
        
        x.res <- sort(as.character(paste(x$ptm.first, ' (', x$residue.first, ')', sep = '')))
        #x.table <- sort(table(x$ptm.first))
        #x.intersect <- names(x.table[x.table==length(input$res)])
        #x <- c(x.intersect, x.res)
        x <- x.res
      } 
    }
    #x <- c(" ", x[x!=''])
  })
  
  ptm.selected <- reactive({
    if (!is.null(input$ptm) & !is.null(input$res)){
      if (!grepl('[A-Z]', input$ptm[1])){
        ptm.sel <- input$ptm
        
        if (length(input$res) > 1){
          ptm.sel <- unname(sapply(ptm.sel, function(x){paste(x," (",input$res[1],")",sep="")}))    
        }
      }else{
        ptm.all <- input$ptm
        ptm.sel <- vector()
        for (res in input$res){
          ptm.sel <- c(ptm.sel, ptm.all[grepl(res, ptm.all)])
        }
        
        if (length(input$res) == 1){
          ptm.sel <- unname(sapply(ptm.sel, function(x) strsplit(x, split=' ')[[1]][1]))
        }
      }
      ptm.sel  
      
    }else{
      input$ptm 
    }
  })
  
  observe({
    updateSelectInput(session, "ptm", choices=ptm(), selected=ptm.selected())
  })
  
  species <- reactive({
    x <- sort(unique(as.character(abs$Host))) 
    x <- x[!x==""]
  })
  
   observe({
       updateCheckboxGroupInput(session, "species", choices=species(), selected=species())
   })
  
  application <- reactive({
    app <- sort(unique(as.character((abs$Application))))
    app <- sapply(app, function(x) strsplit(x, ",")[[1]])
    app <- unique(unname(unlist(app)))
    
  })
  
  observe({
    updateSelectInput(session, "application", choices=application(), selected="")
  })
  
  choose.abs <- reactive({
    antibodies <- abheatmap()[[3]]
    abs.name.catnr <- paste(antibodies$Name, " - ", antibodies$CatNr, sep="")
    abs.name.catnr <- c(" ", abs.name.catnr)
    })
  
  observe({
    updateSelectInput(session, "chooseab", choices=choose.abs(), selected="")
  })
  
  observe({
    if(input$submit_abchoice){     
      updateTabsetPanel(session, "navbar", selected="by Antibody")
    }
  })
  
  
  ## BY ANTIBODY
  ab.name <- reactive({
    x <- unique(as.character(abs$Name))
#     x <- c(" ", x)
  })
  
  observe({
    if(input$submit_mod){
      name <- strsplit(input$chooseab, split=" - ")[[1]][1]
#     }else if(input$submit_comp){
#       name <- strsplit(input$chooseabX, split=" - ")[[1]][1]
    }else{
      name=""
    }
    updateSelectInput(session, "ab.name", choices=ab.name(), selected=name)
  })
  
  ab.catnr <- reactive({
    if (input$ab.name==""){
      x <- unique(as.character(abs$CatNr))
    }else{
      x <- unique(as.character(abs[abs$Name==input$ab.name, ]$CatNr))
    }
#     x <- c(" ", x)
  })

  observe({
    if(input$submit_mod){
      catnr <- strsplit(input$chooseab, split=" - ")[[1]][2]
#     }else if(input$submit_comp){
#       catnr <- strsplit(input$chooseabX, split=" - ")[[1]][2]
    }else{
      catnr=""
    }
    updateSelectInput(session, "ab.catnr", choices=ab.catnr(), selected=catnr)
  })
  
  # COMPARE 
  
#   ##1 (I should look at ways of reducing all these repetitions)
#   observe({ # proteins function in input by modification
#     updateSelectInput(session, "protX1", choices=proteins(), selected=" ")
#   })
#   
#   
#   residuesX1 <- reactive({
#     x <- unique(as.character(mod.list[mod.list$protein==input$protX1, ]$residue.first))
#     x <- c("", x) # otherwise when there is no protein still shows residues...why?
#     
#   })
#   observe({
#     updateSelectInput(session, "resX1", choices=residuesX1())    
#   })
#   
#   ptmX1 <- reactive({
#     
#     if (input$protX1 == " "){
#       x <- sort(unique(as.character(mod.list$ptm.first)))
#       
#     }else if (input$protX1 != " " & is.null(input$resX1)){
#       x <- sort(unique(as.character(mod.list[mod.list$protein==input$protX1, ]$ptm.first)))
#       
#     }else if (input$protX1 != " " & !is.null(input$resX1)){
#       
#       if (length(input$resX1) == 1){
#         x <- sort(unique(as.character(mod.list[mod.list$protein==input$protX1 &
#                                                  mod.list$residue.first==input$resX1, ]$ptm.first)))
#       }else{
#         x <- unique(mod.list[mod.list$protein==input$protX1 &
#                                mod.list$residue.first %in% input$resX1, c('ptm.first', 'residue.first')])
#         
#         x.res <- sort(as.character(paste(x$ptm.first, ' (', x$residue.first, ')', sep = '')))
#         #x.table <- sort(table(x$ptm.first))
#         #x.intersect <- names(x.table[x.table==length(input$res)])
#         #x <- c(x.intersect, x.res)
#         x <- x.res
#       } 
#     }
#     #x <- c(" ", x[x!=''])
#   })
#   
#   ptm.selectedX1 <- reactive({
#     if (!is.null(input$ptmX1) & !is.null(input$resX1)){
#       if (!grepl('[A-Z]', input$ptmX1[1])){
#         ptm.sel <- input$ptmX1
#         
#         if (length(input$resX1) > 1){
#           ptm.sel <- unname(sapply(ptm.sel, function(x){paste(x," (",input$resX1[1],")",sep="")}))    
#         }
#       }else{
#         ptm.all <- input$ptmX1
#         ptm.sel <- vector()
#         for (res in input$resX1){
#           ptm.sel <- c(ptm.sel, ptm.all[grepl(res, ptm.all)])
#         }
#         
#         if (length(input$resX1) == 1){
#           ptm.sel <- unname(sapply(ptm.sel, function(x) strsplit(x, split=' ')[[1]][1]))
#         }
#       }
#       ptm.sel  
#       
#     }else{
#       input$ptmX1 
#     }
#   })
#   
#   observe({
#     updateSelectInput(session, "ptmX1", choices=ptmX1(), selected=ptm.selectedX1())
#   })
#   
#   ##2
#   observe({ # proteins function in input by modification
#     updateSelectInput(session, "protX2", choices=proteins(), selected=" ")
#   })
#   
#   
#   residuesX2 <- reactive({
#     x <- unique(as.character(mod.list[mod.list$protein==input$protX2, ]$residue.first))
#     x <- c("", x) # otherwise when there is no protein still shows residues...why?
#     
#   })
#   observe({
#     updateSelectInput(session, "resX2", choices=residuesX2())    
#   })
#   
#   ptmX2 <- reactive({
#     
#     if (input$protX2 == " "){
#       x <- sort(unique(as.character(mod.list$ptm.first)))
#       
#     }else if (input$protX2 != " " & is.null(input$resX2)){
#       x <- sort(unique(as.character(mod.list[mod.list$protein==input$protX2, ]$ptm.first)))
#       
#     }else if (input$protX2 != " " & !is.null(input$resX2)){
#       
#       if (length(input$resX2) == 1){
#         x <- sort(unique(as.character(mod.list[mod.list$protein==input$protX2 &
#                                                  mod.list$residue.first==input$resX2, ]$ptm.first)))
#       }else{
#         x <- unique(mod.list[mod.list$protein==input$protX2 &
#                                mod.list$residue.first %in% input$resX2, c('ptm.first', 'residue.first')])
#         
#         x.res <- sort(as.character(paste(x$ptm.first, ' (', x$residue.first, ')', sep = '')))
#         #x.table <- sort(table(x$ptm.first))
#         #x.intersect <- names(x.table[x.table==length(input$res)])
#         #x <- c(x.intersect, x.res)
#         x <- x.res
#       } 
#     }
#     #x <- c(" ", x[x!=''])
#   })
#   
#   ptm.selectedX2 <- reactive({
#     if (!is.null(input$ptmX2) & !is.null(input$resX2)){
#       if (!grepl('[A-Z]', input$ptmX2[1])){
#         ptm.sel <- input$ptmX2
#         
#         if (length(input$resX2) > 1){
#           ptm.sel <- unname(sapply(ptm.sel, function(x){paste(x," (",input$resX2[1],")",sep="")}))    
#         }
#       }else{
#         ptm.all <- input$ptmX2
#         ptm.sel <- vector()
#         for (res in input$resX2){
#           ptm.sel <- c(ptm.sel, ptm.all[grepl(res, ptm.all)])
#         }
#         
#         if (length(input$resX2) == 1){
#           ptm.sel <- unname(sapply(ptm.sel, function(x) strsplit(x, split=' ')[[1]][1]))
#         }
#       }
#       ptm.sel  
#       
#     }else{
#       input$ptmX2 
#     }
#   })
#   
#   observe({
#     updateSelectInput(session, "ptmX2", choices=ptmX2(), selected=ptm.selectedX2())
#   })
#   
#   ##general
#   
#   observe({
#     updateRadioButtons(session, "speciesX", choices=species(), selected=species())
#   })
#   
#   observe({
#     updateCheckboxGroupInput(session, "applicationX", choices=application(), selected="")
#   })
#   
#   choose.absX <- reactive({
#     if(input$submit_comp){
#       antibodies <- compare()[[2]]
#       abs.name.catnr <- sort(unique(paste(antibodies$Ab_Name, " - ", 
#                                         antibodies$Ab_CatNr, sep="")))
#       abs.name.catnr <- c(" ", abs.name.catnr)
#     }
#   })
#   
#   observe({
#     updateSelectInput(session, "chooseabX", choices=choose.absX(), selected="")
#   })
#   
#   observe({
#     if(input$submit_abchoiceX){     
#       updateTabsetPanel(session, "navbar", selected="by Antibody")
#     }
#   })
#   
#   
  
  #   ptm.error <- reactive({
  #     alert <- 'FALSE'
  #     if (!is.null(input$ptm) & sum(grepl("^[a-z]{2}[0-9]*$", input$ptm)) == 0){      
  #       for (res.input in input$res){
  #         if (sum(grepl(res.input, input$ptm)) == 0){
  #           alert <- TRUE
  #           break
  #         }
  #       }
  #     }
  #     alert
  #   })
  # 
  #   output$ptm.error <- renderText({ 
  #     if (input$submit_mod == 0) return()
  #     
  #     isolate(
  #       if (ptm.error() == TRUE){
  #         'You must select at least one PTM per residue'
  #       }else{""}
  #     )
  #   })
  
  # OUTPUTS

  ## ENTER
  output$cdlogo <- renderImage({
    return(list(src = "CDlab_white.png",
        filetype = "image/png",
        alt = "CDlabs"    
  ))}, deleteFile = FALSE
  )
  

  
  ## HOME
#   output$abclassification <- renderPlot({
#     print(categoryPiechart(abs))
#   })

  output$abclassification <- renderImage({
    list(src="./www/abclassification.jpg")
  }, deleteFile=FALSE)
  
  output$cemmlogo <- renderImage({
    #filename <- normalizePath(file.path('./www','CeMMlogo.jpg'))
    list(src="./www/CeMMlogo.jpg")},
    deleteFile=FALSE)
  
  ## SEARCH: by modification
  
  abheatmap <- reactive({
    if(!input$submit_mod){return()}

    isolate(
      abHeatmap(delfia, pepts, abs, 
                prot=input$prot, res=input$res, ptm=input$ptm, thres=input$threshold,
                 clonality=input$clonality, host=input$species, application=input$application,
                 specificity=input$specificity,target=input$target)
    )
  
  })
  
  output$abheatmap <- renderPlot({ 
    if(!input$submit_mod){return()}
  
    isolate({

      print(abheatmap()[[1]])
    })
  })
  
  output$pept_table <- renderDataTable({
    input$submit_mod
    
    isolate({
      abheatmap()[[2]]
    })
  })
  
  output$res_table <- renderDataTable({
    input$submit_mod
    
    isolate({
      abheatmap()[[4]]
    })
  })
  
  output$ab_table <- renderDataTable({
    input$submit_mod
    
    isolate({
      abheatmap()[[3]]
    })
  })
  
  output$error_mod <- renderText({
    if (!input$submit_mod) return()
    
    isolate("Results not found!")
  })

# TO RENDER THE WHOLE UI
  output$result_mod <- renderUI({
    if (!input$submit_mod) return()
    
    validate(
      need(abheatmap(), "\nSorry, we don't have results for your request :-(")
    )
    
    isolate(
        
      tabsetPanel(
          tabPanel("Heatmap", 
                   downloadButton('download_heatmap', 'Download'),
                   br(), br(),
                   plotOutput("abheatmap", 
                              height=400+nrow(abheatmap()[[2]])*10, 
                              width=450+nrow(abheatmap()[[3]])*15
                              )
                   ),
          tabPanel("Results",
                   downloadButton('download_heat_restable', 'Download'),
                   br(),br(),
                   dataTableOutput("res_table")
                   ),
          
          
          tabPanel("Antibodies", 
                   downloadButton('download_heat_abtable', 'Download'),
                   br(), br(),
                   dataTableOutput("ab_table")),
          
          tabPanel("Peptides",
                   downloadButton('download_heat_peptable', 'Download'),
                   br(), br(),
                   dataTableOutput("pept_table"))
      )
    )
  })
  
  output$download_heat_abtable <- downloadHandler(
  
    filename = function() {paste(input$prot, input$res, input$ptm, input$thres,
                                 'abs.csv', sep='_')},
    content = function(file){
      write.csv(abheatmap()[[3]], file)
    }
  )

  output$download_heat_peptable <- downloadHandler(
  
    filename = function() {paste(input$prot, input$res, input$ptm, input$thres,
                               'pepts.csv', sep='_')},
    content = function(file){
      write.csv(abheatmap()[[2]], file)
    }
  )
  
  output$download_heat_restable <- downloadHandler(
    
    filename = function() {paste(input$prot, input$res, input$ptm, input$thres,
                                 'results.csv', sep='_')},
    content = function(file){
      write.csv(abheatmap()[[4]], file)
    }
  )

  output$download_heatmap <- downloadHandler(
    filename = function() {paste(input$prot, input$res, input$ptm, input$thres,
                                 'heatmap.pdf', sep='_')},
    content = function(file){
      pdf(file, 
          width=(450+nrow(abheatmap()[[3]])*15)/80, 
          height=(400+nrow(abheatmap()[[2]])*10)/80)
      print(abheatmap()[[1]])
      dev.off()
    })

  ## SEARCH: by Antibody

  delf <- reactive({
    if(!input$submit_ab){return()}
    
    isolate(
      delfiaplot(abs, pepts, delfia, mod.list, input$ab.name, input$ab.catnr)
    )
  })
  
  output$error_ab <- renderText({
    if (!input$submit_ab) return()
    
    isolate("Results not found!")
  })
  
  output$delfia <- renderPlot({ 
    if (!input$submit_ab) return()
    
    isolate(
      print(delf()[[1]])     
    )
  })#, width=1000, height=600)
  
  output$delfia_table <- renderDataTable({
    if (!input$submit_ab) return()
    
    isolate(
      delf()[[2]]
    )
  }, options=list(bSortClasses = TRUE))
  
  
  output$download_table <- downloadHandler(
    
    filename = function() {paste(input$ab.name, input$ab.catnr, 'delfia.csv', sep='_')},
    content = function(file){
      write.csv(delf()[[2]], file)
    }
  )
  
  output$download_delfiaplot <- downloadHandler(
    filename = function() {paste(input$ab.name, input$ab.catnr, 'delfia.pdf', sep='_')},
    content = function(file){
      pdf(file, width=15, height=8)
      print(delf()[[1]])
      dev.off()
    })
  
  output$abinfo <- reactive({
    if (!input$submit_ab) return()
    
    isolate(
      HTML(createAbInfoHTML(abs, input$ab.name, input$ab.catnr)))
  })
  
  output$result_ab <- renderUI({    
    
    if (!input$submit_ab) return()
    
    validate(
      need(delf(), "\nYou must choose an antibody ;)")
    )
    
    isolate(
        
      tabsetPanel(
          tabPanel("Delfia",
                   downloadButton("download_delfiaplot", "Download"),
                   br(),
                   br(),
                   plotOutput("delfia", height=600, width="100%")
                   
          ),
          tabPanel("Table",
                   downloadButton('download_table', 'Download'),
                   br(),
                   br(),
                   dataTableOutput("delfia_table")
          ),
          tabPanel("Ab_Info",
                   htmlOutput("abinfo")
                   )
          )
    )
  })
  
  ## COMPARE
  
#   compare <- reactive({
#     
#     if(!input$submit_comp) return()
#     
#     isolate(
#       comparePlot(delfia, pepts, abs, 
#                   input$protX1, input$resX1, input$ptmX1,
#                   input$protX2, input$resX2, input$ptmX2,
#                   input$fold,
#                   input$clonalityX, input$speciesX, input$applicationX,
#                   input$specificityX)) #to add the specificity filter
#   })
#   
# 
#   output$compareplot <- renderPlot({
#     
#     if(!input$submit_comp) return()
#     
#     isolate(
#       print(compare()[[1]]))
#   }
#   )
#   
#   output$compare_table <- renderDataTable({
#     
#     if(!input$submit_comp) return()
#     
#     isolate(
#       compare()[[2]])
#   })
#   
#   output$download_comparetable <- downloadHandler(
#     
#     filename = function() {paste(input$protX1, input$resX1, input$ptmX1, '_vs_',
#                                  input$protX2, input$resX2, input$ptmX2,
#                                  '_', input$fold, 
#                                  '_filt_', input$clonalityX, input$hostX, input$applicationX,
#                                  '.csv', sep='')},
#     content = function(file){
#       write.csv(compare()[[2]], file)
#     }
#   )
#   
#   output$download_compareplot <- downloadHandler(
#     
#     filename = function() {paste(input$protX1, input$resX1, input$ptmX1, '_vs_',
#                                  input$protX2, input$resX2, input$ptmX2,
#                                  '_', input$fold, 
#                                  '_filt_', input$clonalityX, input$hostX, input$applicationX,
#                                  '.pdf', sep='')},
#     content = function(file){
#       pdf(file, 
#           width=(350+20*nrow(unique(compare()[[2]][,c('Pept2_Name', 'Pept2_CatNr')])))/80, 
#           height=((20*nrow(unique(compare()[[2]][,c('Pept1_Name', 'Pept1_CatNr')]))+
#                   350)*nrow(unique(compare()[[2]][,c('Ab_Name', 'Ab_CatNr', 'Ab_Replicate')])))/80
#           )
#       print(compare()[[1]])
#       dev.off()
#     })
#   
#   output$error_comp <- renderText({
#     if (!input$submit_comp) return()
#     
#     isolate("Results not found!")
#   })
#   
#   output$result_comp <- renderUI({    
#     
#     if (!input$submit_comp) return()
#     
#     isolate(
#       if (is.null(compare()[[1]])){ # so far it will never be NULL
#         h2(textOutput("error_comp"),
#            align="left", style="font-size:10px;color:red") 
#       }else{
#         
#         tabsetPanel(
#           tabPanel("Heatmap",
#                    downloadButton("download_compareplot", "Download"),
#                    br(),
#                    br(),
#                    plotOutput("compareplot",
#                               width=350+20*nrow(unique(compare()[[2]][,c('Pept2_Name', 'Pept2_CatNr')])),
#                               height=(350+20*nrow(unique(compare()[[2]][,c('Pept1_Name', 'Pept1_CatNr')])))*
#                                       nrow(unique(compare()[[2]][,c('Ab_Name', 'Ab_CatNr', 'Ab_Replicate')])))
#                    
#           ),
#           tabPanel("Table",
#                    downloadButton("download_comparetable", "Download"),
#                    br(),
#                    br(),
#                    dataTableOutput("compare_table")
#           )
#         )
#       }
#     )
#   })


  # MATERIALS AND METHODS

})









