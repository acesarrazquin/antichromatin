shinyUI(navbarPage(
  
  theme="theme.css",
  title=includeHTML("html/title.html"),

  windowTitle="antiCHROMATIN",
  id="navbar",
  header="",
  inverse=FALSE,
  collapsable=TRUE,
  fluid=TRUE,
  responsive=TRUE,
  
  
  tabPanel(" ",
          tags$div(class="front",
                   includeHTML("html/front.html")
           
           )),
    
  tabPanel("Home", 
           
           id="home",
           #HTML("<img src='CeMMlogo.jpg'/>")),

           
           tags$div(class="nofront",

            fluidRow(
             column(2),
             
             column(8,
                    
                  tags$style(type="text/css", "p{text-align: justify;
                                                font-size: small;
                                                margin:10px 0px 10px 0px;}"),
                  
                  br(), br(),
                  #plotOutput("abclassification")),
                  includeHTML("html/home.html")
             )),

             column(2)
           )),

 

  
  tabPanel("Search antibodies",
      id="search",
     
      tags$div(class="nofront",
               fluidRow(
          column(3,
              style="height:1000px",
                  br(), br(),
                  
              tags$div(
                "Search antibodies:",
                style="text-align:center;
                        font-size:20px;",
                tags$div(
                  style="float:right;",
                  actionLink("help_mod", "", icon= icon("info-circle")))
              ),
              
              #possible to use includeHTML with external html file
              conditionalPanel(
                condition = "input.help_mod % 2 != 0",
                
                includeHTML("html/searchab_help.html")
              ),
              
              br(),
              
              wellPanel(
                      tags$style(type="text/css", "label.radio { display: inline-block; }",
                                            ".radio input[type=\"radio\"] { float: none; }"),
                                 
                                 # style="background-color:#FFFFFF;border-color:#E8E8E8;",
                                 
                                 
                      selectInput("prot", "Protein:", choices=" ", selected=" "),
                                 
                      selectInput("res", "Residue:", choices="", multiple=TRUE),
                                 
                      selectInput("ptm", "Post-translational modifications:", choices="", multiple=TRUE),
                                 
                      #             h6(textOutput("ptm.error"),
                      #                        align="left", style="font-size:10px;color:red"),
                                 
                      sliderInput("threshold", "Signal threshold:",
                                             min=0, max=1, value = 0.25, step=0.05
                      ),
                                 
                      br(),
                                 
                      checkboxInput("target", 
                          "only Abs targeting modification",
                          value=FALSE),
                                 
                      br(),
                                 
                      actionLink("show_filters", "", icon("plus-square")),
                                 
                      br(),
                                 
                      conditionalPanel(
                          condition = "input.show_filters % 2 != 0",
                          br(),
                                   
                          radioButtons("clonality", "Ab clonality:",
                                        choices=c("Any", "Monoclonal", "Polyclonal")
                          ),
                                   
                          br(),
                                   
                          checkboxGroupInput("species", "Ab host species:", 
                                             choices=""
                          ),
                          
                          br(),
                                   
                          selectInput("application", "Application:",
                                       choices="", multiple=TRUE
                          ),
                                   
                          br(),
                                   
                          radioButtons("specificity", "Ab specificity filter:",
                                        choices=c("High", "Moderate", "Low"),
                                        selected="Low"
                          )
                          )
                          ),
                                                            
                     actionButton("submit_mod", "submit"),
                                                         
                     br(), br(),
                               
                     conditionalPanel(
                        condition = "!(input.submit_mod == 0)",
                        br(), br(),
                        selectInput("chooseab", "Choose antibody:", 
                                             choices="", multiple=FALSE),
                        actionButton("submit_abchoice", "submit"),
                        br())
                     ),
          column(9,
                 br(), br(),        
                 uiOutput('result_mod'),
                 br()
            
#           column(8,
#               br(), br(),
#                                
#               #uiOutput("result_mod")
#               conditionalPanel(
#                 condition = "!(input.submit_mod == 0)",
#                 
#                 tabsetPanel(
#                   tabPanel("Heatmap", 
#                       downloadButton('download_heatmap', 'Download'),               
#                       br(), br(),         
#                       #ggvisOutput("ggvis")
#                       plotOutput("abheatmap")
# 
#                   ),
#                   tabPanel("Results",
#                       downloadButton('download_heat_restable', 'Download'),
#                       br(),br(),
#                       dataTableOutput("res_table")
#                    ),
#       
#                   tabPanel("Antibodies", 
#                       downloadButton('download_heat_abtable', 'Download'),
#                       br(), br(),
#                       dataTableOutput("ab_table")),
#                             
#                   tabPanel("Peptides",
#                       downloadButton('download_heat_peptable', 'Download'),
#                       br(), br(),
#                       dataTableOutput("pept_table"))
#                 )
#              )
          )
      )
  )),
             
  tabPanel("Antibody profiles",
      id="ab",
      tags$div(class="nofront",
      fluidRow(
          column(3,
              style="height:1000px",
              br(), br(),
              
              tags$div(
                "Antibody profile",
                style="text-align:center;
                        font-size:20px;",
                tags$div(
                  style="float:right;",
                  actionLink("help_ab", "", icon= icon("info-circle", lib = "font-awesome")))
              ),
              
              conditionalPanel(
                condition = "input.help_ab % 2 != 0",
                
                includeHTML("html/abprofile_help.html")
              ),
              
              br(),
              
              wellPanel(
                  selectInput("ab.name", "Name:", "", selected=" "),
                                 
                  selectInput("ab.catnr", "CatNr:", "",selected=" ")                        
              ),
                               
             actionButton("submit_ab", "submit")),
                        
          column(9,
              br(), br(),
              uiOutput('result_ab')


#               conditionalPanel(
#                 condition = "!(input.submit_ab == 0)",
#                 
#                 br(),br(),
#                 
#                 
#                 tabsetPanel(
#                  
#                   tabPanel("Delfia",
#                            downloadButton("download_delfiaplot", "Download"),
#                            br(),
#                            br(),
# 
#                            plotOutput("delfia", width="100%", height="auto")
#                            
#                   ),
#                   tabPanel("Table",
#                            downloadButton('download_table', 'Download'),
#                            br(),
#                            br(),
#                            dataTableOutput("delfia_table")
#                   ),
#                   tabPanel("Ab_Info",
#                            htmlOutput("abinfo")
#                   )
#                 )
              )
        ))
  ),
             
#    tabPanel("Compare", # Same as "by Modification" but two times, difference in detection, filters
#                                 # Select if only that one above threshold??
#                       id="comp",
#                       
#                       fluidRow(
#                         column(3,
#                           br(), br(),
#                           
#                           wellPanel(
#                              "Modification 1",
#                              
#                              selectInput("protX1", "protein:", choices=" ", selected=" "),
#                                  
#                              selectInput("resX1", "residue:", choices="", multiple=TRUE),
#                                  
#                              selectInput("ptmX1", "post-translational modifications:", choices="", multiple=TRUE),
#                                  
#                              br()),
#                           
#                           wellPanel(
#                             "Modification 2",
#                              
#                              selectInput("protX2", "protein:", choices=" ", selected=" "),
#                                  
#                              selectInput("resX2", "residue:", choices="", multiple=TRUE),
#                                  
#                              selectInput("ptmX2", "post-translational modifications:", choices="", multiple=TRUE),
#                                  
#                              br()),
#                           
#                           
#                           sliderInput("fold", "Signal difference (x-fold):",
#                                         min=0, max=20, value = 5, step=0.5),
#                         
#                           
#                           actionButton("show_filtersX", "more filters"),
#                           
#                           br(),
#                           
#                           conditionalPanel(
#                               condition = "input.show_filtersX % 2 != 0",
#                               br(),
#                             
#                               wellPanel(
#                                 
#                                 radioButtons("clonalityX", "Ab clonality:",
#                                          choices=c("Any", "Monoclonal", "Polyclonal")
#                                 ),
#                             
#                                 br(),
#                             
#                                 checkboxGroupInput("speciesX", "Ab host:", 
#                                                choices=""),
#                                 br(),
#                             
#                                 selectInput("applicationX", "Application:",
#                                           choices="", multiple=TRUE
#                                 ),
#                                 
#                                 br(),
#                                 
#                                 radioButtons("specificityX", "Ab specificity filter:",
#                                              choices=c("High", "Moderate", "Low"), selected="Low")
#                             
#                               )),
#                           
#                           br(),
#                           
#                           actionButton("submit_comp", "submit"),
#                           
#                           br(), br(),
#                           
#                           conditionalPanel(
#                             condition = "!(input.submit_comp == 0)",
#                             br(), br(),
#                             
#                             selectInput("chooseabX", "Choose antibody:", 
#                             choices="", multiple=FALSE),
#                             
#                             actionButton("submit_abchoiceX", "submit"),
#                             br())  
#                         ),
#                         
#                         column(9,
#                                br(), br(),
#                                uiOutput('result_comp')
#                         )
#                       ))
#                                
#   ),
  
  tabPanel("Materials & Methods",
           tags$div(class="nofront",
           fluidRow(
             column(2),
             column(8,
                    br(),br(),
                    includeHTML("html/materials.html")
             ),
             column(2))
  )),
           
  tabPanel("About us",
           tags$div(class="nofront",
           fluidRow(
             column(2),
             column(8,
                    br(), br(),
                    includeHTML('html/about.html')
             ),
             column(2))
  ))

))
