library(shiny)

shinyUI(navbarPage(theme="bootstrap_flatly.css",
  title="Histone Modification Ab Database",
  id="navbar",
  header="",
  inverse=FALSE,
  collapsable=TRUE,
  fluid=TRUE,
  responsive=TRUE,#theme="", # I can define the bootstrap stylesheet theme in www/
  

  
  
  tabPanel("Home", 
           id="home",
           #HTML("<img src='CeMMlogo.jpg'/>")),
           column(3),
           column(12,
                  tags$style(type="text/css", "p{text-align: justify;}"),
                  #imageOutput('cemmlogo', width="100px", height="200px"))
                  br(),
                  #plotOutput("abclassification")),

                  includeHTML("html/home.html"),
                  #img(src="delfia.png", height = 400, width = 300),
                  includeHTML("html/abclassification.html"),
                  br(),
                  #imageOutput("abclassification", height=600, width=800), #make it reactive
                  img(src="abclassification.jpg", height=600, width=800),
                  br(),br(),
                  includeHTML("html/heatmap.html"),
                  br(),
                  img(src="whole_heatmap.jpg", height=800, width=600),
                  br(),br(),
                  includeHTML("html/heatmap123.html"),
                  br(),
                  img(src="heatmap123.jpg", height=800, width=800),
                  br(),br()),
           column(2)),
 

  
  navbarMenu("Search",
             id="search",
             
             tabPanel("by Modification",
                      id="mod",
                      #value="tab.mod",
                      
                      fluidRow(
                        column(3,
                               style="height:1000px",
                               br(), br(),
                               wellPanel(
                                 tags$style(type="text/css", "label.radio { display: inline-block; }",
                                            ".radio input[type=\"radio\"] { float: none; }"),
                                 
                                 # style="background-color:#FFFFFF;border-color:#E8E8E8;",
                                 
                                 
                                 selectInput("prot", "protein:", choices=" ", selected=" "),
                                 
                                 selectInput("res", "residue:", choices="", multiple=TRUE),
                                 
                                 selectInput("ptm", "post-translational modifications:", choices="", multiple=TRUE),
                                 
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
                                 
                                 actionButton("show_filters", "more filters"),
                                 
                                 br(),
                                 
                                 conditionalPanel(
                                   condition = "input.show_filters % 2 != 0",
                                   br(),
                                   
                                   radioButtons("clonality", "Ab clonality:",
                                                choices=c("Any", "Monoclonal", "Polyclonal")
                                   ),
                                   
                                   br(),
                                   
                                   checkboxGroupInput("species", "Ab host:", 
                                                      choices=""),
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
                        )
                      )),
             
             tabPanel("by Antibody",
                      id="ab",
                      fluidRow(
                        column(3,
                               style="height:1000px",
                               br(), br(),
                               wellPanel(
                                 #             #id="ab",
                                 #             #value="tab.ab",
                                 selectInput("ab.name", "Name:", "", selected=" "),
                                 
                                 selectInput("ab.catnr", "CatNr:", "",selected=" ")
                            

                               ),
                               
                               actionButton("submit_ab", "submit")),
                        
                        column(9,
                               br(), br(),
                               uiOutput('result_ab')
                        )
                      )),
             
             tabPanel("Compare", # Same as "by Modification" but two times, difference in detection, filters
                                # Select if only that one above threshold??
                      id="comp",
                      
                      fluidRow(
                        column(3,
                          br(), br(),
                          
                          wellPanel(
                             "Modification 1",
                             
                             selectInput("protX1", "protein:", choices=" ", selected=" "),
                                 
                             selectInput("resX1", "residue:", choices="", multiple=TRUE),
                                 
                             selectInput("ptmX1", "post-translational modifications:", choices="", multiple=TRUE),
                                 
                             br()),
                          
                          wellPanel(
                            "Modification 2",
                             
                             selectInput("protX2", "protein:", choices=" ", selected=" "),
                                 
                             selectInput("resX2", "residue:", choices="", multiple=TRUE),
                                 
                             selectInput("ptmX2", "post-translational modifications:", choices="", multiple=TRUE),
                                 
                             br()),
                          
                          
                          sliderInput("fold", "Signal difference (x-fold):",
                                        min=0, max=20, value = 5, step=0.5),
                        
                          
                          actionButton("show_filtersX", "more filters"),
                          
                          br(),
                          
                          conditionalPanel(
                              condition = "input.show_filtersX % 2 != 0",
                              br(),
                            
                              wellPanel(
                                
                                radioButtons("clonalityX", "Ab clonality:",
                                         choices=c("Any", "Monoclonal", "Polyclonal")
                                ),
                            
                                br(),
                            
                                checkboxGroupInput("speciesX", "Ab host:", 
                                               choices=""),
                                br(),
                            
                                selectInput("applicationX", "Application:",
                                          choices="", multiple=TRUE
                                ),
                                
                                br(),
                                
                                radioButtons("specificityX", "Ab specificity filter:",
                                             choices=c("High", "Moderate", "Low"), selected="Low")
                            
                              )),
                          
                          br(),
                          
                          actionButton("submit_comp", "submit"),
                          
                          br(), br(),
                          
                          conditionalPanel(
                            condition = "!(input.submit_comp == 0)",
                            br(), br(),
                            
                            selectInput("chooseabX", "Choose antibody:", 
                            choices="", multiple=FALSE),
                            
                            actionButton("submit_abchoiceX", "submit"),
                            br())  
                        ),
                        
                        column(9,
                               br(), br(),
                               uiOutput('result_comp')
                        )
                      ))
                               
  ),
  
  tabPanel("Materials&Methods",
           column(2),
           column(12,
            includeHTML("www/materials.htm"),
            br(), br()
           ),
           column(2)),
           
  tabPanel("About us",
           column(2),
           column(12,
                  br(), br(),
                  tags$p('This is us.')
                  ),
           column(2))
))
