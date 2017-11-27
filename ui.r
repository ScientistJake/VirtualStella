library(plotly)
library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(shinyBS)
library(shinythemes)

ui <- navbarPage(theme = shinytheme("flatly"),id = "inTabset", title="Virtual Embryo",
                 tabPanel("Home",
                          mainPanel(
                            h1("Virtual Embryo!"),
                            plotlyOutput("plot", height = "300px"),
                            # verbatimTextOutput("hover"),
                            sliderInput("slicer","Slice Embryo!",min=-1, max=1,value = 1,step = 0.1)
                          ),
                          mainPanel(
                            uiOutput("click"),
                            h4("Gene's expressed in this territory (click one for more info):"),
                            dataTableOutput("table"),
                            #verbatimTextOutput("myText"),
                            bsModal('expression','Expression',trigger="",size='large',
                                    uiOutput("myText"),
                                    h2("Temporal Expression"),
                                    plotOutput("e_plot"),
                                    h2("Spatial Expression"),
                                    uiOutput("images")
                            )
                          )
                 ),
                 tabPanel("About",
                          mainPanel(
                            h2("Data sources"),
                            p("Write this later.")
                          )
                 )
)