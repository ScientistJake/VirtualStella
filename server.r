options(shiny.sanitize.errors = FALSE)
library(plotly)
library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(shinyBS)
library(shinythemes)

load(file="data.RData")

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}
server <- function(input, output, session) {
  
  #output$images <- renderUI({})
  #outputOptions(output, "images", suspendWhenHidden = FALSE)
  
  #make a spherical embryo
  set.seed(101)
  n <- 2000
  theta <- runif(n,0,2*pi)
  u <- runif(n,-1,1)
  x <- sqrt(1-u^2)*cos(theta)
  y <- sqrt(1-u^2)*sin(theta)
  z <- u
  
  #Make some PMCs 'inside'
  nPMC <- 64
  thetaPMC <- runif(nPMC,0,2*pi)
  uPMC <- runif(nPMC,-1,1)
  xPMC <- (sqrt(1-uPMC^2)*cos(thetaPMC))
  yPMC <- (sqrt(1-uPMC^2)*sin(thetaPMC))
  zPMC <- uPMC
  
  embryo <- data.frame(x=x,y=y,z=z)
  #PMC <- data.frame(x=xPMC,y=yPMC,z=zPMC, territory = 4)
  
  #designate territories
  embryo$territory <- ifelse(embryo$z>0.91, 1,0)
  embryo$territory <- ifelse(embryo$z>0.6501 & embryo$z< 0.9, 2, embryo$territory)
  embryo$territory <- ifelse(embryo$z< 0.40 & embryo$z< 0.65,3,embryo$territory)
  
  #add the coordinates
  #embryo <- rbind(embryo,PMC)
  
  #This will subset the coordinates based on the slicer
  embryoPlot <- reactive({
    embryo[which(embryo$x < input$slicer),]
  })
  
  output$plot <- renderPlotly({
    p <- plot_ly(embryoPlot(), 
                 x = ~x, y = ~y, z = ~z, 
                 color = ~territory, 
                 colors = c('#FF7070','#4AC6B7','#BF382A','#0C4B8E'),
                 type = "scatter3d",
                 mode="markers") %>%
      add_markers()
  })
  
  output$click <- renderUI({
    d <- event_data("plotly_click")
    if (!is.null(d)){
      territory <- embryoPlot()$territory[d$pointNumber + 1]
      if (territory == 0){
        return(h2("External Ring"))
      } else if (territory == 1){
        return(h2("Central Domain"))
      } else if (territory == 2){
        return (h2("Central Ring"))
      } else if (territory == 3){
        return (h2("Ectoderm"))
      } else if (territory == 4){
        return (h2("Neural Ectoderm"))
      }
    } else {
      return(h2("Click a territory"))
    }
  })
  
  territory_table <- reactive({
    d <- event_data("plotly_click")
    if (!is.null(d)){
      territory <- embryoPlot()$territory[d$pointNumber + 1]
      if (territory == 0){
        territory_table <- data[which(data$ish_blastula=="external_ring"),c(1,2,3)]
      } else if (territory == 1){
        territory_table <- data[which(data$ish_blastula=="central_domain"),c(1,2,3)]
      } else if (territory == 2){
        territory_table <- data[which(data$ish_blastula=="central_ring"),c(1,2,3)]
      } else if (territory == 3){
        territory_table <- data[which(data$ish_blastula=="ectoderm"),c(1,2,3)]
      } else if (territory == 4){
        territory_table <- data[which(data$ish_blastula=="neural_ectoderm"),c(1,2,3)]
      }
      territory_table
    }
  })
  
  output$table <- DT::renderDataTable({
    req(territory_table)
    out <- territory_table()
    out$Actions = shinyInput(actionButton, length(out[,1]), 'button_', label = "Show Expression!", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
    out[c(1,4,3,2)]
  },server = FALSE, escape = FALSE, selection = 'none')
  
  selectedgene <- eventReactive(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    territory_table()[selectedRow,]
  })
  
  output$e_plot <- renderPlot({
    req(selectedgene())
    
    gene <- unlist(selectedgene()[1])
    
    embryo_plot <- embryo_line[which(embryo_line$gene==gene),]
    fisher_plot <- fisher[which(fisher$gene==gene),]
    helm_plot <- helm[which(helm$gene==gene),]
    warner_plot <- warner[which(warner$gene==gene),]
    
    q<- ggplot(embryo_plot, aes(x=as.numeric(as.character(embryo_plot$variable)), y=value, colour=gene)) + geom_line() +
      scale_shape_discrete(solid=F, name='Dataset') +
      geom_point(data = warner_plot, aes(x=as.numeric(as.character(warner_plot$variable)), y=value, colour=gene, shape='Warner et al. (2017)'), size=2) +
      geom_point(data = fisher_plot, aes(x=as.numeric(as.character(fisher_plot$variable)), y=value, colour=gene, shape='Fischer et al. (2014)'),size=2) +
      geom_point(data = helm_plot, aes(x=as.numeric(as.character(helm_plot$variable)), y=value, colour=gene, shape='Helm et al. (2103)'),size=2) +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
      scale_x_continuous(minor_breaks = NULL, breaks=c(0,6,12,24,48,72,96,120,144,168,192,120,240))+
      ylab("Log2 CPM") +
      xlab("Hours Post Fertilization")
    return(q)
  })
  
  imageslist<- reactive({
    req(selectedgene())
    gene <- as.character(unlist(selectedgene()[1]))
    imagepaths <- imagedata[which(imagedata$gene_common == gene),]
    outputimages <- lapply(imagepaths[,2],function(x){
      title <- paste0("<p>",gsub(".jpg","",x),"</p>")
      link <- paste0("<img src='images/",x,"'></img>")
      rbind(title,link)
    })
    return(outputimages)
  })
  
  output$images<- renderUI({
    HTML(unlist(imageslist()))
  })
  
  observeEvent(input$select_button,{ 
    toggleModal(session, "expression", "open") 
  })
  
  output$myText <- renderUI({
    h1(selectedgene()[1])
  })
  
}