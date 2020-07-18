library(shiny)
library(plotly)
library(RColorBrewer)
library(DT)
library(RRNA)
library(shinydashboard)
library(shinythemes)


# cache computation of the correlation matrix
#load("./data/dataForShiny2.RData")
#load("./data/clusterPositionsWithStructures.RData")
myCol = c("#000000","#000000","#000000","#000000",colorRampPalette(brewer.pal(8,"YlOrRd"))(40))
server <- function(input, output, session) {
  
  
  
  ######################################################
  # Heatmaps - sorting out the zoom of both
  ######################################################
  output$heatSample <- renderPlotly({
    plot_ly(source = "heat_plot") %>%
      add_heatmap(x = row.names( combinedMatListScaleNoReps[["sars"]][["sample"]]), 
                  y = colnames( combinedMatListScaleNoReps[["sars"]][["sample"]]),
                  z = log2( combinedMatListScaleNoReps[["sars"]][["sample"]]+1),
                  colors = myCol)
  })
  
  output$heatControl <- renderPlotly({
    
    # if there is no click data, render the normal plot
    clickData <- event_data("plotly_relayout", source = "heat_plot")
    clickDataOld = clickData
    
    if (is.null(clickData)){
      clickData = list()
      clickData$`xaxis.range[0]` = 0
      clickData$`xaxis.range[1]` = 3000
      clickData$`yaxis.range[0]` = 0
      clickData$`yaxis.range[1]` = 3000
      # return(NULL)
    }else{
      # must change the coordinates to indexes
      factor = 3000 /29883
      clickData$`xaxis.range[0]` = round(clickData$`xaxis.range[0]` * factor)
      clickData$`xaxis.range[1]` = round(clickData$`xaxis.range[1]` * factor)
      clickData$`yaxis.range[0]` = round(clickData$`yaxis.range[0]` * factor)
      clickData$`yaxis.range[1]` = round(clickData$`yaxis.range[1]` * factor)
    }
    print(clickDataOld)
    print( clickData)
    control = combinedMatListScaleNoReps[["mers"]][["sample"]][clickData$`yaxis.range[0]`:clickData$`yaxis.range[1]`,
                                                               clickData$`xaxis.range[0]`: clickData$`xaxis.range[1]` ]
    # scatterplot with fitted line
    plot_ly() %>%
      add_heatmap(x = row.names(control),
                  y = colnames(control),
                  z = log2(control+1),
                  colors = myCol)
  })
  
  ######################################################
  # END - Heatmaps
  ######################################################
  
  
  ######################################################
  # Render Data Table 
  ######################################################
  output$clusterTable <- renderDataTable({
    
    #get zoom parameters
    clickData <- event_data("plotly_relayout", source = "heat_plot")
    clickDataOld = clickData
    if (is.null(clickData)){
      clickData = list()
      clickData$`xaxis.range[0]` = 0
      clickData$`xaxis.range[1]` = 30000
      clickData$`yaxis.range[0]` = 0
      clickData$`yaxis.range[1]` = 30000
    }
    # return(NULL)
    t = rbind.data.frame(clusterPositionsListTrimmed[["sars"]][["sample"]],
                         clusterPositionsListTrimmed[["sars"]][["control"]],
                         stringsAsFactors = F)
    
    row.names(t) = NULL
    
    t = t[t$ls >= clickData$`xaxis.range[0]` & t$le <= clickData$`xaxis.range[1]` & t$rs >= clickData$`yaxis.range[0]` & t$re <= clickData$`yaxis.range[1]`,]  
    
    t[order(t$size.x, decreasing = T),]
    
  })
  
  ######################################################
  # END - data table 
  ######################################################
  
  ######################################################
  # Sequence
  ######################################################
  
  
  output$Sequences <- renderUI({
    clickData <- event_data("plotly_relayout", source = "heat_plot")
    clickDataOld = clickData
    if (is.null(clickData)){
      clickData = list()
      clickData$`xaxis.range[0]` = 0
      clickData$`xaxis.range[1]` = 30000
      clickData$`yaxis.range[0]` = 0
      clickData$`yaxis.range[1]` = 30000
    }
    # return(NULL)
    t = rbind.data.frame(clusterPositionsListTrimmed[["sars"]][["sample"]],
                         clusterPositionsListTrimmed[["sars"]][["control"]],
                         stringsAsFactors = F)
    row.names(t) = NULL
    t = t[t$ls >= clickData$`xaxis.range[0]` &
            t$le <= clickData$`xaxis.range[1]` &
            t$rs >= clickData$`yaxis.range[0]` &
            t$re <= clickData$`yaxis.range[1]`,]  
    t = t[order(t$size.x, decreasing = T),]
    
    index = input$clusterTable_row_last_clicked
    index = t[index,]
    print(index$ls)
    #print(index)
    #now get the sequence that relates to the chosen cluster
    seq = rnaRefs[["sars"]]
    seq1 = paste(seq[[1]][index$ls:index$le], collapse = "")
    name1 = paste(">",names(seq),"-",index$id,"-",index$ls,"-",index$le,sep="")
    seq2 = paste(seq[[1]][index$rs:index$re], collapse = "")
    name2 = paste(">",names(seq),"-",index$id,"-",index$rs,"-",index$re,sep="")
    seq3 = paste(seq[[1]][index$ls:index$re], collapse = "")
    name3 = paste(">",names(seq),"-",index$id,"-",index$ls,"-",index$re,"-","FULLSEQ",sep="")
    
    HTML(paste(name1, toupper(seq1),name2,toupper(seq2),name3,toupper(seq3), sep = '<br/>'))
    
  })
  
  
  output$structure <- renderImage({
    clickData <- event_data("plotly_relayout", source = "heat_plot")
    clickDataOld = clickData
    if (is.null(clickData)){
      clickData = list()
      clickData$`xaxis.range[0]` = 0
      clickData$`xaxis.range[1]` = 30000
      clickData$`yaxis.range[0]` = 0
      clickData$`yaxis.range[1]` = 30000
    }
    # return(NULL)
    t = clusterPositionsListTrimmedSarsCombinedWithStructures
    
    row.names(t) = NULL
    t = t[t$ls >= clickData$`xaxis.range[0]` &
            t$le <= clickData$`xaxis.range[1]` &
            t$rs >= clickData$`yaxis.range[0]` &
            t$re <= clickData$`yaxis.range[1]`,]  
    t = t[order(t$size.x, decreasing = T),]
    
    index = input$clusterTable_row_last_clicked
    index = t[index,]
    print(index$ls)
    
    outfile <- tempfile(fileext = '.vienna')
    writeLines(c(">x",paste(index[,"seq1new"],index[,"seq2new"], sep ="  "),
                 sub("&","  ",index[,"vienna"])), "./programs/vienna.vienna")
    
    print(c(">x",paste(index[,"seq1new"],index[,"seq2new"], sep =""),
            sub("&","",index[,"vienna"])))
    annotString = paste("-annotations ",'"',(index[,"ls"]+9),":anchor=",10,';',index[,"re"],":anchor=",nchar(sub("&","  ",index[,"vienna"])),'"', sep = "")
    command = paste("java -cp ./programs/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -i ./programs/vienna.vienna -o ./programs/output.svg", annotString)
    print(command)
    x = system(command,intern = T)
    
    list(src = "./programs/output.svg",
         contentType = "image/svg+xml",
         height = 1600,
         width = 1200)
    
  }, deleteFile = F)
  
  
}

#shinyApp(ui = ui, server = server)

