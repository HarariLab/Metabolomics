# UI function
boxplots_UI <- function(id, metab_meta, Origscale, subject_status) {
  ns <- NS(id)
  tabItem(
    tabName = "boxplots",
    h3("Choose a variable and a metabolite to view reading distributions"),
    fluidRow(
      box(
        dataTableOutput(ns('table')),
        width = 12)),
    fluidRow(
      box(id = ns("boxplot"), 
          title = "Metabolite Distribution", 
          solidHeader = TRUE, 
          status = "primary",           
          plotlyOutput(ns("p")), 
          width = 8),
      # Info on metab selected in table
      box(id = ns("metabinfobox"), 
          title = "Metabolite Details", 
          solidHeader = TRUE, 
          status = "primary", 
          htmlOutput(ns("metabinfo")), 
          width = 4)
    ))
}

# Server function
boxplots <- function(input, output, session, box_var) {
  
  Origscale_clean_transformed_t <- readRDS("data/Origscale_clean_transformed_t.rds")
  
  output$table <- renderDataTable(
    datatable(
      metab_meta,
      selection = list(mode = "single", selected = 1),
      rownames = FALSE,
      options = list(
        autoWidth = TRUE,
        columnDefs = list(list(width = '100px', targets = c(8:16)), list(width = '120px', targets = c(0:7))),
        scrollX = TRUE,
        rowCallback = JS(
          "function(row, data) {",
          "for (i = 8; i < data.length; i++) {",
          "$('td:eq('+i+')', row).html(data[i].toExponential(2));",
          "}",
          "}")
        ),
      filter = "top",
      escape = FALSE
      )
  )
  
  # Prepare data for boxplots
  metab_by_subject <- reactive({
    validate(
      need(
        input$table_rows_selected,
        "Select a metabolite to view plots"
      )
    )
    selected <- input$table_rows_selected
    metab_counts <- data.frame(
      colnames(Origscale_clean_transformed_t)[1:454],
      unlist(c(Origscale_clean_transformed_t[Origscale_clean_transformed_t$BIOCHEMICAL == metab_meta[selected,1],1:454]))
    )
    colnames(metab_counts) <- c("id", "count")
    #metab_counts$count <- log10(metab_counts$count)
    
    merge(subject_status, metab_counts, by = "id")
  })
  
  values <- reactiveValues()
  
  observeEvent(input$table_rows_selected, {
    values$metab_by_subject <- metab_by_subject
  })
  
  # Reading distribution plot ouput
  output$p <- renderPlotly({
    req(values$metab_by_subject)
    # Plot by CACO status
    if (box_var() == "Status") {
      plot_ly(
        data = metab_by_subject()[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD"),],
        y = ~count,
        color = ~CACO,
        colors = "Dark2",
        type = "box",
        boxpoints = "all",
        pointpos = 0,
        text = paste("Transformed Reading: ", 
                     format(metab_by_subject()$count[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD")], digits = 3)),
        hoverinfo = list("median", "text")
      ) %>%
        layout(
          title = metab_meta[input$table_rows_selected,1],
          yaxis = list(title = "Transformed Reading",
                       hoverformat = ".2f"),
          showlegend = FALSE
        ) %>%
        config(displayModeBar = FALSE)
    }
    
    # Plot by Sex
    else if (box_var() == "Sex") {
      plot_ly(
        data = metab_by_subject(),
        y = ~count,
        color = ~SEX,
        colors = "Dark2",
        type = "box",
        boxpoints = "all",
        pointpos = 0,
        text = paste("Log10 Reading: ", 
                     format(metab_by_subject()$count, digits = 3)),
        hoverinfo = list("median", "text")
        ) %>%
        layout(
          title = metab_meta[input$table_rows_selected,1],
          yaxis = list(title = "Log10 Reading"),
          showlegend = FALSE
        ) %>%
        config(displayModeBar = FALSE)
    }
  })
  
  # Information on metabolite selected on volcano plot
  output$metabinfo <- renderText({
    
    metab_name <- as.character(metab_meta[input$table_rows_selected,1])
    
    paste(
      "Metabolite:", metab_name, "</br>",
      "Super Pathway:", metab_meta$`Super Pathway`[metab_meta$Metabolite == metab_name], "</br>",
      "Sub Pathway:", metab_meta$`Sub Pathway`[metab_meta$Metabolite == metab_name], 
      # Link to HMDB, if ID is available
      ifelse(
        !is.na(metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]),
        paste(
          "</br>HMDB ID: ",
          metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]
        ),
        ""
      ), 
      # Link to PubChem, if ID is available
      ifelse(
        !is.na(metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]),
        paste("</br>PubChem ID: ",
              metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]
        ), 
        ""
      )
    )
  })
}




