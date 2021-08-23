# "Volcano" tab: shows volcano plots for comparisons between
# ADAD, CA, CO, and TREM2.  When a point is selected on the
# volcano plot, reading distributions and information on that metabolite
# are shown below the plot.

# UI function
volcano_UI <- function(id, metab_meta, Origscale, subject_status) {
  ns <- NS(id)
  tabItem(
    tabName = "volcano",
    h4("Select a comparison to view volcano plots."),
    h4("Click on plot points to view reading distributions and information 
       on specific metabolites."),
    fluidRow(
      # Volcano plot
      box(plotlyOutput(ns("p")),
          title = "Volcano Plot", 
          solidHeader = TRUE, 
          status = "primary", 
          width = 12)
    ),
    fluidRow(
      # Boxplot for distribution of selected metab by CACO group
      box(id = ns("boxplot"), 
          title = "Metabolite Distribution", 
          solidHeader = TRUE, 
          status = "primary",           
          plotlyOutput(ns("bp")), 
          width = 8),
      # Info on metab selected on volcano plot
      box(id = ns("metabinfobox"), 
          title = "Metabolite Details", 
          solidHeader = TRUE, 
          status = "primary", 
          htmlOutput(ns("metabinfo")), 
          width = 4)
    )
    )
}

# Server function
volcano <- function(input, output, session, compare) {
  
  # Data for boxplots
  Origscale_clean_transformed_t <- readRDS("data/Origscale_clean_transformed_t.rds")
  
  # Read in data formatted for volcano plot
  volcano_data_adadvsco <- readRDS("data/volcano_data_adadvsco.rds")
  volcano_data_trem2vsco <- readRDS("data/volcano_data_trem2vsco.rds")
  volcano_data_cavsco <- readRDS("data/volcano_data_cavsco.rds")
  # volcano_data_adadvsca <- readRDS("data/volcano_data_adadvsca.rds")
  
  # Choose which comparison to display volcano plot for
  volcano_data <- reactive({
    if(req(compare()) == "ADAD vs Control") dat <- volcano_data_adadvsco
    else if (compare() == "TREM2 vs Control") dat <- volcano_data_trem2vsco
    else if (compare() == "AD vs Control") dat <- volcano_data_cavsco
    # else if (compare() == "ADAD vs AD") dat <- volcano_data_adadvsca
    dat$f <- 1
    dat$point_color <- recode(as.character(dat$significant), "TRUE" = "#FF7F00", "FALSE" = "#1F78B4")
    dat
  })
  
  # Store data for animation
  my <- reactive({
    reactiveValues(
      olddata = volcano_data(),
      newdata = volcano_data()
    )
  })
  
  # Volcano plot
  output$p <- renderPlotly(
    plot_ly(
      x = volcano_data()$effects,
      y = -log10(volcano_data()$pvals),
      name = "FDR > 0.05",
      type = "scatter",
      showlegend = FALSE,
      mode = "markers",
      # Hovertext
      text = paste(volcano_data()$metabs,
                   "</br></br>Beta: ",
                   format(volcano_data()$effects, digits = 3, scientific = TRUE),
                   "</br>Q-value: ",
                   format(volcano_data()$pvals_adj, digits = 3, scientific = TRUE)),
      hoverinfo = "text",
      color = ~I(volcano_data()$point_color)) %>%
      # Adding markers for a custom legend.  Technically,
      # the entire volcano plot trace is colored blue,
      # but we need a legend to indicate the meaning of the orange points,
      # so we add traces with orange and blue and relabel.
      # Works better for animation and plotly_click purposes.
      
      # Blue/not significant
      add_markers(x= 0.8, y = 6.5, color = I("#1F78B4"), showlegend = FALSE, hoverinfo = "skip") %>%
      add_annotations(x=0.8, y=6.5, xref = "x", yref = "y", text = "FDR > 0.05", 
                      xanchor = 'left', showarrow = F, xshift = 10) %>%
      # Orange/significant
      add_markers(x= 0.8, y = 7, color = I("#FF7F00"), showlegend = FALSE, hoverinfo = "skip") %>%
      add_annotations(x=0.8, y=7, xref = "x", yref = "y", text = "FDR < 0.05", 
                      xanchor = 'left', showarrow = F, xshift = 10) %>% 
      
      layout(
        title = compare(),
        xaxis = list(title = "Effect (Beta)", range = c(-1, 1)),
        yaxis = list(title = "-log10 p-value", range = c(0, 7.25))
      ) %>%
      # Disable the legend click since our traces do not correspond to the
      # actual legend labels
      onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
      animation_opts(frame=500, transition=500, redraw=TRUE, mode = "next")  %>%
      config(displayModeBar = FALSE)
  )
  
  # Animate points when selection changes
  # List structure in plotlyProxyInvoke is important
  observeEvent(compare(), {
    req(my()$newdata)
    myswitch <- list()
    myswitch$olddata <- my()$newdata # save old data
    myswitch$newdata <- volcano_data() %>% # generate new data
      mutate(f = myswitch$olddata$f + 1) # Advance frame
    plotlyProxy("p", session=session, deferUntilFlush=FALSE) %>%
      plotlyProxyInvoke("animate",
                        # frameOrGroupNameOrFrameList
                        list(
                          data = 
                            list(
                              list(
                                x = myswitch$newdata$effects,
                                y = -log10(myswitch$newdata$pvals),
                                frame = myswitch$newdata$f
                              )
                            ),
                          layout = list(frame = 500, transition = 500)
                        ),
                        # animationAttributes
                        list()
      )# plotlyProxyInvoke
  })
  
  ### Set up data for boxplots
  metab_by_subject <- reactive({
    selected <- event_data(event = "plotly_click")$pointNumber + 1
    metab_counts <- data.frame(
      colnames(Origscale_clean_transformed_t)[1:454],
      unlist(c(Origscale_clean_transformed_t[Origscale_clean_transformed_t$BIOCHEMICAL == volcano_data()[selected,1], 1:454]))
    )
    colnames(metab_counts) <- c("id", "count")
    #metab_counts$count <- log10(metab_counts$count)
    
    merge(subject_status, metab_counts, by = "id")
  })
  
  output$bp <- renderPlotly({
    validate(
      need(
        event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
        "Select a metabolite on the volcano plot above to view distribution among statuses."
      )
    )
    # Plot by CACO status
    # CO
    plot_ly(
      data = metab_by_subject()[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD"),],
      y = ~count,
      color = ~CACO,
      type = "box",
      boxpoints = "all",
      pointpos = 0,
      text = paste("Transformed Reading: ",
                   format(metab_by_subject()$count[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD")], digits = 3)),
      hoverinfo = list("median", "text"),
      colors = "Dark2"
    ) %>%
      layout(
        title = volcano_data()[event_data(event = "plotly_click")$pointNumber + 1,1],
        yaxis = list(title = "Transformed Reading",
                     hoverformat = ".2f"),
        showlegend = FALSE
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  # Information on metabolite selected on volcano plot
  output$metabinfo <- renderText({
    validate(
      need(
        event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
        "Select a metabolite on the volcano plot above to view distribution among statuses."
      )
    )
    metab_name <- volcano_data()[event_data(event = "plotly_click")$pointNumber + 1, 1]
    
    paste(
      "<b>Metabolite</b>: ", metab_name, "</br>",
      "<b>Super Pathway</b>: ", metab_meta$`Super Pathway`[metab_meta$Metabolite == metab_name], "</br>",
      "<b>Sub Pathway:</b> ", metab_meta$`Sub Pathway`[metab_meta$Metabolite == metab_name], 
      # Link to HMDB, if ID is available
      ifelse(
        !is.na(metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]),
        paste(
          "</br><b>HMDB ID:</b> ",
          metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]
        ),
        ""
      ), 
      # Link to PubChem, if ID is available
      ifelse(
        !is.na(metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]),
        paste("</br><b>PubChem ID:</b> ",
              metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]
        ), 
        ""
      ),
      ifelse(
        !is.na(metab_name),
        paste(
          "</br><b>ADAD vs CO effect (q):</b> &nbsp&nbsp&nbsp", 
          format(volcano_data_adadvsco$effects[volcano_data_adadvsco$metabs == metab_name], digits = 3, scientific = TRUE), 
          " (", format(volcano_data_adadvsco$pvals_adj[volcano_data_adadvsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
          "<b>AD vs CO effect (q)</b>: &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp", 
          format(volcano_data_cavsco$effects[volcano_data_cavsco$metabs == metab_name], digits = 3, scientific = TRUE), 
          " (", format(volcano_data_cavsco$pvals_adj[volcano_data_cavsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
          "<b>TREM2 vs CO effect (q)</b>: &nbsp", 
          format(volcano_data_trem2vsco$effects[volcano_data_trem2vsco$metabs == metab_name], digits = 3, scientific = TRUE), 
          " (", format(volcano_data_trem2vsco$pvals_adj[volcano_data_trem2vsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
          # "<b>ADAD vs AD effect (q)</b>: &nbsp&nbsp&nbsp", 
          # format(volcano_data_adadvsca$effects[volcano_data_adadvsca$metabs == metab_name], digits = 3, scientific = TRUE), 
          # " (", format(volcano_data_adadvsca$pvals_adj[volcano_data_adadvsca$metabs == metab_name], digits = 3, scientific = TRUE), ")",
          sep = ""),
        ""
      )
    )
  })
}
