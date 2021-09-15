library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(DT)
library(shinyjs)
library(htmlwidgets)
library(RColorBrewer)
library(kableExtra)
library(dplyr)
library(waiter)

# Load in modules for tabs and sidebar
module_list <- sprintf("modules/%s", c("home.R", "boxplots.R", "volcano.R"))
for (module in module_list) source(module, local=TRUE)

# Data for boxplots tab table
metab_meta <- readRDS("data/metab_meta.rds")

# For boxplots
subject_status <- readRDS("data/subject_status.rds")

# Sidebar
sidebar <- dashboardSidebar(
  useShinyjs(),
  tags$head(tags$style(HTML('.logo {
                              background-color: #a51417 !important;
                              }
                              .navbar {
                              background-color: #a51417 !important;
                              }
                              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              border-left-color: #a51417;
                              }
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
                              border-left-color: #a51417;
                              }  
                              .box.box-solid.box-primary>.box-header {
                              color:#000000;
                              background:#c8c8c8
                              }

                              .box.box-solid.box-primary{
                              border-bottom-color:#c8c8c8;
                              border-left-color:#c8c8c8;
                              border-right-color:#c8c8c8;
                              border-top-color:#c8c8c8;
                              }
                              .control-sidebar-dark {
                              color: #b8c7ce;
                              padding-left: 10px;
                              padding-top: 60px;
                              }
                            
                              .checkbox label, .radio label {
                              padding-left:30px
                              }'))),
  sidebarMenu (
    id = "tabs",
    menuItem("Home", icon = icon("home"), tabName = "home"),
    menuItem("Reading Distributions", icon = icon("chart-bar"), tabName = "boxplots"),
    menuItem("Volcano Plots", icon = icon("chart-bar"), tabName = "volcano")
  )
)

# Tabs
body <- dashboardBody(
  # Resize plots when changing tabs (prevents plots running over
  # boxes when right sidebar expands)
  tags$script('
      // Bind function to sidebar buttons
      $(".sidebar-menu").on("click",function(){
        $(window).trigger("resize"); // Trigger resize event
      })'
  ),
  tabItems(
    home_UI("home_ns"),
    boxplots_UI("boxplot_ns", metab_meta, subject_status),
    volcano_UI("volcano_ns", subject_status)
  )
)

rightsidebar <- dashboardControlbar(
  id = "controlbar",
  uiOutput("right_sidebar")
)

# UI functions for tabs
ui <- dashboardPage(
  title = "ONTIME Metabolomics",
  header = dashboardHeader(
    title = "ONTIME",
    controlbarIcon = shiny::icon("bars")),
  sidebar = sidebar,
  body = body,
  controlbar = rightsidebar,
  preloader = list(html = spin_solar(), color = "#333e48"),
  footer = NULL,
  skin = "blue"
)

# Server functions for tabs
server <- function(input, output, session) {
  
  # Right sidebar should be shown for the volcano and boxplot tabs only
  # and have inputs that are specific to which tab is selected.
  observe({
    # Custom right sidebar for volcano tab
    if (req(input$tabs) == "volcano") {
      # Show right sidebar toggle button
      # show(selector = "body > div.wrapper > header > nav > div:nth-child(4) > ul")
      hide(selector = "body > div.wrapper > header > nav > div:nth-child(4) > ul")
      # Show right sidebar
      addClass(selector = "body", class = "control-sidebar-open")
      # Inputs for which comparison to show on volcano plot
      output$right_sidebar <- renderUI(
        controlbarMenu(
          id = "volcano_sidebar",
        controlbarItem(
          "Options",
          icon = icon("gears"),
          radioButtons("comparison", 
                       "Choose a comparison", 
                       c("ADAD vs Control", 
                         "AD vs Control",
                         "TREM2 vs Control")
          )
        )
      )
      )
    }    
    
    # Custom right sidebar for boxplots tab
    else if (req(input$tabs) == "boxplots") {
      # Show right sidebar toggle button
      # show(selector = "body > div.wrapper > header > nav > div:nth-child(4) > ul")
      hide(selector = "body > div.wrapper > header > nav > div:nth-child(4) > ul")
      # Show right sidebar
      addClass(selector = "body", class = "control-sidebar-open")
      
      # Inputs for showing boxplots by status or sex
      output$right_sidebar <- renderUI(
        controlbarMenu(
          id = "boxplot_sidebar",
          controlbarItem(
            "Options",
            icon = shiny::icon("gears"),
            radioButtons("var", 
                         "Boxplot Variable",
                         c("Status",
                           "Sex"))
          )
        )
      )
    } else {
      output$right_sidebar <- NULL
      # On homepage, hide both the right sidebar and the toggle button
      hide(selector = "body > div.wrapper > header > nav > div:nth-child(4) > ul")
      removeClass(selector = "body", class = "control-sidebar-open")
    }
  })
  callModule(home, "home_ns", pheno_avg = pheno_avg)
  callModule(boxplots, "boxplot_ns", box_var = reactive(input$var))
  callModule(volcano, "volcano_ns", compare = reactive(input$comparison))
}

shinyApp(ui, server)
