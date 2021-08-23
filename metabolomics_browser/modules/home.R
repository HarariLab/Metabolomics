home_UI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "home",
    h1("ONTIME Metabolomics Browser"),
    htmlOutput(ns("description"))
  )
}

home <- function(input, output, session, pheno_avg) {
  
  abbreviations <- c("CO", "AD", "ADAD", "TREM2") #, "FTD", "OT", "Presymptomatic")
  full <- c("Control", "Sporadic AD Case", "Autosomal Dominant AD", "Carrier of TREM2 Variant") #,
            #"Frontotemporal Dementia", "Other Neurodegenerative Disease", "Presymptomatic")
  n <- c( table(pheno_avg$Final_CACO)["CO"],  table(pheno_avg$Final_CACO)["CA"], 
          table(pheno_avg$Final_CACO)["ADAD"],  table(pheno_avg$Final_CACO)["TREM2"]) #,
          # table(pheno_avg$Final_CACO)["FTD"],  table(pheno_avg$Final_CACO)["OT"],
          # table(pheno_avg$Final_CACO)["Presymptomatic"])
  pct_female <- c(66.67, 61.64, 32.00, 52.38) #, 25.00, 64.44, 60.00)
  age <- c("88.00 ± 8.99", "83.72 ± 8.81", "54.2 ± 13.87", "83.52 ± 7.48") #),
           # "71.44 ± 10.71", "83.36 ± 11.95", "86.87 ± 9.71")
  
  abbr_key <- 
    data.frame(
      Abbreviation = abbreviations,
      `Full Name` = full,
      n = n,
      Age = age,
      `% Female` = pct_female,
      check.names = FALSE
    )
  
  output$description <- renderText(
    paste(
      "Here you can browse metabolomics data and analyses generated
      from the DIAN and Knight-ADRC cohorts. The \"Reading Distributions\" tab
      allows viewing of detailed information on metabolites, as well as their
      reading distributions between disease statuses and between sexes.  The 
      \"Volcano Plots\" tab shows volcano plots for linear models of differential
      metabolites between the case, control, ADAD, and TREM2 groups, as well as 
      information on the various metabolites shown in the plots.  A summary of 
      sample phenotypes can be found below.</br>",
      kable(abbr_key, row.names = FALSE) %>%
        kable_styling("striped", full_width = FALSE, position = "left")
      )
  )
}


