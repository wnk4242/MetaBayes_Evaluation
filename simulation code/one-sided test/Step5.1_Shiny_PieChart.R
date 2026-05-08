# There are two Shiny apps in this script, based on Step4.0_Data_Visualization_PieChart2 and Step4.0_Data_Visualization_PieChart1
# The first Shiny don't include anecdotal evidence so you can compare FEMA with MABF methods
# The second Shiny include anecdotal evidence so you only compare MABF methods with themselves. Use the second version so you can 
# select BF cutoffs.
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

#######Not including anecdotal evidence ################
#######Prepare for original p-cutoff = 0.01##############
# This Shiny is based on Step4.0_Data_Visualization_PieChart2
# This Shiny creates pie plots for MABF methods and FE meta-analysis
# Import data sets
# Set the path to the directory containing RDS files
folder_path0.01 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c0/" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files0.01 <- list.files(folder_path0.01, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files0.01) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
remove(rates_FEMA_0.2null_c0);remove(rates_FEMA_0.5null_c0);
# Import FEMA rates when pcutoff_r = 0.05 is used (c1)
rates_FEMA_0.2null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c1/rates_FEMA_0.2null_c1.RDS")
rates_FEMA_0.5null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c1/rates_FEMA_0.5null_c1.RDS")

# Add a method column containing the method name in each dataset as identifier
rates_BFbMA_0.2null_c0 <- rates_BFbMA_0.2null_c0 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.2null_c0 <- rates_EUBF_0.2null_c0 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.2null_c0 <- rates_FEMABF_0.2null_c0 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.2null_c0 <- rates_iBF_0.2null_c0 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.2null_c1 <- rates_FEMA_0.2null_c1 %>% 
  mutate(method = rep("FEMA", 81)) %>% 
  relocate(method, .before = everything())
rates_BFbMA_0.5null_c0 <- rates_BFbMA_0.5null_c0 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.5null_c0 <- rates_EUBF_0.5null_c0 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.5null_c0 <- rates_FEMABF_0.5null_c0 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.5null_c0 <- rates_iBF_0.5null_c0 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.5null_c1 <- rates_FEMA_0.5null_c1 %>% 
  mutate(method = rep("FEMA", 81)) %>% 
  relocate(method, .before = everything())

# Combine MABF and FEMA datasets
rates_MABFEMA_0.2null0.01 <- rbind(rates_FEMA_0.2null_c1,rates_BFbMA_0.2null_c0,rates_EUBF_0.2null_c0,rates_FEMABF_0.2null_c0,rates_iBF_0.2null_c0)
rates_MABFEMA_0.5null0.01 <- rbind(rates_FEMA_0.5null_c1,rates_BFbMA_0.5null_c0,rates_EUBF_0.5null_c0,rates_FEMABF_0.5null_c0,rates_iBF_0.5null_c0)

rates_MABFEMA_0.2null0.01 <- rates_MABFEMA_0.2null0.01 %>%
  mutate(p_cutoff = 0.01) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rates_MABFEMA_0.5null0.01 <- rates_MABFEMA_0.5null0.01 %>%
  mutate(p_cutoff = 0.01) %>% 
  mutate(true_es = 0.5) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rm(list = setdiff(ls(), c("rates_MABFEMA_0.2null0.01", "rates_MABFEMA_0.5null0.01")))


#######Prepare for original p-cutoff = 0.05##############
# Import data sets
# Set the path to the directory containing RDS files
folder_path0.05 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c0/" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files0.05 <- list.files(folder_path0.05, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files0.05) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
remove(rates_FEMA_0.2null_c0);remove(rates_FEMA_0.5null_c0);
# Import FEMA rates when pcutoff_r = 0.05 is used (c1)
rates_FEMA_0.2null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.2null_c1.RDS")
rates_FEMA_0.5null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.5null_c1.RDS")

# Add a method column containing the method name in each dataset as identifier
rates_BFbMA_0.2null_c0 <- rates_BFbMA_0.2null_c0 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.2null_c0 <- rates_EUBF_0.2null_c0 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.2null_c0 <- rates_FEMABF_0.2null_c0 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.2null_c0 <- rates_iBF_0.2null_c0 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.2null_c1 <- rates_FEMA_0.2null_c1 %>% 
  mutate(method = rep("FEMA", 81)) %>% 
  relocate(method, .before = everything())
rates_BFbMA_0.5null_c0 <- rates_BFbMA_0.5null_c0 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_EUBF_0.5null_c0 <- rates_EUBF_0.5null_c0 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMABF_0.5null_c0 <- rates_FEMABF_0.5null_c0 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_iBF_0.5null_c0 <- rates_iBF_0.5null_c0 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) %>% 
  select(-ADRNE,-ADRTE,-ADR,-AD,-ADTE,-ADNE)
rates_FEMA_0.5null_c1 <- rates_FEMA_0.5null_c1 %>% 
  mutate(method = rep("FEMA", 81)) %>% 
  relocate(method, .before = everything())

# Combine MABF and FEMA datasets
rates_MABFEMA_0.2null0.05 <- rbind(rates_FEMA_0.2null_c1,rates_BFbMA_0.2null_c0,rates_EUBF_0.2null_c0,rates_FEMABF_0.2null_c0,rates_iBF_0.2null_c0)
rates_MABFEMA_0.5null0.05 <- rbind(rates_FEMA_0.5null_c1,rates_BFbMA_0.5null_c0,rates_EUBF_0.5null_c0,rates_FEMABF_0.5null_c0,rates_iBF_0.5null_c0)

rates_MABFEMA_0.2null0.05 <- rates_MABFEMA_0.2null0.05 %>%
  mutate(p_cutoff = 0.05) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rates_MABFEMA_0.5null0.05 <- rates_MABFEMA_0.5null0.05 %>%
  mutate(p_cutoff = 0.05) %>% 
  mutate(true_es = 0.5) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

# Remove data that won't be used
rm(list = setdiff(ls(), c(
  "rates_MABFEMA_0.2null0.01",
  "rates_MABFEMA_0.5null0.01",
  "rates_MABFEMA_0.2null0.05",
  "rates_MABFEMA_0.5null0.05"
)))

# Combine final dataset
rates_all <- bind_rows(
  rates_MABFEMA_0.2null0.01,
  rates_MABFEMA_0.5null0.01,
  rates_MABFEMA_0.2null0.05,
  rates_MABFEMA_0.5null0.05
)


######## Last step before Shiny ##########
# Turn wide-format into long-format for pie chart visualization
prepare_pie_data <- function(data, true_es_value, use_TSn = "TS1") {
  metric_cols <- if (use_TSn == "TS1") {
    c("TS1", "FF1", "TF1", "FS1")
  } else {
    c("TS2", "FF2", "TF2", "FS2")
  }
  
  data %>%
    filter(true_es == true_es_value) %>%
    mutate(
      orig.n = factor(orig.n, levels = c(20, 50, 200)),
      rep.number = factor(rep.number, levels = c(2, 5, 10),
                          labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
      rep.n = factor(rep.n, levels = c(40, 100, 400),
                     labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
      bias.level = factor(censorFunc, levels = c("low", "medium", "high")),
      p_cutoff = factor(p_cutoff, levels = c(0.01, 0.05)),
      true_es = factor(true_es, levels = c(0.2, 0.5)),
      method = factor(method, levels = c("FEMA", "BFbMA", "EUBF", "FEMABF", "iBF"))
    ) %>%
    relocate(bias.level, .after = orig.n) %>%
    select(method, rep.n, rep.number, orig.n, bias.level, p_cutoff, true_es, all_of(metric_cols)) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = recode(metric,
                            TS1 = "true success", FF1 = "false failure", TF1 = "true failure", FS1 = "false success",
                            TS2 = "true success", FF2 = "false failure", TF2 = "true failure", FS2 = "false success")
    ) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))
}


# plot_data_null <- prepare_pie_data(rates_all, true_es_value = 0.2, use_TSn = "TS2")
# plot_data_0.2  <- prepare_pie_data(rates_all, true_es_value = 0.2, use_TSn = "TS1")
# plot_data_0.5  <- prepare_pie_data(rates_all, true_es_value = 0.5, use_TSn = "TS1")



#################Shiny App#####################

# UI
ui <- fluidPage(
  titlePanel("MABF Pie Chart Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null (TS2)" = "null", "0.2", "0.5"), selected = "null"),
      selectInput("p_cutoff", "Original p-value Cutoff",
                  choices = c("All", "0.01", "0.05"), selected = "All"),
      selectInput("orig_n", "Original Sample Size",
                  choices = c("All", "20", "50", "200"), selected = "All"),
      selectInput("bias_level", "Bias Level",
                  choices = c("All", "low", "medium", "high"), selected = "All")
    ),
    mainPanel(
      plotOutput("piePlot", height = "800px")
    )
  )
)

# Server
server <- function(input, output) {
  
  filtered_data <- reactive({
    use_TSn <- if (input$dataset_choice == "null") "TS2" else "TS1"
    true_es_num <- if (input$dataset_choice == "null") 0.2 else as.numeric(input$dataset_choice)
    
    metric_cols <- if (use_TSn == "TS1") {
      c("TS1", "FF1", "TF1", "FS1")
    } else {
      c("TS2", "FF2", "TF2", "FS2")
    }
    
    df <- rates_all %>%
      # Always filter true_es based on dataset_choice
      filter(true_es == true_es_num) %>%
      # Conditionally filter p_cutoff
      { if (input$p_cutoff != "All") filter(., p_cutoff == as.numeric(input$p_cutoff)) else . } %>%
      # Conditionally filter orig.n
      { if (input$orig_n != "All") filter(., orig.n == as.numeric(input$orig_n)) else . } %>%
      # Conditionally filter bias level
      { if (input$bias_level != "All") filter(., censorFunc == input$bias_level) else . } %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10),
                            labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
        rep.n = factor(rep.n, levels = c(40, 100, 400),
                       labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
        method = factor(method, levels = c("FEMA", "BFbMA", "EUBF", "FEMABF", "iBF"))
      ) %>%
      select(method, rep.n, rep.number, all_of(metric_cols)) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(
        metric_label = dplyr::recode(metric,
                                     TS1 = "true success", FF1 = "false failure", TF1 = "true failure", FS1 = "false success",
                                     TS2 = "true success", FF2 = "false failure", TF2 = "true failure", FS2 = "false success"
        ),
        metric_label = factor(metric_label, levels = c(
          "true success",
          "false failure",
          "true failure",
          "false success"
        ))
      ) %>%
      group_by(method, rep.n, rep.number, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value))
    
    
    return(df)
  })
  
  output$piePlot <- renderPlot({
    plot_data <- filtered_data()
    
    ggplot(plot_data, aes(x = "", y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = c("true success" = "#56B4E9",
                                   "false failure" = "#F0E442",
                                   "true failure" = "#0072B2",
                                   "false success" = "#009E73")) +
      labs(
        x = "", y = "Proportion", fill = "Evidence Category",
        title = paste0("Pie Charts of Evidence Categories\n(data = ", input$dataset_choice,
                       ", p_cutoff = ", input$p_cutoff,
                       ", orig.n = ", input$orig_n,
                       ", bias.level = ", input$bias_level, ")")
      ) +
      facet_grid(method ~ rep.n + rep.number, scales = "free",
                 labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  })
}

# Run the app
shinyApp(ui = ui, server = server)






#########################################################
####### Including anecdotal evidence ####################
#######Prepare for original p-cutoff = 0.01##############
# The following script is based on Step4.0_Data_Visualization_PieChart1
# Replace all c1 to c2 to increase BFcutoff = 1/3 and 3 to BFcutoff = 1/10 and 10
# Import data sets
# Set the path to the directory containing RDS files
folder_path0.01 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.01, EScutoff_o=0/c1/" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files0.01 <- list.files(folder_path0.01, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files0.01) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
remove(rates_FEMA_0.2null_c1);remove(rates_FEMA_0.5null_c1);

# Add a method column containing the method name in each dataset as identifier
rates_BFbMA_0.2null_c1 <- rates_BFbMA_0.2null_c1 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) 

rates_EUBF_0.2null_c1 <- rates_EUBF_0.2null_c1 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_FEMABF_0.2null_c1 <- rates_FEMABF_0.2null_c1 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_iBF_0.2null_c1 <- rates_iBF_0.2null_c1 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_BFbMA_0.5null_c1 <- rates_BFbMA_0.5null_c1 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) 

rates_EUBF_0.5null_c1 <- rates_EUBF_0.5null_c1 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_FEMABF_0.5null_c1 <- rates_FEMABF_0.5null_c1 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_iBF_0.5null_c1 <- rates_iBF_0.5null_c1 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) 



# Combine MABF and FEMA datasets
rates_MABF_0.2null0.01 <- rbind(rates_BFbMA_0.2null_c1,rates_EUBF_0.2null_c1,rates_FEMABF_0.2null_c1,rates_iBF_0.2null_c1)
rates_MABF_0.5null0.01 <- rbind(rates_BFbMA_0.5null_c1,rates_EUBF_0.5null_c1,rates_FEMABF_0.5null_c1,rates_iBF_0.5null_c1)

rates_MABF_0.2null0.01 <- rates_MABF_0.2null0.01 %>%
  mutate(p_cutoff = 0.01) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rates_MABF_0.5null0.01 <- rates_MABF_0.5null0.01 %>%
  mutate(p_cutoff = 0.01) %>% 
  mutate(true_es = 0.5) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rm(list = setdiff(ls(), c("rates_MABF_0.2null0.01", "rates_MABF_0.5null0.01")))


#######Prepare for original p-cutoff = 0.05##############
# Import data sets
# Set the path to the directory containing RDS files
folder_path0.05 <- "./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/" #These dataset contain original study p and ES for calculating TSR, FSR, etc.
# List all RDS files in the directory
rds_files0.05 <- list.files(folder_path0.05, pattern = "\\.RDS$", full.names = TRUE)
# Read in all RDS files into the workspace
for (file_path in rds_files0.05) {
  # Extract the base name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Create a variable with the name of the file and assign the data from readRDS
  assign(file_name, readRDS(file_path), envir = .GlobalEnv)
}
remove(rates_FEMA_0.2null_c1);remove(rates_FEMA_0.5null_c1);
# Import FEMA rates when pcutoff_r = 0.05 is used (c1)
rates_FEMA_0.2null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.2null_c1.RDS")
rates_FEMA_0.5null_c1 <- readRDS("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1/rates_FEMA_0.5null_c1.RDS")

# Add a method column containing the method name in each dataset as identifier
rates_BFbMA_0.2null_c1 <- rates_BFbMA_0.2null_c1 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) 

rates_EUBF_0.2null_c1 <- rates_EUBF_0.2null_c1 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_FEMABF_0.2null_c1 <- rates_FEMABF_0.2null_c1 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_iBF_0.2null_c1 <- rates_iBF_0.2null_c1 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_BFbMA_0.5null_c1 <- rates_BFbMA_0.5null_c1 %>% 
  mutate(method = rep("BFbMA", 81)) %>% 
  relocate(method, .before = everything()) 

rates_EUBF_0.5null_c1 <- rates_EUBF_0.5null_c1 %>% 
  mutate(method = rep("EUBF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_FEMABF_0.5null_c1 <- rates_FEMABF_0.5null_c1 %>% 
  mutate(method = rep("FEMABF", 81)) %>% 
  relocate(method, .before = everything()) 

rates_iBF_0.5null_c1 <- rates_iBF_0.5null_c1 %>% 
  mutate(method = rep("iBF", 81)) %>% 
  relocate(method, .before = everything()) 



# Combine MABF datasets
rates_MABF_0.2null0.05 <- rbind(rates_BFbMA_0.2null_c1,rates_EUBF_0.2null_c1,rates_FEMABF_0.2null_c1,rates_iBF_0.2null_c1)
rates_MABF_0.5null0.05 <- rbind(rates_BFbMA_0.5null_c1,rates_EUBF_0.5null_c1,rates_FEMABF_0.5null_c1,rates_iBF_0.5null_c1)

rates_MABF_0.2null0.05 <- rates_MABF_0.2null0.05 %>%
  mutate(p_cutoff = 0.05) %>% 
  mutate(true_es = 0.2) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

rates_MABF_0.5null0.05 <- rates_MABF_0.5null0.05 %>%
  mutate(p_cutoff = 0.05) %>% 
  mutate(true_es = 0.5) %>% 
  relocate(p_cutoff, .before = method) %>% 
  relocate(true_es, .before = p_cutoff)

# Remove data that won't be used
rm(list = setdiff(ls(), c(
  "rates_MABF_0.2null0.01",
  "rates_MABF_0.5null0.01",
  "rates_MABF_0.2null0.05",
  "rates_MABF_0.5null0.05"
)))

# Combine final dataset
rates_all <- bind_rows(
  rates_MABF_0.2null0.01,
  rates_MABF_0.5null0.01,
  rates_MABF_0.2null0.05,
  rates_MABF_0.5null0.05
)


######## Last step before Shiny ##########
# Turn wide-format into long-format for pie chart visualization
prepare_pie_data <- function(data, true_es_value, use_TSn = "TS1") {
  metric_cols <- if (use_TSn == "TS1") {
    c("ADTE", "TS1", "FF1", "TF1", "FS1")
  } else {
    c("ADNE", "TS2", "FF2", "TF2", "FS2")
  }
  
  data %>%
    filter(true_es == true_es_value) %>%
    mutate(
      orig.n = factor(orig.n, levels = c(20, 50, 200)),
      rep.number = factor(rep.number, levels = c(2, 5, 10),
                          labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
      rep.n = factor(rep.n, levels = c(40, 100, 400),
                     labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
      bias.level = factor(censorFunc, levels = c("low", "medium", "high")),
      p_cutoff = factor(p_cutoff, levels = c(0.01, 0.05)),
      true_es = factor(true_es, levels = c(0.2, 0.5)),
      method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
    ) %>%
    relocate(bias.level, .after = orig.n) %>%
    select(method, rep.n, rep.number, orig.n, bias.level, p_cutoff, true_es, all_of(metric_cols)) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = recode(metric,
                            ADTE = "anecdotal evidence", TS1 = "true success", FF1 = "false failure", TF1 = "true failure", FS1 = "false success", 
                            ADNE = "anecdotal evidence", TS2 = "true success", FF2 = "false failure", TF2 = "true failure", FS2 = "false success")
    ) %>%
    group_by(method, rep.n, rep.number) %>%
    mutate(proportion = value / sum(value))
}


# plot_data_null <- prepare_pie_data(rates_all, true_es_value = 0.2, use_TSn = "TS2")
# plot_data_0.2  <- prepare_pie_data(rates_all, true_es_value = 0.2, use_TSn = "TS1")
# plot_data_0.5  <- prepare_pie_data(rates_all, true_es_value = 0.5, use_TSn = "TS1")



#################Shiny App#####################
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)

# UI
ui <- fluidPage(
  titlePanel("MABF Pie Chart Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null (TS2)" = "null", "0.2", "0.5"), selected = "null"),
      selectInput("p_cutoff", "Original p-value Cutoff",
                  choices = c("All", "0.01", "0.05"), selected = "All"),
      selectInput("orig_n", "Original Sample Size",
                  choices = c("All", "20", "50", "200"), selected = "All"),
      selectInput("bias_level", "Bias Level",
                  choices = c("All", "low", "medium", "high"), selected = "All")
    ),
    mainPanel(
      plotOutput("piePlot", height = "800px")
    )
  )
)

# Server
server <- function(input, output) {
  
  filtered_data <- reactive({
    use_TSn <- if (input$dataset_choice == "null") "TS2" else "TS1"
    true_es_num <- if (input$dataset_choice == "null") 0.2 else as.numeric(input$dataset_choice)
    
    metric_cols <- if (use_TSn == "TS1") {
      c("ADTE", "TS1", "FF1", "TF1", "FS1")
    } else {
      c("ADNE", "TS2", "FF2", "TF2", "FS2")
    }
    
    df <- rates_all %>%
      # Always filter true_es based on dataset_choice
      filter(true_es == true_es_num) %>%
      # Conditionally filter p_cutoff
      { if (input$p_cutoff != "All") filter(., p_cutoff == as.numeric(input$p_cutoff)) else . } %>%
      # Conditionally filter orig.n
      { if (input$orig_n != "All") filter(., orig.n == as.numeric(input$orig_n)) else . } %>%
      # Conditionally filter bias level
      { if (input$bias_level != "All") filter(., censorFunc == input$bias_level) else . } %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10),
                            labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
        rep.n = factor(rep.n, levels = c(40, 100, 400),
                       labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      ) %>%
      select(method, rep.n, rep.number, all_of(metric_cols)) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(
        metric_label = dplyr::recode(metric,
                                     ADTE = "anecdotal evidence", TS1 = "true success", FF1 = "false failure", TF1 = "true failure", FS1 = "false success",
                                     ADNE = "anecdotal evidence", TS2 = "true success", FF2 = "false failure", TF2 = "true failure", FS2 = "false success"
        ),
        metric_label = factor(metric_label, levels = c(
          "anecdotal evidence",
          "true success",
          "false failure",
          "true failure",
          "false success"
          
        ))
      ) %>%
      group_by(method, rep.n, rep.number, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value))
    
    
    return(df)
  })
  
  output$piePlot <- renderPlot({
    plot_data <- filtered_data()
    
    ggplot(plot_data, aes(x = "", y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = c("anecdotal evidence" = "#E69F00",
                                   "true success" = "#56B4E9",
                                   "false failure" = "#F0E442",
                                   "true failure" = "#0072B2",
                                   "false success" = "#009E73"
                                   
      )) +
      labs(
        x = "", y = "Proportion", fill = "Evidence Category",
        title = paste0("Pie Charts of Evidence Categories\n(data = ", input$dataset_choice,
                       ", p_cutoff = ", input$p_cutoff,
                       ", orig.n = ", input$orig_n,
                       ", bias.level = ", input$bias_level, ")")
      ) +
      facet_grid(method ~ rep.n + rep.number, scales = "free",
                 labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  })
}

# Run the app
shinyApp(ui = ui, server = server)



###################
# The following version allows you to select BF cutoff (3 or 10)
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)

# === Data Preparation Function ===
load_and_prepare_data <- function(cutoff_folder = "c1") {
  cutoff_label <- ifelse(cutoff_folder == "c1", 3, 10)
  
  load_group <- function(folder_path) {
    rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
    lapply(rds_files, function(file_path) {
      data <- readRDS(file_path)
      name <- tools::file_path_sans_ext(basename(file_path))
      assign(name, data, envir = .GlobalEnv)
    })
  }
  
  get_data_combined <- function(pcutoff, true_es) {
    folder <- paste0("./MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=",
                     pcutoff, ", EScutoff_o=0/", cutoff_folder, "/")
    
    load_group(folder)
    
    method_list <- list(
      BFbMA = get(paste0("rates_BFbMA_", true_es, "null_", cutoff_folder)),
      EUBF = get(paste0("rates_EUBF_", true_es, "null_", cutoff_folder)),
      FEMABF = get(paste0("rates_FEMABF_", true_es, "null_", cutoff_folder)),
      iBF = get(paste0("rates_iBF_", true_es, "null_", cutoff_folder))
    )
    
    combined <- bind_rows(lapply(names(method_list), function(method) {
      method_list[[method]] %>%
        mutate(method = method) %>%
        relocate(method, .before = everything())
    }))
    
    combined %>%
      mutate(p_cutoff = as.numeric(pcutoff),
             true_es = as.numeric(true_es)) %>%
      relocate(true_es, p_cutoff)
  }
  
  df_all <- bind_rows(
    get_data_combined("0.01", "0.2"),
    get_data_combined("0.01", "0.5"),
    get_data_combined("0.05", "0.2"),
    get_data_combined("0.05", "0.5")
  )
  
  return(df_all)
}

# === UI ===
ui <- fluidPage(
  titlePanel("MABF Pie Chart Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("bf_cutoff", "Bayes Factor Cutoff",
                  choices = c("3" = "c1", "10" = "c2"), selected = "c1"),
      selectInput("dataset_choice", "Data Scenario",
                  choices = c("Null (TS2)" = "null", "0.2", "0.5"), selected = "null"),
      selectInput("p_cutoff", "Original p-value Cutoff",
                  choices = c("All", "0.01", "0.05"), selected = "All"),
      selectInput("orig_n", "Original Sample Size",
                  choices = c("All", "20", "50", "200"), selected = "All"),
      selectInput("bias_level", "Bias Level",
                  choices = c("All", "low", "medium", "high"), selected = "All")
    ),
    mainPanel(
      plotOutput("piePlot", height = "800px")
    )
  )
)

# === Server ===
server <- function(input, output) {
  # Reactive dataset
  rates_all <- reactive({
    load_and_prepare_data(input$bf_cutoff)
  })
  
  # Reactive pie chart data
  filtered_data <- reactive({
    use_TSn <- if (input$dataset_choice == "null") "TS2" else "TS1"
    true_es_num <- if (input$dataset_choice == "null") 0.2 else as.numeric(input$dataset_choice)
    
    metric_cols <- if (use_TSn == "TS1") {
      c("ADTE", "TS1", "FF1", "TF1", "FS1")
    } else {
      c("ADNE", "TS2", "FF2", "TF2", "FS2")
    }
    
    df <- rates_all() %>%
      filter(true_es == true_es_num) %>%
      { if (input$p_cutoff != "All") filter(., p_cutoff == as.numeric(input$p_cutoff)) else . } %>%
      { if (input$orig_n != "All") filter(., orig.n == as.numeric(input$orig_n)) else . } %>%
      { if (input$bias_level != "All") filter(., censorFunc == input$bias_level) else . } %>%
      mutate(
        rep.number = factor(rep.number, levels = c(2, 5, 10),
                            labels = c("N[rep] == 2", "N[rep] == 5", "N[rep] == 10")),
        rep.n = factor(rep.n, levels = c(40, 100, 400),
                       labels = c("n[rep] == 40", "n[rep] == 100", "n[rep] == 400")),
        method = factor(method, levels = c("BFbMA", "EUBF", "FEMABF", "iBF"))
      ) %>%
      select(method, rep.n, rep.number, all_of(metric_cols)) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(
        metric_label = dplyr::recode(metric,
                                     ADTE = "anecdotal evidence", TS1 = "true success", FF1 = "false failure", TF1 = "true failure", FS1 = "false success",
                                     ADNE = "anecdotal evidence", TS2 = "true success", FF2 = "false failure", TF2 = "true failure", FS2 = "false success"),
        metric_label = factor(metric_label, levels = c("anecdotal evidence", "true success", "false failure", "true failure", "false success"))
      ) %>%
      group_by(method, rep.n, rep.number, metric_label) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(method, rep.n, rep.number) %>%
      mutate(proportion = value / sum(value))
    
    return(df)
  })
  
  output$piePlot <- renderPlot({
    plot_data <- filtered_data()
    
    ggplot(plot_data, aes(x = "", y = proportion, fill = metric_label)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = c(
        "anecdotal evidence" = "#E69F00",
        "true success" = "#56B4E9",
        "false failure" = "#F0E442",
        "true failure" = "#0072B2",
        "false success" = "#009E73"
      )) +
      labs(
        x = "", y = "Proportion", fill = "Evidence Category",
        title = paste0("Pie Charts of Evidence Categories\n(data = ", input$dataset_choice,
                       ", BF cutoff = ", ifelse(input$bf_cutoff == "c1", 3, 10),
                       ", p_cutoff = ", input$p_cutoff,
                       ", orig.n = ", input$orig_n,
                       ", bias.level = ", input$bias_level, ")")
      ) +
      facet_grid(method ~ rep.n + rep.number, scales = "free",
                 labeller = labeller(rep.number = label_parsed, rep.n = label_parsed)) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  })
}

# === Run the App ===
shinyApp(ui = ui, server = server)


