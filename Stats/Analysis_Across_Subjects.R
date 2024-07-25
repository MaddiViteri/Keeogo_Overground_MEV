# Code to run statistical tests for Keeogo analysis
# Change change the path for when it is a different slope
# in "Section 2: Load excel files"
# Change

# Section 1: Load and set packages ----
# Un-comment the install lines if you do not have the packages
# downloaded already

#install.packages("pacman")

# Then load the package by using either of the following:
# require(pacman)  # Gives a confirmation message.
# pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
#                ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
#                stringr, tidyr) 
# install.packages("readxl")
# install.packages("purrr")
# install.packages("car")
# install.packages("dplyr")
# install.packages("openxlsx")


# Or, by using "pacman::p_load" you can use the p_load
# function from pacman without actually loading pacman.

library(pacman)
library(dplyr)
library(readxl)
library(purrr)
library(ggplot2)
library(car)
library(openxlsx)

# Section 2: Load excel files ----

# Set folder path where excel sheets contain
folder_path <- "C:/Users/MII36/University of Pittsburgh/Keeogo Project - Documents/Data Processing/Keeogo Data Analysis/Excel_Files/Uphill"

# Obtain files and file names
file_list <- list.files(folder_path, pattern = "*.xlsx", full.names = TRUE)
file_names <- list.files(folder_path, pattern = "*.xlsx", full.names = FALSE)

# Initialize an empty list to store the data
all_data <- list()

# Loop through each file
for (file in file_list) {
  # Get sheet names
  sheet_names <- excel_sheets(file)
  
  # Initialize a list to store the data for each sheet in the current file
  file_data <- list()
  
  # Loop through each sheet and read the data
  for (sheet in sheet_names) {
    data <- read_excel(file, sheet = sheet)
    file_data[[sheet]] <- data
  }
  
  # Store the file data in the main list
  all_data[[file]] <- file_data
}

# Section 3: Create variables for more/less comparisons ----

# Initialize more/less severe data
more_severe_compare_data <- list()
less_severe_compare_data <- list()
conditions <- c("control", "keeogo")

# Initialize the nested lists
for (sheet_index in 1:length(sheet_names)){
  sheet <- sheet_names[sheet_index]
  # Initialize for each metric
  more_severe_compare_data[[sheet]] <- list()
  less_severe_compare_data[[sheet]] <- list()
  
  # Initialize for splitting with condition
  for (condition in conditions){
    more_severe_compare_data[[sheet]][[condition]] <- matrix(NA, nrow = length(file_list),ncol = 1)
    less_severe_compare_data[[sheet]][[condition]] <- matrix(NA, nrow = length(file_list),ncol = 1)
  }
}

# Aggregate data for severe and condition
for (file_index in 1:length(file_list)){
  file <- file_list[file_index]
  for (sheet_index in 1:length(sheet_names)){
    sheet <- sheet_names[sheet_index]
    more_severe_control <- all_data[[file]][[sheet]][["Mean"]][1]
    more_severe_keeogo <- all_data[[file]][[sheet]][["Mean"]][2]
    more_severe_compare_data[[sheet]][[conditions[1]]][file_index] <-
      more_severe_control
    more_severe_compare_data[[sheet]][[conditions[2]]][file_index] <-
      more_severe_keeogo
    
    less_severe_control <- all_data[[file]][[sheet]][["Mean"]][3]
    less_severe_keeogo <- all_data[[file]][[sheet]][["Mean"]][4]
    less_severe_compare_data[[sheet]][[conditions[1]]][file_index] <- 
      less_severe_control
    less_severe_compare_data[[sheet]][[conditions[2]]][file_index]<-
      less_severe_keeogo
  }
}

# Section 4: More Severe Statistical tests----

# Initialize data frames
shapiro_results_df <- data.frame(
  set_name = character(),
  W = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
t_tests_result_df <- data.frame(
  metric = character(),
  t = numeric(),
  p_value = numeric(),
  conf_int_low = numeric(),
  conf_int_high = numeric(),
  mean_diff = numeric(),
  std_error = numeric(),
  stringsAsFactors = FALSE
)
wilcoxon_result_df <- data.frame(
  metric = character(),
  v = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
levene_results_df <- data.frame(
  metric = character(),
  F = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Initialize severity
severity <- list(
  more_severe_compare_data = more_severe_compare_data,
  less_severe_compare_data = less_severe_compare_data
  )

# Test for Normality and Homogeneity of Variance for more Severe
for (sheet in sheet_names){
  
  # Set the two conditions and find differences
  condition1 <- more_severe_compare_data[[sheet]][conditions[1]]
  condition2 <- more_severe_compare_data[[sheet]][conditions[2]]
  differences <- condition1$control - condition2$keeogo
  
  # Test for normality
  shapiro_test <- shapiro.test(differences)
  
  # Store normality test results
  shapiro_results_df <- rbind(shapiro_results_df,data.frame(
    set_name = sheet,
    W = shapiro_test$statistic,
    p_value = shapiro_test$p.value
  ))
  
  # Combine the data into a single vector
  combined_data <- c(condition1$control, condition2$keeogo)
  
  # Create a corresponding group vector
  group <- factor(rep(c("control", "keeogo"), each = length(condition1$control)))
  
  # Perform Levene's test
  levene_test <- leveneTest(combined_data ~ group)
  
  # Store the results in the list
  levene_results_df <- rbind(levene_results_df,data.frame(
    metric = sheet,
    F = levene_test$`F value`[1],
    p_value = levene_test$`Pr(>F)`[1]
  ))
  
  # Determine whether to use paired t-test or Wilcoxon-signed rank test
  if (shapiro_test$p.value > 0.05 && levene_test$`Pr(>F)`[1]>0.05 && length(condition1$control)>10){
    print("The differences are normally distributed. Use paired t-test.")
    
    # Conduct t-test
    t_test <- t.test(more_severe_compare_data[[sheet]][["control"]],
                     more_severe_compare_data[[sheet]][["keeogo"]], 
                     alternative = "two.sided",paired = TRUE)
    
    # Print if t-test results are statistically significant
    if (t_test$p.value < 0.05){
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are statistically significant")
    } else {
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are NOT statistically significant")
    }
    # Store t-test results
    t_tests_result_df <- rbind(t_tests_result_df,data.frame(
      metric = sheet,
      t = t_test$statistic,
      p_value = t_test$p.value,
      conf_int_low = t_test$conf.int[1],
      conf_int_high = t_test$conf.int[2],
      mean_diff = t_test$estimate,
      std_error = t_test$stderr
    ))
    
  } else {
    print("The differences are not normally distributed. Use Wilcoxon signed-rank test.")
    
    # Perform Wilcoxon signed-rank test
    wilcoxon_test <- wilcox.test(more_severe_compare_data[[sheet]][["control"]],
                                 more_severe_compare_data[[sheet]][["keeogo"]], 
                                 alternative = "two.sided",paired = TRUE)
    
    # Print if statistically significant
    if (wilcoxon_test$p.value < 0.05){
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are statistically significant")
    } else {
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are NOT statistically significant")
    }
    
    # Store data into data frame
    wilcoxon_result_df <- rbind(wilcoxon_result_df,data.frame(
      metric = sheet,
      v = wilcoxon_test$statistic,
      p_value = wilcoxon_test$p.value
    ))
  }
}

# Section 5: Less Severe Statistical tests----

# Initialize data frames
shapiro_results_less_severe_df <- data.frame(
  set_name = character(),
  W = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
t_tests_result_less_severe_df <- data.frame(
  metric = character(),
  t = numeric(),
  p_value = numeric(),
  conf_int_low = numeric(),
  conf_int_high = numeric(),
  mean_diff = numeric(),
  std_error = numeric(),
  stringsAsFactors = FALSE
)
wilcoxon_result_less_severe_df <- data.frame(
  metric = character(),
  v = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
levene_results_less_severe_df <- data.frame(
  metric = character(),
  F = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Test for Normality and Homogeneity of Variance for more Severe
for (sheet in sheet_names){
  if (sheet == "speed"){
    next
  }
  # Set the two conditions and find differences
  condition1 <- less_severe_compare_data[[sheet]][conditions[1]]
  condition2 <- less_severe_compare_data[[sheet]][conditions[2]]
  differences <- condition1$control - condition2$keeogo
  
  # Test for normality
  shapiro_test_less_severe <- shapiro.test(differences)
  
  # Store normality test results
  shapiro_results_less_severe_df <- rbind(shapiro_results_less_severe_df,data.frame(
    set_name = sheet,
    W = shapiro_test_less_severe$statistic,
    p_value = shapiro_test_less_severe$p.value
  ))
  
  # Combine the data into a single vector
  combined_data <- c(condition1$control, condition2$keeogo)
  
  # Create a corresponding group vector
  group <- factor(rep(c("control", "keeogo"), each = length(condition1$control)))
  
  # Perform Levene's test
  levene_test_less_severe <- leveneTest(combined_data ~ group)
  
  # Store the results in the list
  levene_results_less_severe_df <- rbind(levene_results_less_severe_df,data.frame(
    metric = sheet,
    F = levene_test_less_severe$`F value`[1],
    p_value = levene_test_less_severe$`Pr(>F)`[1]
  ))
  
  # Determine whether to use paired t-test or Wilcoxon-signed rank test
  if (shapiro_test_less_severe$p.value > 0.05 && levene_test_less_severe$`Pr(>F)`[1]>0.05 && length(condition1$control)>10){
    print("The differences are normally distributed. Use paired t-test.")
    
    # Conduct t-test
    t_test <- t.test(less_severe_compare_data[[sheet]][["control"]],
                     less_severe_compare_data[[sheet]][["keeogo"]], 
                     alternative = "two.sided",paired = TRUE)
    
    # Print if t-test results are statistically significant
    if (t_test$p.value < 0.05){
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are statistically significant")
    } else {
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are NOT statistically significant")
    }
    # Store t-test results
    t_tests_result_less_severe_df <- rbind(t_tests_result_less_severe_df,data.frame(
      metric = sheet,
      t = t_test$statistic,
      p_value = t_test$p.value,
      conf_int_low = t_test$conf.int[1],
      conf_int_high = t_test$conf.int[2],
      mean_diff = t_test$estimate,
      std_error = t_test$stderr
    ))
    
  } else {
    print("The differences are not normally distributed. Use Wilcoxon signed-rank test.")
    
    # Perform Wilcoxon signed-rank test
    wilcoxon_test_less_severe <- wilcox.test(less_severe_compare_data[[sheet]][["control"]],
                                 less_severe_compare_data[[sheet]][["keeogo"]], 
                                 alternative = "two.sided",paired = TRUE)
    
    # Print if statistically significant
    if (wilcoxon_test_less_severe$p.value < 0.05){
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are statistically significant")
    } else {
      paste("The differences between the control and Keeogo conditions for",
            sheet, "are NOT statistically significant")
    }
    
    # Store data into data frame
    wilcoxon_result_less_severe_df <- rbind(wilcoxon_result_less_severe_df,data.frame(
      metric = sheet,
      v = wilcoxon_test_less_severe$statistic,
      p_value = wilcoxon_test_less_severe$p.value
    ))
  }
}
#
# Section 6: Save data to excel file ----
# Create a new workbook
# Function to write nested list data frames to the workbook
write_nested_list <- function(wb,new_sheet_name,nested_list) {
  current_col <- 1
  addWorksheet(wb, new_sheet_name)
  for (outer_name in names(nested_list)) {
    for (inner_name in names(nested_list[[outer_name]])) {
      # sheet_name <- paste(outer_name, inner_name, sep = "_")
      
      writeData(wb, new_sheet_name, paste(outer_name, inner_name, sep = "_"), startCol = current_col, startRow = 1)
      writeData(wb, new_sheet_name, nested_list[[outer_name]][[inner_name]],
                startCol = current_col,startRow = 2, colNames =TRUE)

      # Move to the next column group
      current_col <- current_col + 1
    }
  }
}

wb <-createWorkbook()

# Add data frames to the workbook
addWorksheet(wb, "More Severe T-Test Results")
writeData(wb, "More Severe T-Test Results",t_tests_result_df)

addWorksheet(wb,"More Severe Wilcoxon Results")
writeData(wb,"More Severe Wilcoxon Results", wilcoxon_result_df)

write_nested_list(wb,"More Severe Data",more_severe_compare_data)

addWorksheet(wb, "Less Severe T-Test Results")
writeData(wb, "Less Severe T-Test Results",t_tests_result_less_severe_df)

addWorksheet(wb,"Less Severe Wilcoxon Results")
writeData(wb,"Less Severe Wilcoxon Results", wilcoxon_result_less_severe_df)

write_nested_list(wb,"Less Severe Data",less_severe_compare_data)

output_directory <- "C:/Users/MII36/University of Pittsburgh/Keeogo Project - Documents/Data Processing/Keeogo Data Analysis/Stats"
if (grepl("neutral",file_names[1])){
  append_name <- "neutral"
} else if (grepl("uphill",file_names[1])){
  append_name <- "uphill"
} else if (grepl("downhill",file_names[1])){
  append_name <- "downhill"
}
output_file <- file.path(output_directory, 
                         gsub(" ", "",
                              paste("test_results_", append_name,
                                    ".xlsx",collapse = NULL)))
# Save the workbook
saveWorkbook(wb, file = output_file, overwrite = TRUE)
