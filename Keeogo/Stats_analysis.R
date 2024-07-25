# install.packages("pacman")

# Then load the package by using either of the following:
# require(pacman)  # Gives a confirmation message.
library(pacman)  # No message.

# Or, by using "pacman::p_load" you can use the p_load
# function from pacman without actually loading pacman.
# These are packages I load every time.
# pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
#                ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
#                stringr, tidyr) 
# install.packages("readxl")
# install.packages("purrr")
install.packages("dplyr")
library(dplyr)
library(readxl)
library(purrr)
library(ggplot2)


file_path <- "C:/Users/MII36/University of Pittsburgh/Keeogo Project - Documents/Data Processing/Keeogo Data Analysis/Excel_Files"
# file_path <- "C:/Users/MII36/University of Pittsburgh/Keeogo Project - Documents/Data Processing/Local/test_data1.xlsx"
sheet_names <- excel_sheets(file_path)
sheet_names <- trimws(sheet_names)
all_data <- set_names(map(sheet_names, ~ read_excel(file_path, sheet = .x)),sheet_names)
# ROM <- all_data[["ROM"]]
# head(ROM)
test_data <- all_data[["grade"]]
CCR_Loading_control_r <- all_data[["CC R_Loading_control_r"]]
CCR_Loading_keeogo_r <- all_data[["CCR_Loading_keeogo_r"]]
CCR_Early <- all_data[["CCR_Early"]]
CCR_Mid <- all_data[["CCR_Mid"]]
CCR_Late <- all_data[["CCR_Late"]]
CCI_Loading <- all_data[["CCI_Loading"]]
CCI_Early <- all_data[["CCI_Early"]]
CCI_Mid <- all_data[["CCI_Mid"]]
CCI_Late <- all_data[["CCI_Late"]]
CCI_Combined <- all_data[["CCI_Combined"]]
## EMG
# Early_Flexion <- all_data[["Early_Flexion"]]
# Mid_Flexion <- all_data[["Mid_Flexion"]]
# Late_Flexion <- all_data[["Late_Flexion"]]
# Early_Adduction <- all_data[["Early_Adduction"]]
# Mid_Adduction <- all_data[["Mid_Adduction"]]
# Late_Adduction <- all_data[["Late_Adduction"]]
## ROM
# boxplot(ROM_right ~ Condition_r, data=ROM)
# boxplot(ROM_left ~ Condition_l, data=ROM)
# t.test(ROM_right ~ Condition_r, data=ROM)
# t.test(ROM_left ~ Condition_l, data=ROM)
t.test(CCR_Loading_control_r, CCR_Loading_keeogo_r)
#CCR Loading
boxplot(CCR_zero_r_loading ~ Condition_r, data=CCR_Loading)
boxplot(CCR_zero_l_loading ~ Condition_l, data=CCR_Loading)
t.test(CCR_zero_r_loading ~ Condition_r, data=CCR_Loading)
t.test(CCR_zero_l_loading ~ Condition_l, data=CCR_Loading)
#CCR Early
boxplot(CCR_zero_r_early ~ Condition_r, data=CCR_Early)
boxplot(CCR_zero_l_early ~ Condition_l, data=CCR_Early)
t.test(CCR_zero_r_early ~ Condition_r, data=CCR_Early)
t.test(CCR_zero_l_early ~ Condition_l, data=CCR_Early)
#CCR Mid
boxplot(CCR_zero_r_mid ~ Condition_r, data=CCR_Mid)
boxplot(CCR_zero_l_mid ~ Condition_l, data=CCR_Mid)
t.test(CCR_zero_r_mid ~ Condition_r, data=CCR_Mid)
t.test(CCR_zero_l_mid ~ Condition_l, data=CCR_Mid)
#CCR Late
boxplot(CCR_zero_r_late ~ Condition_r, data=CCR_Late)
boxplot(CCR_zero_l_late ~ Condition_l, data=CCR_Late)
t.test(CCR_zero_r_late ~ Condition_r, data=CCR_Late)
t.test(CCR_zero_l_late ~ Condition_l, data=CCR_Late)

#CCI Loading
boxplot(CCI_r_loading ~ Condition_r, data=CCI_Loading)
boxplot(CCI_l_loading ~ Condition_l, data=CCI_Loading)
t.test(CCI_r_loading ~ Condition_r, data=CCI_Loading)
t.test(CCI_l_loading ~ Condition_l, data=CCI_Loading)
#CCI Early
boxplot(CCI_r_early ~ Condition_r, data=CCI_Early)
boxplot(CCI_l_early ~ Condition_l, data=CCI_Early)
t.test(CCI_r_early ~ Condition_r, data=CCI_Early)
t.test(CCI_l_early ~ Condition_l, data=CCI_Early)
#CCI Mid
boxplot(CCI_r_mid ~ Condition_r, data=CCI_Mid)
boxplot(CCI_l_mid ~ Condition_l, data=CCI_Mid)
t.test(CCI_r_mid ~ Condition_r, data=CCI_Mid)
t.test(CCI_l_mid ~ Condition_l, data=CCI_Mid)
#CCI Late
boxplot(CCI_r_late ~ Condition_r, data=CCI_Late)
boxplot(CCI_l_late ~ Condition_l, data=CCI_Late)
t.test(CCI_r_late ~ Condition_r, data=CCI_Late)
t.test(CCI_l_late ~ Condition_l, data=CCI_Late)

summary_data <- CCI_Combined%>%
  group_by(Phase_r,Condition_r)%>%
  summarise(
    Mean = mean(CCI_r),
    SD = sd(CCI_r)
  )
custom_colors <- c("blue","red")
ggplot(summary_data, aes(x = Phase_r,y=Mean,fill=Condition_r)) +
  geom_bar(stat = "identity", position = "dodge", width=0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.7), width = 0.2) +
  labs(x = "Phase", y = "CCI Full") +
  ggtitle("CCI Full") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for better readability
  scale_fill_manual(values = custom_colors)

counts <- table(CCI_Combined$CCI_r,CCI_Combined$Condition_r)
barplot(counts,
        main = "CCI Full",
        xlab="Phase",col=c("blue","red"),
        legend=rownames(counts), beside=TRUE)
barplot()
## 
#Early Flexion
# boxplot(Flexion_r ~ Condition_r, data=Early_Flexion)
# boxplot(Flexion_l ~ Condition_l, data=Early_Flexion)
# t.test(Flexion_r ~ Condition_r, data=Early_Flexion)
# t.test(Flexion_l ~ Condition_l, data=Early_Flexion)
# #Mid Flexion
# boxplot(flexion_r ~ Condition_r, data=Mid_Flexion)
# boxplot(flexion_l ~ Condition_l, data=Mid_Flexion)
# t.test(flexion_r ~ Condition_r, data=Mid_Flexion)
# t.test(flexion_l ~ Condition_l, data=Mid_Flexion)
# #Late Flexion
# boxplot(flexion_r ~ Condition_r, data=Late_Flexion)
# boxplot(flexion_l ~ Condition_l, data=Late_Flexion)
# t.test(flexion_r ~ Condition_r, data=Late_Flexion)
# t.test(flexion_l ~ Condition_l, data=Late_Flexion)
# #
# # Early Adduction
# boxplot(adduction_r ~ Condition_r, data=Early_Adduction)
# boxplot(adduction_l ~ Condition_l, data=Early_Adduction)
# t.test(adduction_r ~ Condition_r, data=Early_Adduction)
# t.test(adduction_l ~ Condition_l, data=Early_Adduction)
# # Mid Adduction
# boxplot(adduction_r ~ Condition_r, data=Mid_Adduction)
# boxplot(adduction_l ~ Condition_l, data=Mid_Adduction)
# t.test(adduction_r ~ Condition_r, data=Mid_Adduction)
# t.test(adduction_l ~ Condition_l, data=Mid_Adduction)
# # Late Adduction
# boxplot(adduction_r ~ Condition_r, data=Late_Adduction)
# boxplot(adduction_l ~ Condition_l, data=Late_Adduction)
# t.test(adduction_r ~ Condition_r, data=Late_Adduction)
# t.test(adduction_l ~ Condition_l, data=Late_Adduction)

# tryCatch({
#   library(readxl)
#   library(purrr)
#   
#   # file_path <- "C:/Path/To/Your/Excel/File.xlsx"
#   sheet_names <- excel_sheets(file_path)
#   
#   all_data <- set_names(map(sheet_names, ~ read_excel(file_path, sheet = .x)), sheet_names)
# }, error = function(e) {
#   print(paste("Error:", e))
# })
