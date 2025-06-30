getwd()
setwd("#Location where you have the supplementary data from Kader et al paper")

#======================================================================================================
#Load the counts data from supplementary data
#=======================================================================================================


#Check if directory exists
dir.exists("#Location where you have the supplementary data from Kader et al paper")
# Remove ALL objects, including hidden ones
rm(list = ls(all.names = TRUE))
list.files()

#Load a csv file
library(readr) 
library(readxl)
target_counts <- read_excel("Kader et al_SuppTable S8_Q3 norm_countmatrix.xlsx",
                            sheet = "TargetCountMatrix")

#head(target_counts)
#write_csv(target_counts, "TargetCountMatrix.csv")

#==============================================================================================================
#Filter only the Full ROI samples
#==============================================================================================================


#Remove all other columns except | Full ROI since our target samples are from Full ROI
library(dplyr)
FullROI_data <- target_counts %>%
  select(TargetName, contains("| Full ROI"))
head(FullROI_data)
#write_csv(FullROI_data, "FullROI_only.csv")


# Read data_clinical_sample.csv 
clinical_sample <- read_excel("your file path/ovary_geomx_gray_foundation_2024/data_clinical_sample.xlsx")
head(clinical_sample)
colnames(clinical_sample)


#Extract sample Id
#sample_ids <- clinical_sample$`#Sample Identifier`

#Replace " - " to " | " in the #Sample Identifier to avoid discrepencies
library(stringr)
clinical_sample <- clinical_sample %>%
  mutate(
    `#Sample Identifier` = str_replace_all(
      `#Sample Identifier`, 
      pattern = "-", 
      replacement = " | "
    )
  )

#Replace Full | ROI with Full ROI
clinical_sample <- clinical_sample %>%
  mutate(
    `#Sample Identifier` = str_replace(
      `#Sample Identifier`, 
      pattern = " \\| (ROI)$",  # More explicit pattern
      replacement = " \\1"       # Keeps ROI but removes pipe
    )
  )


#========================================================================================================
# Filter samples for three lesion types in Epithelial cells
#========================================================================================================

filtered_data <- clinical_sample %>%
  filter(
    #`Lesions Types` %in% c("Incidental Fimbriae", "Incidental STIC", "Incidental p53sig"),
    `Lesions Types` %in% c("Incidental Fimbriae", "Incidental STIC"),
    `ROI Types` == "Epithelial"
  ) %>%
  mutate(status = case_when(
    `Lesions Types` == "Incidental Fimbriae" ~ "normal",
    #`Lesions Types` %in% c("Incidental STIC", "Incidental p53sig") ~ "transformed"
    `Lesions Types` %in% c("Incidental STIC") ~ "transformed"
  ))

# Count normal vs. transformed samples
sample_counts <- filtered_data %>%
  count(status) %>%
  rename(Number_of_Samples = n)

# Print the counts (optional)
print(sample_counts)

# Get ALL sample IDs we want to keep (combined normal + transformed)
matched_columns <- filtered_data %>% 
  pull(`#Sample Identifier`) %>% 
  intersect(colnames(FullROI_data))

# Verify matches
print(paste("Total matching samples found:", length(matched_columns)))

# Select operation with safety check
FullROI_filtered <- if(length(matched_columns) > 0) {
  FullROI_data %>%
    select(
      TargetName,  
      all_of(matched_columns)  # Now using the properly defined variable
    )
} else {
  stop("No matching columns found between filtered_data and FullROI_data")
}
# Filter only rows with CCN2 data and add status information
library(tidyr)
ccn2_data <- FullROI_filtered %>% 
  filter(TargetName == "CCN2") %>%
  pivot_longer(
    cols = -TargetName,
    names_to = "SampleID",
    values_to = "Expression"
  ) %>%
  left_join(
    filtered_data %>% 
      select(`#Sample Identifier`, status) %>% 
      rename(SampleID = `#Sample Identifier`),
    by = "SampleID"
  )

# Write to CSV
#write_csv(ccn2_data, "ccn2_filtered.csv")

#======================================================================================================================
#Start with creating the plot
#======================================================================================================================

library(tidyverse)
# Corrected pipeline
ccn2_long <- FullROI_filtered %>% 
  filter(TargetName == "CCN2") %>%
  pivot_longer(
    cols = -TargetName,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  # Extract Patient ID from sample name
  mutate(PatientID = str_extract(Sample, "LSP\\d+")) %>%
  # Add status information
  left_join(
    filtered_data %>% 
      select(`#Sample Identifier`, status) %>% 
      rename(Sample = `#Sample Identifier`),
    by = "Sample"
  )


# Calculate t-test p-value (parametric)
t_test_result <- ccn2_long %>% 
  t.test(Expression ~ status, data = .) %>% 
  broom::tidy() %>% 
  mutate(p_label = paste0("p = ", format.pval(p.value, digits = 2)))

# Create the boxplot with modified theme
ccn2_long %>%
  mutate(status = factor(status, levels = c("normal", "transformed"))) %>% 
  ggplot(aes(x = status, y = Expression, fill = status)) +
  
  # Boxplot and points
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  
  # Add t-test p-value
  geom_text(
    data = t_test_result,
    aes(x = 1.5, y = max(ccn2_long$Expression) * 1.1, label = p_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  # Colors and labels
  scale_fill_manual(values = c("normal" = "#1b9e77", "transformed" = "#d95f02")) +
  labs(
    title = "CCN2 Expression: Normal vs Transformed (t-test)",
    x = NULL,
    y = "CCN2 Expression"
  ) +
  
  # Modified theme - removes gridlines but keeps axes
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),  # Removes major gridlines
    panel.grid.minor = element_blank(),  # Removes minor gridlines
    axis.line = element_line(color = "black")  # Ensures axis lines remain
  )

# Create and store the plot
ccn2_plot <- ccn2_long %>%
  mutate(status = factor(status, levels = c("normal", "transformed"))) %>% 
  ggplot(aes(x = status, y = Expression, fill = status)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  geom_text(
    data = t_test_result,
    aes(x = 1.5, y = max(ccn2_long$Expression) * 1.1, label = p_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(values = c("normal" = "#1b9e77", "transformed" = "#d95f02")) +
  #labs(title = "CCN2 Expression: Normal vs Transformed(STIC and p53)", x = NULL, y = "CCN2 Expression") +
  labs(title = "CCN2 Expression: Normal vs Transformed(STIC)", x = NULL, y = "CCN2 Expression") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# Now save it
ggsave("CCN2_expression_boxplot (STIC).png", 
#ggsave("CCN2_expression_boxplot (STIC and p53).png",
       plot = ccn2_plot,
       width = 8, height = 6, dpi = 300, bg = "white")








