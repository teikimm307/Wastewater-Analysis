library(dplyr)
library(ggplot2)
library(readxl)
library(VennDiagram)
library(grid)
# reading in the data
feature_table <- read_excel("~/Downloads/Wastewater_FeatureTable5.xlsx")
mapfile <- read_excel("~/Downloads/Wastewater_Mapfile3.xlsx")

# beginning of SECTION 1

# function that will be used to append .mzXML to the end of column names
id1_function <- function(x, type, ext) {
  x <- mapfile$File.Name[mapfile$Sample_type == type]
  x <-paste0(x, ext)
}
# samples are divided into study samples and water blanks
study_sample_ids <- id1_function(study_sample_ids, "Study_Sample", ".mzXML")
water_blank_ids <- id1_function(water_blank_ids, "Water", ".mzXML")

# function to calculate the adjusted average intensity of water blanks
calculate_adjusted_average_intensity <- function(row, water_blank_ids) {
  water_blank_intensities <- row[water_blank_ids]
  average_intensity <- mean(water_blank_intensities, na.rm = TRUE)
  return(average_intensity)
}

# function to calculate the maximum intensity of study samples
calculate_max_intensity <- function(row, study_sample_ids) {
  study_sample_intensities <- row[study_sample_ids]
  max(study_sample_intensities, na.rm = TRUE)
}

# apply the average intensity function to each row
adjusted_average_intensities <- apply(feature_table[, water_blank_ids], 1, calculate_adjusted_average_intensity, water_blank_ids)
rows_with_zero_intensity <- adjusted_average_intensities == 0
new_table <- feature_table[rows_with_zero_intensity, ]
rows_with_zero_intensity_df <- new_table
# replace 0s with 1s for fold change calculations 
adjusted_average_intensities[rows_with_zero_intensity] <- 1
# apply the max intensity function to each row 
max_intensities <- apply(feature_table[, study_sample_ids], 1, calculate_max_intensity, study_sample_ids)

# calculate fold changes
fold_changes <- max_intensities / adjusted_average_intensities
feature_table$FoldChange <- fold_changes 
feature_table$adjusted_average_intensities <- adjusted_average_intensities
filtered_feature_table <- feature_table %>%
  # filter for rows with fold change values greater than or equal to 5
  filter(FoldChange>=5)

# create a modified filtered feature table that only contains rows with average intensity of 1 for working with sample sites
modified_filtered_feature_table <- filtered_feature_table %>%
  filter(adjusted_average_intensities == 1)
View(modified_filtered_feature_table)

# function for adding .mzXML to mapfile  
add_suffix <- function(x, column_name, ext) {
  x[[column_name]] <- paste0(x[[column_name]], ext)
  return(x)
}

# the modified mapfile with column names appended with .mzXML
mapfile2 <- add_suffix(mapfile, "File.Name", ".mzXML")

get_columns_for_site <- function(mapfile, site_number) {
  if(site_number < 10) {
    site_number <- paste0("0", site_number)
  }
  file_names <- mapfile$File.Name[mapfile$Sample_site == site_number]
  return(file_names)
}

#  find the number of chemicals that are found at each site from 1 to 28 
counters <- numeric(28)
for(site in 1:28) {
  site_columns <- get_columns_for_site(mapfile2, site)
  counters[site] <- sum(apply(modified_filtered_feature_table[site_columns], 1, function(row) {
    any(row != 0)
  }))
}

# table containing the site number and corresponding number of chemical compounds 
counters_df <- data.frame(Site = 1:28, Count = counters)

# bar graph of the site number and number of chemicals detected  
ggplot(counters_df, aes(x = Site, y = Count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Number of Chemicals Detected At Each Site", x = "Site Number", y = "Count") +
  theme(text = element_text(family = "Times New Roman"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30, by=5), limits = c(0,30))

# end of SECTION 1 

# beginning of SECTION 2

# function to filter out columns by sample time (A,B,C) and append ".mzXML" to columns 
id2_function <- function(x, type, ext) {
  x <- mapfile$File.Name[mapfile$Sample_time == type]
  x <-paste0(x, ext)
}
sample_time_A_ids <- id2_function(sample_time_A_ids, "A", ".mzXML")
sample_time_B_ids <- id2_function(sample_time_B_ids, "B", ".mzXML")
sample_time_C_ids <- id2_function(sample_time_C_ids, "C", ".mzXML")

# find the index of the starting column
start_col_index <- which(colnames(filtered_feature_table) == "R230705_CLU0002_C18neg_001.mzXML")

# create a new table with columns from the starting column onwards (easier to work with)
new_feature_table <- filtered_feature_table[, start_col_index:ncol(filtered_feature_table)]

# regression analysis for the three different sample times 
perform_regression <- function(row, ids) {
  # replace 0s with 1s and log-transform
  intensities <- log2(ifelse(row[ids] == 0, 1, row[ids]))
  
  # prepare data for regression
  data_for_regression <- data.frame(Intensity = intensities, SampleTime = 1:length(ids))
  
  # perform linear regression
  lm_model <- lm(Intensity ~ SampleTime, data = data_for_regression)
  summary_model <- summary(lm_model)
  
  # extract and return necessary values
  c(
    PValue = coef(summary(lm_model))["SampleTime", "Pr(>|t|)"],
    AdjustedR2 = summary_model$adj.r.squared,
    slope = coef(lm_model)["SampleTime"]  
  )
}

# apply the linear regression function to each row for each sample time
results_A <- t(apply(new_feature_table, 1, perform_regression, ids = sample_time_A_ids))
results_B <- t(apply(new_feature_table, 1, perform_regression, ids = sample_time_B_ids))
results_C <- t(apply(new_feature_table, 1, perform_regression, ids = sample_time_C_ids))

# convert them into data frames 
results_A_df <- as.data.frame(results_A)
colnames(results_A_df) <- c("PValue_A", "AdjustedR2_A", "Slope_A")
results_B_df <- as.data.frame(results_B)
colnames(results_B_df) <- c("PValue_B", "AdjustedR2_B", "Slope_B")
results_C_df <- as.data.frame(results_C)
colnames(results_C_df) <- c("PValue_C", "AdjustedR2_C", "Slope_C")

# the combined feature table contains information from the first 26 columns and the regression analysis for A,B, and C
mod_feature_table <- filtered_feature_table[,c(1:26,ncol(filtered_feature_table))]
combined_feature_table <- cbind(mod_feature_table, results_A_df, results_B_df, results_C_df)

# Manhattan plot function 
create_manhattan_plot <- function(data, retention_time_col_index, p_value_col_index, title) {
  # create a new data frame specifically for plotting
  plot_data <- data.frame(
    retention_time = data[, retention_time_col_index],
    p_values = -log10(data[, p_value_col_index])
  )
  
  # Create the Manhattan Plot using ggplot
  ggplot(plot_data, aes(x = retention_time, y = p_values)) +
    geom_point(aes(color = p_values > -log10(0.05))) +  # Plot points
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") + # Significance line at p-value = 0..05 
    ggtitle(title) +
    theme_minimal() +
    labs(x = "Retention Time", y = "-log10(P-Value)") +
    theme(text = element_text(family = "Times New Roman"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), plot.title = element_text(hjust = 0.5))  # Rotate x-axis labels for better readability
}

# Call the create_manhattan_plot() function 
plot_A <- create_manhattan_plot(combined_feature_table, 6, 28, "Manhattan Plot for Sample Time A")
plot_B <- create_manhattan_plot(combined_feature_table, 6, 31, "Manhattan Plot for Sample Time B")
plot_C <- create_manhattan_plot(combined_feature_table, 6, 34, "Manhattan Plot for Sample Time C")
#print(plot_A)
#print(plot_B)
#print(plot_C)

# end of SECTION 2

#beginning of SECTION 3

# function to create tables containing rows with p-value less than 0.05 
filter_by_p_value <- function(data, p_value_col_index) {
  filtered_data <- subset(data, data[p_value_col_index] < 0.05)
  if (p_value_col_index == 28) { 
    filtered_data <- filtered_data[, -c(31:36)]
  }
  else if (p_value_col_index == 31) { 
    filtered_data <- filtered_data[, -c(28:30, 34:36)]
  }
  else {
  filtered_data <- filtered_data[, -c(28:33)]
  }
  return(filtered_data)
}
plot_A_significant <- filter_by_p_value(combined_feature_table, 28)
plot_B_significant <- filter_by_p_value(combined_feature_table, 31)
plot_C_significant <- filter_by_p_value(combined_feature_table, 34)

# order and sort function that will be used to filter the 20 chemicals with lowest p-values 
order_and_sort <- function(data, col_name, num) {
  data1 <- data.frame(data)
  data1_mod1 <- data1[order(data1[[col_name]]),]
  data1_mod2 <- data1_mod1[1:num, ]
  return(data1_mod2)
}

# Extract unique chemical_IDs from each table
chemical_ids_A <- unique(plot_A_significant$chemical_ID)
chemical_ids_B <- unique(plot_B_significant$chemical_ID)
chemical_ids_C <- unique(plot_C_significant$chemical_ID)

# Find intersection between A and B, B and C, A and C
only_A = setdiff(chemical_ids_A, union(chemical_ids_B, chemical_ids_C))
only_B = setdiff(chemical_ids_B, union(chemical_ids_A, chemical_ids_C))
only_C = setdiff(chemical_ids_C, union(chemical_ids_A, chemical_ids_B))
intersection_AB_all = intersect(chemical_ids_A, chemical_ids_B)
intersection_BC_all = intersect(chemical_ids_B, chemical_ids_C)
intersection_AC_all = intersect(chemical_ids_A, chemical_ids_C)
common_intersection_ABC = intersect(intersect(intersection_AB_all, intersection_BC_all), intersection_AC_all)
intersection_AB = setdiff(intersection_AB_all, common_intersection_ABC)
intersection_BC = setdiff(intersection_BC_all, common_intersection_ABC)
intersection_AC = setdiff(intersection_AC_all, common_intersection_ABC)


# Convert to a dataframe or table for viewing or analysis
check_for_ID <- function(data, chemical_id, vec) {
  data1 <- data.frame(data)
  data_filtered <- data1[data1[[chemical_id]] %in% vec, ]
  return(data_filtered)
}
remove_duplicates <- function(data, col) {
  data1 <- data.frame(data)
  data1_unique <- data1[!duplicated(data1[[col]]), ]
  return(data1_unique)
}

# converting the data into data frames
plot_A_significant_unique <- remove_duplicates(plot_A_significant, "chemical_ID")
plot_B_significant_unique <- remove_duplicates(plot_B_significant, "chemical_ID")
plot_C_significant_unique <- remove_duplicates(plot_C_significant, "chemical_ID")
only_A_df <- check_for_ID(plot_A_significant_unique, "chemical_ID", only_A)
only_B_df <- check_for_ID(plot_B_significant_unique, "chemical_ID", only_B)
only_C_df <- check_for_ID(plot_C_significant_unique, "chemical_ID", only_C)
intersection_AB_df <- check_for_ID(plot_A_significant_unique, "chemical_ID", intersection_AB)
intersection_BC_df <- check_for_ID(plot_B_significant_unique, "chemical_ID", intersection_BC)
intersection_AC_df <- check_for_ID(plot_C_significant_unique, "chemical_ID", intersection_AC)
common_intersection_A_ABC_df <-  check_for_ID(plot_A_significant_unique, "chemical_ID", common_intersection_ABC)
common_intersection_B_ABC_df <-  check_for_ID(plot_B_significant_unique, "chemical_ID", common_intersection_ABC)
common_intersection_C_ABC_df <-  check_for_ID(plot_C_significant_unique, "chemical_ID", common_intersection_ABC)


# extracting the top 10 chemicals for each sample site based on p-value
only_A_df_2 <- order_and_sort(only_A_df, "PValue_A", 10)
only_A_df_top10 <- only_A_df_2 %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_A, AdjustedR2_A, Slope_A)
only_B_df_2 <- order_and_sort(only_B_df, "PValue_B", 10)
only_B_df_top10 <- only_B_df_2 %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_B, AdjustedR2_B, Slope_B)
only_C_df_2 <- order_and_sort(only_C_df, "PValue_C", 10)
only_C_df_top10 <- only_C_df_2 %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_C, AdjustedR2_C, Slope_C)

# extracting top 5 chemicals from each sample site for samples found in all sample sites. 
intersect_A_ABC <- order_and_sort(common_intersection_A_ABC_df, "PValue_A", 5)
intersect_A_ABC_top5 <- intersect_A_ABC %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_A, AdjustedR2_A, Slope_A)
intersect_B_ABC <- order_and_sort(common_intersection_B_ABC_df, "PValue_B", 5)
intersect_B_ABC_top5 <- intersect_B_ABC %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_B, AdjustedR2_B, Slope_B)
intersect_C_ABC <- order_and_sort(common_intersection_C_ABC_df, "PValue_C", 5)
intersect_C_ABC_top5 <- intersect_C_ABC %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue_C, AdjustedR2_C, Slope_C)

# reformatting the column headers for data presentation 
intersect_A_ABC_top5 <- intersect_A_ABC %>%
  rename(PValue = PValue_A, AdjustedR2 = AdjustedR2_A, Slope = Slope_A) %>%
  mutate(Sample_Time = 'A') %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue, AdjustedR2, Slope, Sample_Time)
intersect_B_ABC_top5 <- intersect_B_ABC %>%
  rename(PValue = PValue_B, AdjustedR2 = AdjustedR2_B, Slope = Slope_B) %>%
  mutate(Sample_Time = 'B') %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue, AdjustedR2, Slope, Sample_Time)
intersect_C_ABC_top5 <- intersect_C_ABC %>%
  rename(PValue = PValue_C, AdjustedR2 = AdjustedR2_C, Slope = Slope_C) %>%
  mutate(Sample_Time = 'C') %>%
  select(chemical_ID, Name, Formula, Adduct, mz_1, time_2, PValue, AdjustedR2, Slope, Sample_Time)

combined_df <- bind_rows(intersect_A_ABC_top5, intersect_B_ABC_top5, intersect_C_ABC_top5)

# Create a list of these sets
list_of_ids <- list(A = chemical_ids_A, B = chemical_ids_B, C = chemical_ids_C)

# Generate Venn Diagram
venn.plot <- venn.diagram(
  x = list_of_ids,
  category.names = c("Total Number of Features For Sample Time A", "Total Number of Features For Sample Time B", "Total Number of Features For Sample Time C"),
  fill = c("green", "purple", "orange"),
  cex = 2.5,
  cat.cex = 1.25,
  cat.dist = 0.1,
  cat.pos = c(330, 30, 180),
  margin = 0.1,
  filename = NULL,
  output = TRUE
)
# Display the plot
grid.newpage()
grid.draw(venn.plot)

write.csv(counters_df, "counters_df.csv", row.names = FALSE)
write.csv(only_A_df_top10, "only_A_df_top10(off).csv", row.names = FALSE)
write.csv(only_B_df_top10, "only_B_df_top10(off).csv", row.names = FALSE)
write.csv(only_C_df_top10, "only_C_df_top10(off).csv", row.names = FALSE)
write.csv(combined_df, "combined_df(off).csv", row.names = FALSE)





