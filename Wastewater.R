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

#  find the number of chemicals that were found at each site from 1 to 28 
counters <- numeric(28)
site_data_frames <- list()

for(site in 1:28) {
  # Get the file names for the site using the existing function
  site_files <- get_columns_for_site(mapfile2, site)
  
  # Determine the column names in modified_filtered_feature_table that match the site files
  site_cols <- which(colnames(modified_filtered_feature_table) %in% site_files)
  
  # Apply a function across the rows to check for any non-zero values in the site's columns
  # and sum these to count the number of non-zero rows for the site
  counters[site] <- sum(apply(modified_filtered_feature_table[site_cols], 1, function(row) any(row != 0)))
  
  # Subset the data for rows where any column for the site has a non-zero value
  site_data <- modified_filtered_feature_table[apply(modified_filtered_feature_table[site_cols], 1, function(row) any(row != 0)), ]
  
  # Store the subsetted data frame in the list with a name
  site_data_frames[[paste0("found_site_", site)]] <- site_data
}

#BEGINNING OF SECTION INTERSECTION

# Find the common chemical IDs between the two sites
# BEGINNING OF Site 27
common_rows_list <- list()

# Loop through sites 1 to 4 (extend this range as needed)
for(site in 1:5) {
  # Calculate the intersecting chemical IDs between the current site and site 27
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_27"]][["chemical_ID"]])
  
  # Subset the site 27 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_27"]][site_data_frames[["found_site_27"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_27")]] <- common_rows
}

for(site in 7:14) {
  # Calculate the intersecting chemical IDs between the current site and site 27
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_27"]][["chemical_ID"]])
  
  # Subset the site 27 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_27"]][site_data_frames[["found_site_27"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_27")]] <- common_rows
}

for(site in 25:26) {
  # Calculate the intersecting chemical IDs between the current site and site 27
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_27"]][["chemical_ID"]])
  
  # Subset the site 27 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_27"]][site_data_frames[["found_site_27"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_27")]] <- common_rows
}

for(site in 28) {
  # Calculate the intersecting chemical IDs between the current site and site 27
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_27"]][["chemical_ID"]])
  
  # Subset the site 27 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_27"]][site_data_frames[["found_site_27"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_27")]] <- common_rows
}
# END OF Site 27

# BEGINNING OF Site 15

for(site in 16) {
  # Calculate the intersecting chemical IDs between the current site and site 16
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_15"]][["chemical_ID"]])
  
  # Subset the site 15 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_15"]][site_data_frames[["found_site_15"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_15")]] <- common_rows
}
#END OF Site 15


#BEGINNING OF Site 21
for(site in 18) {
  # Calculate the intersecting chemical IDs between the current site and site 21
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_21"]][["chemical_ID"]])
  
  # Subset the site 21 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_21"]][site_data_frames[["found_site_21"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_21")]] <- common_rows
}

for(site in 20) {
  # Calculate the intersecting chemical IDs between the current site and site 21
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_21"]][["chemical_ID"]])
  
  # Subset the site 21 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_21"]][site_data_frames[["found_site_21"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_21")]] <- common_rows
}

for(site in 22:24) {
  # Calculate the intersecting chemical IDs between the current site and site 21
  common_chemical_ids <- intersect(site_data_frames[[paste0("found_site_", site)]][["chemical_ID"]],
                                   site_data_frames[["found_site_21"]][["chemical_ID"]])
  
  # Subset the site 21 data frame for rows with the common chemical IDs
  common_rows <- site_data_frames[["found_site_21"]][site_data_frames[["found_site_21"]][["chemical_ID"]] %in% common_chemical_ids, ]
  
  # Store the result in the list using the site number as the name
  common_rows_list[[paste0("site_", site, "_and_21")]] <- common_rows
}

# END OF SECTION INTERSECTION



# table containing the site number and corresponding number of chemical compounds 
counters_df <- data.frame(Site = 1:28, Count = counters)

# bar graph of the site number and number of chemicals detected  
ggplot(counters_df, aes(x = Site, y = Count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(title = "Number of Chemicals Detected At Each Site", x = "Site Number", y = "Count") +
  theme(text = element_text(family = "Times New Roman"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30, by=5), limits = c(0,30))


#modified bar graph that groups bars together based on what site they converge to
site_order <- c(
  as.character(1:5), as.character(7:14), as.character(25:26), "28", "27", # Group 1 with 27 last
  "16", "15",                                                             # Group 2 with 15 last
  "20", as.character(22:24), "21", "6"                                    # Group 3 with 21 last
)

# Set the levels of the Site factor according to the custom order
counters_df$Site <- factor(counters_df$Site, levels = site_order)

counters_df$Group <- ifelse(counters_df$Site %in% c(1:5, 7:14, 25:26, 28, 27), "Group 1",
                            ifelse(counters_df$Site %in% c(15,16), "Group 2",
                                   ifelse(counters_df$Site %in% c(20, 22:24, 21), "Group 3", NA)))

# Define a new column for color (factor variable for differentiating specific bars)
counters_df$Color <- ifelse(counters_df$Site == 27, "Special",
                            ifelse(counters_df$Site == 15, "Special",
                                   ifelse(counters_df$Site == 21, "Special", "Standard")))

# Plotting the bar graph
ggplot(counters_df, aes(x = Group, y = Count, fill = Color, group = Site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.45), width = 0.37) +
  theme_minimal() +
  labs(title = "Number of Chemicals Detected At Sites Compared to Their Corresponding WQTC", x = "Group", y = "Count") +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Standard" = "blue", "Special" = "orange")) +
  scale_x_discrete(limits = c("Group 1", "Group 2", "Group 3"))

# end of SECTION 1 


# beginning of SECTION 2

# function to filter out columns by sample time (A,B,C) and append ".mzXML" to columns 
# Initial setup: Extracting sample time IDs
id2_function <- function(x, type, ext) {
  x <- mapfile$File.Name[mapfile$Sample_time == type]
  x <- paste0(x, ext)
}
sample_time_A_ids <- id2_function(sample_time_A_ids, "A", ".mzXML")
sample_time_B_ids <- id2_function(sample_time_B_ids, "B", ".mzXML")
sample_time_C_ids <- id2_function(sample_time_C_ids, "C", ".mzXML")

# Find the index of the starting column
start_col_index <- which(colnames(filtered_feature_table) == "R230705_CLU0002_C18neg_001.mzXML")

# Create a new table with columns from the starting column onwards
new_feature_table <- filtered_feature_table[, start_col_index:ncol(filtered_feature_table)]

# Regression analysis across sample times
perform_regression <- function(row, sample_time_ids) {
  # Create a numeric vector representing the sample times
  # Assuming 'A' corresponds to 1, 'B' to 2, and 'C' to 3
  sample_times <- rep(1:3, times = sapply(sample_time_ids, length))
  intensities <- unlist(lapply(sample_time_ids, function(ids) row[ids]))
  intensities <- log2(ifelse(intensities == 0, 1, intensities))
  
  # Prepare and perform the regression
  data_for_regression <- data.frame(Intensity = intensities, SampleTime = sample_times)
  lm_model <- lm(Intensity ~ SampleTime, data = data_for_regression)
  
  # Extract and return the p-value for the linear trend
  summary_lm_model <- summary(lm_model)
  p_value_linear_trend <- coef(summary_lm_model)[2, "Pr(>|t|)"]
  
  return(p_value_linear_trend)
}

# Apply the regression to each row (feature)
p_values_linear <- apply(new_feature_table, 1, perform_regression, sample_time_ids = list(sample_time_A_ids, sample_time_B_ids, sample_time_C_ids))
results_df <- data.frame(P_Value_Linear_Trend = p_values_linear)
results_df$mz_1 <- filtered_feature_table$mz_1 #adding a column for mass to charge ratio
results_df$time_1 <- filtered_feature_table$time_1 #adding a column for retention time

# Manhattan plot functions
create_manhattan_plot_feature_id <- function(results_df) { #Manhattan plot in terms of Feature ID's 
  plot_data <- data.frame(
    FeatureID = 1:nrow(results_df),
    p_values = -log10(results_df$P_Value_Linear_Trend)
  )
  
  ggplot(plot_data, aes(x = FeatureID)) +
    geom_point(aes(y = p_values, color = p_values > -log10(0.05)), alpha = 1) +  
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) +
    ggtitle("Manhattan Plot for Changes Across Feature ID's") +
    xlab("Feature ID") +
    ylab("-log(p-value)") +
    theme(text = element_text(family = "Times New Roman"), plot.title = element_text(hjust=0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

create_manhattan_plot_mz <- function(results_df) { #Manhattan plot in terms of mass-to-charge
  plot_data <- data.frame(
    mz = results_df$mz_1,
    p_values = -log10(results_df$P_Value_Linear_Trend)
  )
  
  ggplot(plot_data, aes(x = mz, y = p_values)) +
    geom_point(aes(color = p_values > -log10(0.05)), alpha = 1, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) +
    ggtitle("Manhattan Plot for Changes Across Mass-To-Charge Ratios") +
    xlab("Mass-To-Charge Ratios") +
    ylab("-log(p-value)") +
    theme(text = element_text(family = "Times New Roman"), plot.title = element_text(hjust=0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none")
}

create_manhattan_plot_rt <- function(results_df) { #Manhattan plot in terms of retention time
  plot_data <- data.frame(
    rt = results_df$time_1,
    p_values = -log10(results_df$P_Value_Linear_Trend)
  )
  
  ggplot(plot_data, aes(x = rt, y = p_values)) +
    geom_point(aes(color = p_values > -log10(0.05)), alpha = 1, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) +
    ggtitle("Manhattan Plot for Changes Across Retention Times") +
    xlab("Retention Times") +
    ylab("-log(p-value)") +
    theme(text = element_text(family = "Times New Roman"), plot.title = element_text(hjust=0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none")
}

# Generate and display the Manhattan plots
create_manhattan_plot_feature_id(results_df)
create_manhattan_plot_mz(results_df)
create_manhattan_plot_rt(results_df)

# End of SECTION 2


#Beginning of SECTION 3
#PART 1: SCATTERPLOTS FOR INCOME, POPULATION, AND AREA
incomes <- c(0,34490, 36812, 37059, 40168, 31050, 27752, 27517, 65791, 56672, 27517, 74500, 101140, 86478, 81994, 108021, 55433, 76796, 72401, 55346, 61923, 45794, 53857, 61081, 28054, 52751, 20000)
population <- c(0, 8258, 3459, 1610, 3587, 10949, 9073, 23751, 145346, 8838, 23751, 95603, 11444, 40824, 5781, 37193, 78206, 60885, 24969, 309998, 45148, 37972, 23135, 309184, 41777, 350766, 74)
area <- c(0.04, 4, 1, 1, 3, 5, 5, 12, 112, 3, 12, 80, 12, 67, 11, 88, 55, 80, 23, 332, 37, 28, 21, 242, 20, 280, 3)

#creating the data table from the paper
mod_counters_df <- counters_df[-6, ]
mod_counters_df$Income <- incomes
mod_counters_df$Population <- population
mod_counters_df$Area <- area

mod_perform_regression <- function(data, independent, dependent) {
  formula <- as.formula(paste(dependent, "~", independent))
  lm_model <- lm(formula, data = data)
  summary_model <- summary(lm_model)
  
  # Extracting the coefficients table
  coefficients_table <- summary_model$coefficients
  
  # Adding R-squared and Adjusted R-squared
  summary_model$R_squared <- summary_model$r.squared
  summary_model$Adjusted_R_squared <- summary_model$adj.r.squared
  
  print(summary_model)
  
  # Renaming the 'Pr(>|t|)' column for clarity
  colnames(coefficients_table)[4] <- "PValue"
  
  conf_intervals <- confint(lm_model, level = 0.95)
  print(conf_intervals)
  
  # Return the summary table with all the necessary statistics
  return(coefficients_table)
}

plot_regression <- function(data, independent, dependent, model) {
  # Create a scatter plot of the independent variable vs. the dependent variable
  plot(data[[independent]], data[[dependent]], main = paste("Regression of", dependent, "on", independent),
       xlab = independent, ylab = dependent, pch = 19)
  
  # Add the regression line
  abline(model, col = "blue")
}

income_results <- mod_perform_regression(mod_counters_df, "Income", "Count")
population_results <- mod_perform_regression(mod_counters_df, "Population", "Count")
area_results <- mod_perform_regression(mod_counters_df, "Area", "Count")

income_model <- lm(Count ~ Income, data = mod_counters_df)
plot_regression(mod_counters_df, "Income", "Count", income_model)
population_model <- lm(Count ~ Population, data = mod_counters_df)
plot_regression(mod_counters_df, "Population", "Count", population_model)
area_model <- lm(Count ~ Area, data = mod_counters_df)
plot_regression(mod_counters_df, "Area", "Count", area_model)

#PART 2: BOXPLOTS

#Working with nested sewersheds vs. treatment plants
count_sewer <- c(counters_df[1:5, "Count"], counters_df[7:14, "Count"], counters_df[16, "Count"], counters_df[18, "Count"], counters_df[20, "Count"], counters_df[22:26, "Count"], counters_df[28, "Count"])
count_plant <- c(counters_df[15, "Count"], counters_df[17, "Count"], counters_df[19, "Count"], counters_df[21, "Count"], counters_df[27, "Count"])
data_sample_sites <- data.frame(Group = rep(c("Sewershed", "WQTC"), c(length(count_sewer), length(count_plant))), Count = c(count_sewer, count_plant))
ggplot(data_sample_sites, aes(x = Group, y = Count)) + 
  geom_boxplot() + xlab("Group") + ylab("Count") + 
  geom_jitter(width = 0, alpha = 0.5, color = "blue") +
  ggtitle("Chemicals Found in Sewersheds and WQTC") +
  theme(text = element_text(family = "Times New Roman"), plot.title = element_text(hjust=0.5))

t_test_result <- t.test(Count ~ Group, data = data_sample_sites)

count_combined_overflow_yes <- c(counters_df[1:5, "Count"], counters_df[7:14, "Count"], counters_df[25:26, "Count"], counters_df[28, "Count"])
count_combined_overflow_no <- c(counters_df[16, "Count"], counters_df[18, "Count"], counters_df[20, "Count"], counters_df[22:24, "Count"])
data_combined_overflow <- data.frame(Group = rep(c("Combined Overflow", "No Combined Overflow"), c(length(count_combined_overflow_yes), length(count_combined_overflow_no))), Count = c(count_combined_overflow_yes, count_combined_overflow_no))
ggplot(data_combined_overflow, aes(x = Group, y = Count)) + 
  geom_boxplot() + xlab("Group") + ylab("Count") + 
  geom_jitter(width = 0, alpha = 0.5, color = "blue") +
  ggtitle("Chemicals Found in Sewersheds With/Without Combined Sewer Overflow") +
  theme(text = element_text(family = "Times New Roman"), plot.title = element_text(hjust=0.5))

t_test_result_1 <- t.test(Count ~ Group, data = data_combined_overflow)

print(t_test_result)
print(t_test_result_1)

#End of SECTION 3


#Beginning of SECTION 4

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

# order and sort function that will be used to filter the 10 chemicals with lowest p-values 
order_and_sort <- function(data, col_name, num) {
  data1 <- data.frame(data)
  data1_mod1 <- data1[order(data1[[col_name]]),]
  data1_mod2 <- data1_mod1[1:num, ]
  return(data1_mod2)
}

# Extract chemical_IDs from each table
chemical_ids_A <- plot_A_significant$chemical_ID
chemical_ids_B <- plot_B_significant$chemical_ID
chemical_ids_C <- plot_C_significant$chemical_ID

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

# converting the data into data frames
plot_A_significant_unique <- plot_A_significant
plot_B_significant_unique <- plot_B_significant
plot_C_significant_unique <- plot_C_significant
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

#write.csv(counters_df, "counters_df.csv", row.names = FALSE)
#write.csv(only_A_df_top10, "only_A_df_top10(off).csv", row.names = FALSE)
#write.csv(only_B_df_top10, "only_B_df_top10(off2).csv", row.names = FALSE)
#write.csv(only_C_df_top10, "only_C_df_top10(off).csv", row.names = FALSE)
#write.csv(combined_df, "combined_df(off).csv", row.names = FALSE)




