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

#  find the chemicals and number of chemicals that were found at each site from 1 to 28 
counters <- numeric(28)
site_data_frames <- list()

for(site in 1:28) {
  # Get the file names for the site using the existing function
  site_files <- get_columns_for_site(mapfile2, site)
  
  # Determine the column names in filtered_feature_table that match the site files
  site_cols <- which(colnames(filtered_feature_table) %in% site_files) #FROM THIS POINT ONWARDS, ALL MENTIONS OF FILTERED FEATURE TABLE CAN BE CONVERTED TO MODIFIED FILTERED FEATURE TABLE IF NEEDED
  
  # Apply a function across the rows to check for any non-zero values in the site's columns
  # and sum these to count the number of non-zero rows for the site
  counters[site] <- sum(apply(filtered_feature_table[site_cols], 1, function(row) any(row != 0)))
  
  # Subset the data for rows where any column for the site has a non-zero value
  site_data <- filtered_feature_table[apply(filtered_feature_table[site_cols], 1, function(row) any(row != 0)), ]
  
  # Store the subseted data frame in the list with a name
  site_data_frames[[paste0("found_site_", site)]] <- site_data
}

#BEGINNING OF SUBSECTION FOR SITE NUMBER AND SAMPLETIME
get_columns_for_site_and_time <- function(mapfile, site_number, sample_time) {
  if(site_number < 10) {
    site_number <- paste0("0", site_number)
  }
  file_names <- mapfile$File.Name[mapfile$Sample_site == site_number & mapfile$Sample_time == sample_time]
  return(file_names)
}

#  Find the chemicals and number of chemicals that were found at each site from 1 to 28 
site_and_time_data_frames <- list()

for(site in 1:28) {
  if (site == 20) {
    time_points <- c("A", "B")
  } 
  else {
    time_points <- c("A", "B", "C")
  }
  for(time in time_points) {
      # Get the file names for the site and timepoint using the existing function
      site_files <- get_columns_for_site_and_time(mapfile2, site, time)
      
      # Determine the column names in filtered_feature_table that match the site files
      site_cols <- which(colnames(filtered_feature_table) %in% site_files)
      
      # Subset the data for rows where any column for the site has a non-zero value
      site_data <- filtered_feature_table[apply(filtered_feature_table[site_cols], 1, function(row) any(row != 0)), ]
      
      # Store the subseted data frame in the list with a name
      site_and_time_data_frames[[paste0("found_site_", site, "_and_time_", time)]] <- site_data
    }
}
#END OF SUBSECTION FOR SITE NUMBER AND SAMPLETIME

#BEGINNING OF SUBSECTION LINEAR CORRELATION BETWEEN SITES 
percent_match <- function(data, filename1, filename2, colname) {
  data_current <- data[[filename1]]
  data_compare <- data[[filename2]]
  
  id_current <- data_current[[colname]]
  id_compare <- data_compare[[colname]] 
  match <- sum(id_current %in% id_compare)
  percentage <- (match/length(id_compare))*100
  return(percentage)
}


sites_15 <- c(16)
results_for_site_15 <- list() 
percentage_match_site_15 <- list()
for (site_number in sites_15) {
  for (time in c("A", "B", "C")) {
    site_str <- as.character(site_number)
    current_file_name <- mapfile2 %>%
      filter(Sample_site == site_str, Sample_time == time) %>%
      pull(File.Name)
    
    file_name_15 <- mapfile2 %>%
      filter(Sample_site == "15", Sample_time == time) %>%
      pull(File.Name)
    
    r_value <- cor(filtered_feature_table[[current_file_name]], filtered_feature_table[[file_name_15]])
    results_for_site_15[[paste("Site", site_number, "vs 15 at time", time)]] <- r_value
    
   percentage_key <- paste("Site", site_number, "vs 15 at time", time, "% Match")
   percentage_match_site_15[[percentage_key]] <- percent_match(site_and_time_data_frames, paste0("found_site_", site_number, "_and_time_", time), paste0("found_site_15_and_time_", time), "chemical_ID")
  }
}

sites_21 <- c(18, 20, 22:24)
results_for_site_21 <- list() 
percentage_match_site_21 <- list()
for (site_number in sites_21) {
  if (site_number == 20) {
    time_points <- c("A", "B")
  } else {
    time_points <- c("A", "B", "C")
  }
  site_str <- as.character(site_number)
  for (time in time_points) {
    current_file_name <- mapfile2 %>%
    filter(Sample_site == site_str, Sample_time == time) %>%
    pull(File.Name) 
      file_name_21 <- mapfile2 %>%
        filter(Sample_site == "21", Sample_time == time) %>%
        pull(File.Name)
      r_value <- cor(filtered_feature_table[[current_file_name]], filtered_feature_table[[file_name_21]])
      results_for_site_21[[paste("Site", site, "vs. 21 at time", time)]] <- r_value
      
      percentage_key <- paste("Site", site_number, "vs 21 at time", time, "% Match")
      percentage_match_site_21[[percentage_key]] <- percent_match(site_and_time_data_frames, paste0("found_site_", site_number, "_and_time_", time), paste0("found_site_21_and_time_", time), "chemical_ID")
    }
  }


sites_27 <- c(1:5, 7:14, 25, 26, 28)
percentage_match_site_27 <- list()
results_for_site_27 <- list() 
for (site_number in sites_27) {
  if (site_number < 10) {
    site_str <- paste0("0", site_number)
  }
  else {
    site_str <- as.character(site_number)
  }
  for (time in c("A", "B", "C")) {
    current_file_name <- mapfile2 %>%
      filter(Sample_site == site_str, Sample_time == time) %>%
      pull(File.Name)
    
    file_name_27 <- mapfile2 %>%
      filter(Sample_site == "27", Sample_time == time) %>%
      pull(File.Name)
    
    r_value <- cor(filtered_feature_table[[current_file_name]], filtered_feature_table[[file_name_27]])
    results_for_site_27[[paste("Site", site_number, "vs 27 at time", time)]] <- r_value
    
    percentage_key <- paste("Site", site_number, "vs 27 at time", time, "% Match")
    percentage_match_site_27[[percentage_key]] <- percent_match(site_and_time_data_frames, paste0("found_site_", site_number, "_and_time_", time), paste0("found_site_27_and_time_", time), "chemical_ID")
  }
}

#END OF SUBSECTION LIENAR CORRELATION BETWEEN SITES

#BEGINNING OF SUBSECTION CREATING BAR GRAPHS COMPARING SITES AND WQTC SITES
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
  "16", "15",                                                             
  "18", "20", as.character(22:24), "21", "6",
  as.character(1:5), as.character(7:14), as.character(25:26), "28", "27"
)

# Set the levels of the Site factor according to the custom order
counters_df$Site <- factor(counters_df$Site, levels = site_order)

counters_df$Group <- ifelse(counters_df$Site %in% c(1:5, 7:14, 25:26, 28, 27), "Site 27",
                            ifelse(counters_df$Site %in% c(15,16), "Site 15",
                                   ifelse(counters_df$Site %in% c(18, 20, 22:24, 21), "Site 21", NA)))

# Define a new column for color (factor variable for differentiating specific bars)
counters_df$Color <- ifelse(counters_df$Site %in% c(15, 21, 27), "WQTC", "Nested")

# Plotting the bar graph
ggplot(counters_df, aes(x = Group, y = Count, fill = Color, group = Site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.45), width = 0.37) +
  theme_minimal() +
  labs(title = "Number of Chemicals Detected At Sites Compared to Their Corresponding WQTC", x = "Group", y = "Count") +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Nested" = "blue", "WQTC" = "orange")) +
  scale_x_discrete(limits = c("Site 15", "Site 21", "Site 27"))

#END OF SUBSECTION CREATING BAR GRAPHS COMPARING SITES AND WQTC SITES
# END of SECTION 1 


# BEGINNING of SECTION 2

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

# END of SECTION 2


#BEGINNING of SECTION 3
#SUBSECTION 1 of SECTION 3: SCATTERPLOTS FOR INCOME, POPULATION, AND AREA
incomes_total <- c(0,34490, 36812, 37059, 40168, 31050, 27752, 27517, 65791, 56672, 27517, 74500, 101140, 86478, 81994, 108021, 55433, 76796, 72401, 55346, 61923, 45794, 53857, 61081, 28054, 52751, 20000)
population_total <- c(0, 8258, 3459, 1610, 3587, 10949, 9073, 23751, 145346, 8838, 23751, 95603, 11444, 40824, 5781, 37193, 78206, 60885, 24969, 309998, 45148, 37972, 23135, 309184, 41777, 350766, 74)
area_total <- c(0.04, 4, 1, 1, 3, 5, 5, 12, 112, 3, 12, 80, 12, 67, 11, 88, 55, 80, 23, 332, 37, 28, 21, 242, 20, 280, 3)

incomes <- c(0,34490, 36812, 37059, 40168, 31050, 27752, 27517, 65791, 56672, 27517, 74500, 101140, 81994, 55433, 72401, 61923, 45794, 53857, 61081, 28054, 20000)
population <- c(0, 8258, 3459, 1610, 3587, 10949, 9073, 23751, 145346, 8838, 23751, 95603, 11444, 5781, 78206, 24969, 45148, 37972, 23135, 309184, 41777, 74)
area <- c(0.04, 4, 1, 1, 3, 5, 5, 12, 112, 3, 12, 80, 12, 11, 55, 23, 37, 28, 21, 242, 20, 3)

#incomes excludes WQTC's, incomes_total WQTC's

#creating the data table from the paper provided by Dr. Walker
mod_counters_df <- counters_df[-c(6,15,17,19,21,27), ] #Remove site 6 because it was a blank in the paper provided by Dr. Walker
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

#SUBSECTION 2 OF SECTION 3: BOXPLOTS

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

#END of SECTION 3





