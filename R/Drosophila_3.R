#ATTRIBUTIONS
## Original Script & Project Creation: Sabrina
## Edits and modifications: Nadira

# Load in necessary libraries
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer("dplyr","select")
library(viridis)
theme_set(theme_light())
library("vegan")
library(MASS)
library(ggplot2)
### EDIT: #install.packages("ggVennDiagram")
library(ggVennDiagram)

# Clear the workspace

rm(list = ls())

# Load and filter raw data ----

dros_data <- read_tsv("../data/drosophila_data.tsv")
names(dros_data)

# Clean the data to remove N/A and duplicate samples entries in the columns bin_uri & country

df_usa <- dros_data %>%
  distinct(processid, .keep_all = TRUE) %>%
  filter(!is.na(`country/ocean`) &
    !is.na(`province/state`) &
    !is.na(bin_uri) &
    str_detect(`country/ocean`, regex("united states", ignore_case = TRUE)))

# Split USA into Mainland and Hawaii ----

df_hawaii_na <- df_usa %>%
  filter(str_detect(`province/state`, regex("hawaii", ignore_case = TRUE)))

df_mainland <- df_usa %>%
  filter(str_to_lower(`province/state`) != "hawaii")

# Identify regions of Hawaii
table(df_hawaii_na$region)

# Assign island geological ages based on region

df_hawaii <- df_hawaii_na %>%
  mutate(
    island_age = case_when(
      region %in% c("Captian Cook", "Hawaii", "Hilo", "Kamuela", "Maunakea") ~ 0.5,
      region %in% c("Kauai") ~ 5.1,
      region %in% c("Honolulu", "O'ahu", "Oahu") ~ 3.0,
      region %in% c("Maui", "Lanai") ~ 1.3,
      region %in% c("Molokai") ~ 1.9
    )
  )

# Identify any NAs

sum(is.na(df_hawaii$island_age))

# There are 8, I need to figure out from where

df_test <- df_hawaii %>%
  filter(is.na(island_age)) %>%
  select(`country/ocean`, `province/state`, region, island_age)

# They are from the missing region information, I will convert these NAs to Unknowns and use this information for the Mainland vs Hawaii, but filter out these entries for the Hawaii Island Analysis

df_hawaii <- df_hawaii %>%
  mutate(region = if_else(is.na(region), "Unknown", region))

# Now this output should be 0
sum(is.na(df_hawaii$region))

# Remove unnecessary variables from workspace
rm(df_hawaii_na, df_test)

# EDIT #1: Set custom theme to avoid repeating theme for plots ----

custom_theme <- theme_light() +
  theme(
    plot.background = element_rect(fill = "snow2", colour = "snow2"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(color = "gray90")
  )
theme_set(custom_theme)

# EDIT #2: Replacement Venn Diagram for BIN Overlap between Hawaii and Mainland USA ----

# Get Unique Bins for Hawaii and Mainland

hawaii_bins <- unique(df_hawaii$bin_uri)
mainland_bins <- unique(df_mainland$bin_uri)

# Create a list for the unique bins to graph

bin_sets <- list(
  Hawaii = hawaii_bins,
  Mainland = mainland_bins
)

# Use Venn Diagram package to plot

ggVennDiagram(bin_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Venn Diagram: BIN Overlap Between Hawaii and Mainland USA")

# Accumulation Curve to demonstrate the bin richness in the varying areas ----

df_H <- df_hawaii %>%
  count(sampleid, bin_uri, name = "n")
df_H <- pivot_wider(df_H, names_from = bin_uri, values_from = n, values_fill = list(n = 0))
df_H <- df_H %>%
  column_to_rownames("sampleid") %>%
  as.matrix()

df_M <- df_mainland %>%
  count(sampleid, bin_uri, name = "n")
df_M <- pivot_wider(df_M, names_from = bin_uri, values_from = n, values_fill = list(n = 0))
df_M <- df_M %>%
  column_to_rownames("sampleid") %>%
  as.matrix()

df_all <- df_usa %>%
  count(sampleid, bin_uri, name = "n")
df_all <- pivot_wider(df_all, names_from = bin_uri, values_from = n, values_fill = list(n = 0))
df_all <- df_all %>%
  column_to_rownames("sampleid") %>%
  as.matrix()

# Run specaccum for each region

acc_H <- specaccum(df_H, method = "random", permutations = 1000)
acc_M <- specaccum(df_M, method = "random", permutations = 1000)
acc_ALL <- specaccum(df_all, method = "random", permutations = 1000)

# Convert specaccum objects to data frames

df_H <- data.frame(
  samples = acc_H$sites,
  richness = acc_H$richness,
  group = "Hawaii"
)

df_M <- data.frame(
  samples = acc_M$sites,
  richness = acc_M$richness,
  group = "Mainland USA"
)

df_ALL <- data.frame(
  samples = acc_ALL$sites,
  richness = acc_ALL$richness,
  group = "All (USA+Hawaii)"
)

# Bind for plotting
df_plot_hm <- bind_rows(df_H, df_M)
df_plot_all <- bind_rows(df_H, df_M, df_ALL)

# Simple gglpot for Hawaii and Mainland
ggplot(df_plot_hm, aes(x = samples, y = richness, color = group, linetype = group)) +
  geom_line(linewidth = 1) +
  labs(
    title = "BIN Accumulation: Hawaii vs Mainland USA",
    x = "Number of Samples", y = "BIN Richness", color = "BIN Region",
    linetype = "BIN Region"
  ) +
  scale_color_manual(
    values = c("Hawaii" = "blue", "Mainland USA" = "orange2")
  )

ggsave(filename = "../figs/hm_bin_accumulation.png", width = 6, height = 4, dpi = 300)

# Simple ggplot for ALL
ggplot(df_plot_all, aes(x = samples, y = richness, color = group, linetype = group)) + # box graph color = fill
  labs(
    title = "BIN Accumulation: Hawaii vs Mainland USA vs USA",
    x = "Number of Samples", y = "BIN Richness", color = "BIN Region",
    linetype = "BIN Region"
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Hawaii" = "blue", "Mainland USA" = "orange2", "All (USA+Hawaii)" = "purple4")
  )

ggsave(filename = "../figs/hmusa_bin_accumulation.png", width = 6, height = 4, dpi = 300)

# Hawaii Age Island Comparison of BIN richness ----

dfh_island <- df_hawaii[df_hawaii$region != "Unknown", ]

# Filter the data further by finding the unique BINs per Hawaii regions
dfh_bin_regions <- dfh_island %>%
  group_by(island_age, bin_uri) %>%
  distinct(bin_uri) %>%
  count() %>%
  group_by(island_age) %>%
  summarise(bins = sum(n))

ggplot(dfh_bin_regions, aes(x = factor(island_age), y = bins, fill = island_age)) +
  geom_col() +
  geom_text(aes(label = bins), vjust = -0.5, size = 4) +
  scale_fill_gradient(low = "blue", high = "orange2") +
  ylim(0, 60) +
  labs(
    title = "Distinct BINs by Island Age",
    x = "Island Age (millions of years)",
    y = "Number of Distinct BINs",
    fill = "Island Age"
  )

ggsave(filename = "../figs/h_age_bar.png", width = 6, height = 4, dpi = 300)

# EDIT #3: Replacement loop for checking for overlap between islands ----

# Make vector for unique ages of island

ages <- c(0.5, 1.3, 1.9, 3.0, 5.1)

# Use a loop to fill the islands list

islands <- list()

for (age in ages) {
  subset <- dfh_island[dfh_island$island_age == age, ]
  islands[[as.character(age)]] <- subset$bin_uri
}

# use a similarity matrix to show how many BINs overlap between islands

overlap_matrix <- sapply(islands, function(x) sapply(islands, function(y) length(intersect(x, y))))

overlap_matrix

capture.output(overlap_matrix, file = "../output/overlap_matrix.txt")

# it has been shown that only a couple islands share one BIN, demonstrating the unique diversity of the Hawaii Islands

# it is evident that our data isn't normally distributed, and therefore we should calcuate the mean and variance to determine which statistical test is appropriate
# comparing mean and variance to determine regression choice

mean_bins <- mean(dfh_bin_regions$bins)
var_bins <- var(dfh_bin_regions$bins)

mean_bins
var_bins

# the variance is drastically larger than the mean and therefore the negative bionomial regression is the most appropriate choice for this count data

# Negative Bionomial Regression chosen for count data and overdispersion shown above ----

nb_model <- glm.nb(bins ~ island_age, data = dfh_bin_regions)
summary(nb_model)

capture.output(summary(nb_model), file = "../output/nb_model_summary.txt")

# Plot regression result

ggplot(dfh_bin_regions, aes(x = island_age, y = bins)) +
  geom_point(size = 3, color = "orange3") +
  stat_smooth(
    method = MASS::glm.nb,
    formula = y ~ x,
    color = "darkblue",
    se = TRUE
  ) +
  labs(
    title = "Negative Binomial Regression: Island Age vs BIN Diversity",
    x = "Island Age (millions of years)",
    y = "Distinct BIN Count"
  )

ggsave(filename = "../figs/h_age_nbr_bin_regression.png", width = 6, height = 4, dpi = 300)
