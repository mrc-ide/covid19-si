library(dplyr)
relaxed <- readRDS("processed_stanfits/relaxed_priors/nf_overall_table2.rds")
s3s4mix <- readRDS("processed_stanfits/s3s4mix/nf_overall_table2.rds")
discrete <- readRDS("processed_stanfits/discrete_pairs/nf_overall_table2.rds")
x <- list(
  relaxed = relaxed, s3s4mix = s3s4mix,
  discrete_pairs = discrete
)
x <- bind_rows(x, .id = "meta-model")

## same underlying data
same_data <- list(
  relaxed = relaxed, s3s4mix = s3s4mix
) %>% bind_rows(.id = "meta-model")
same_data <- arrange(same_data, DIC)
